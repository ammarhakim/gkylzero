#!/usr/bin/env bash
# Shared Jenkins HTTP client. Source after setting JENKINS_URL, JENKINS_JOB,
# JENKINS_CLI_AUTH_FILE and JENKINS_CLI_JAR. Callers provide die(),
# java_command(), format_builds() and format_queue_items(); sourcing has no
# network or controller side effects. prepare_auth installs credential cleanup.

path_for() {
    local job="$1" path='' part
    local -a parts
    IFS=/ read -r -a parts <<< "$job"
    for part in "${parts[@]}"; do path+="/job/$part"; done
    printf '%s' "$path"
}
job_path() { path_for "$JENKINS_JOB"; }

ci_build_url() { printf '%s%s/%s' "$JENKINS_URL" "$(job_path)" "$1"; }

prepare_auth() {
    # Cluster queue-to-build status dispatch may prepare authentication twice.
    [[ -n "${CURL_CONFIG:-}" && -f "$CURL_CONFIG" ]] && return 0
    [[ -f "$JENKINS_CLI_AUTH_FILE" ]] || die "credential file is missing: $JENKINS_CLI_AUTH_FILE"
    [[ -O "$JENKINS_CLI_AUTH_FILE" ]] || die "credential file is not owned by $USER"
    [[ "$(stat -c '%a' "$JENKINS_CLI_AUTH_FILE" 2>/dev/null || stat -f '%Lp' "$JENKINS_CLI_AUTH_FILE")" == 600 ]] || die 'credential file must have mode 600'
    local credential; credential="$(<"$JENKINS_CLI_AUTH_FILE")"
    [[ "$credential" =~ ^[^[:space:]:]+:[^[:space:]:]+$ ]] || die 'credential file must contain user:api-token'
    [[ "$credential" != *'"'* && "$credential" != *'\'* ]] || die 'credential file contains an unsupported character'
    CURL_CONFIG="$(mktemp "${TMPDIR:-/tmp}/gkeyll-jenkins.XXXXXX")"; chmod 600 "$CURL_CONFIG"
    printf 'user = "%s"\n' "$credential" > "$CURL_CONFIG"
    trap 'rm -f "$CURL_CONFIG"' EXIT
}

curl_auth() { curl --fail --silent --show-error --globoff --config "$CURL_CONFIG" "$@"; }

curl_auth_quiet() { curl --fail --silent --globoff --config "$CURL_CONFIG" "$@"; }

positive() { [[ "$2" =~ ^[1-9][0-9]*$ ]] || die "$1 must be a positive integer"; }

download_cli() {
    if [[ ! -s "$JENKINS_CLI_JAR" ]]; then
        local directory temporary
        directory="$(dirname "$JENKINS_CLI_JAR")"
        mkdir -p "$directory"
        temporary="$(mktemp "$directory/gkeyll-jenkins-cli.XXXXXX")"
        curl --fail --silent --show-error --output "$temporary" "$JENKINS_URL/jnlpJars/jenkins-cli.jar"
        mv "$temporary" "$JENKINS_CLI_JAR"
    fi
}

state() {
    local payload
    payload="$(curl_auth "${JENKINS_URL}$(job_path)/$1/api/json")" || return 75
    python3 -c 'import json,sys; x=json.load(sys.stdin); print("BUILDING" if x["building"] else (x.get("result") or "UNKNOWN"))' <<< "$payload"
}

queue_state() {
    local queue_id="$1" payload
    payload="$(curl_auth_quiet "${JENKINS_URL}/queue/item/$queue_id/api/json")" || return 1
    python3 -c '
import json,sys
item=json.load(sys.stdin); executable=item.get("executable") or {}
print("{}|{}|{}".format(executable.get("number", ""), "true" if item.get("cancelled", False) else "false", " ".join((item.get("why") or "").split())))
' <<< "$payload"
}

build_for_queue() {
    local queue_id="$1" payload
    payload="$(curl_auth "${JENKINS_URL}$(job_path)/api/json?tree=builds[number,queueId]")" || return 1
    python3 -c '
import json,sys
queue_id=int(sys.argv[1])
for build in json.load(sys.stdin).get("builds",[]):
    if build.get("queueId") == queue_id:
        print(build["number"]); break
' "$queue_id" <<< "$payload"
}

wait_for_build_number() {
    local queue_id="$1" record number cancelled why previous='' attempt
    while :; do
        if ! record="$(queue_state "$queue_id")"; then
            for attempt in {1..6}; do
                if number="$(build_for_queue "$queue_id")" && [[ -n "$number" ]]; then
                    RESOLVED_BUILD_NUMBER="$number"
                    return 0
                fi
                sleep 2
            done
            die "Queue item $queue_id is unavailable and no matching Jenkins build was found"
        fi
        IFS='|' read -r number cancelled why <<< "$record"
        # Jenkins can report cancelled=true briefly while the item becomes a build.
        if [[ "$cancelled" == true ]]; then for attempt in {1..12}; do if number="$(build_for_queue "$queue_id")" && [[ "$number" =~ ^[1-9][0-9]*$ ]]; then RESOLVED_BUILD_NUMBER="$number"; return 0; fi; sleep 5; done; die "Queue item $queue_id was cancelled"; fi
        if [[ -n "$number" ]]; then
            RESOLVED_BUILD_NUMBER="$number"
            return 0
        fi
        if [[ "$why" != "$previous" ]]; then
            echo "Queue item $queue_id: ${why:-waiting}" >&2
            previous="$why"
        fi
        sleep 5
    done
}

status_build() {
    positive 'build number' "$1"
    local s
    s="$(state "$1")" || die "Build #$1 is not available yet. Retry shortly."
    echo "Build #$1: $s"
    [[ "$s" == SUCCESS ]]
}

status_command() {
    [[ $# == 2 ]] || die 'usage: status --queue ID | status --build NUMBER'
    case "$1" in
        --queue)
            positive 'queue ID' "$2"
            local record number cancelled why
            if ! record="$(queue_state "$2")"; then
                number="$(build_for_queue "$2")" || die "Queue item $2 is unavailable"
                [[ -n "$number" ]] || die "Queue item $2 is unavailable"
                echo "Queue item $2 is build #$number"
                status_build "$number"
                return
            fi
            IFS='|' read -r number cancelled why <<< "$record"
            if [[ -n "$number" ]]; then
                echo "Queue item $2 is build #$number"
                status_build "$number"
            elif [[ "$cancelled" == true ]]; then
                echo "Queue item $2: CANCELLED"
                return 1
            else
                echo "Queue item $2: ${why:-waiting}"
            fi
            ;;
        --build) status_build "$2" ;;
        *) die 'usage: status --queue ID | status --build NUMBER' ;;
    esac
}

abort() {
    [[ $# == 2 ]] || die 'usage: abort --queue ID | abort --build NUMBER'
    case "$1" in
        --queue)
            positive 'queue ID' "$2"
            curl_auth --output /dev/null --request POST --data-urlencode "id=$2" \
                "${JENKINS_URL}/queue/cancelItem"
            echo "Requested cancellation of queue item $2"
            ;;
        --build)
            positive 'build number' "$2"
            curl_auth --output /dev/null --request POST \
                "${JENKINS_URL}$(job_path)/$2/stop"
            echo "Requested cancellation of $JENKINS_JOB #$2"
            ;;
        *) die 'usage: abort --queue ID | abort --build NUMBER' ;;
    esac
}

follow() {
    positive 'build number' "$1"
    local java
    java="$(java_command)"
    download_cli
    echo "Following $JENKINS_JOB #$1"
    echo "Build URL: ${JENKINS_URL}$(job_path)/$1/"
    trap "echo 'Stopped following build #$1'; exit 130" INT TERM HUP
    if ! "$java" -jar "$JENKINS_CLI_JAR" -s "$JENKINS_URL" -http \
        -auth "@$JENKINS_CLI_AUTH_FILE" console "$JENKINS_JOB" "$1" -f; then
        echo 'Console follower ended before Jenkins reported a terminal result; checking build status.' >&2
    fi
    trap - INT TERM HUP
    while :; do
        local s
        if ! s="$(state "$1")"; then
            echo "Build #$1 is not available yet; retrying." >&2
            sleep 2
            continue
        fi
        [[ "$s" == BUILDING ]] || { echo "Build #$1: $s"; [[ "$s" == SUCCESS ]]; return; }
        sleep 5
    done
}

follow_command() {
    [[ $# == 2 ]] || die 'usage: follow --queue ID | follow --build NUMBER'
    case "$1" in
        --queue) positive 'queue ID' "$2"; wait_for_build_number "$2"; [[ "$RESOLVED_BUILD_NUMBER" =~ ^[1-9][0-9]*$ ]] || die "Queue item $2 did not resolve to a build number (got '$RESOLVED_BUILD_NUMBER')"; follow "$RESOLVED_BUILD_NUMBER" ;;
        --build) follow "$2" ;;
        *) die 'usage: follow --queue ID | follow --build NUMBER' ;;
    esac
}

job_builds() { curl_auth "${JENKINS_URL}$(job_path)/api/json?tree=builds[number,queueId,building,result,timestamp,url,actions[parameters[name,value]]]"; }

queue_items() { curl_auth "${JENKINS_URL}/queue/api/json?tree=items[id,task[name],why,cancelled,executable[number],actions[parameters[name,value]]]"; }

active() {
    local builds queues output
    builds="$(job_builds)"; queues="$(queue_items)"
    output="$(format_builds 100 true "$builds")"
    output+=$'\n'"$(format_queue_items "$queues")"
    [[ -n "${output//$'\n'/}" ]] && printf '%s\n' "$output" || echo "No queued or running builds for $JENKINS_JOB"
}

recent() {
    local limit="${1:-10}" output
    positive '--limit' "$limit"
    output="$(format_builds "$limit" false "$(job_builds)")"
    [[ -n "$output" ]] && printf '%s\n' "$output" || echo "No retained builds for $JENKINS_JOB"
}
