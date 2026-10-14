#!/usr/bin/env bash
# Submit and inspect the local personal Jenkins Pipeline. Jenkins itself must
# already be running (for example as a Homebrew or systemd service).
set -euo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
JENKINS_URL="${JENKINS_URL:-http://127.0.0.1:8080}"
JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-personal}"
JENKINS_CLI_AUTH_FILE="${JENKINS_CLI_AUTH_FILE:-$HOME/.config/gkeyll/jenkins/personal.auth}"
JENKINS_CLI_JAR="${JENKINS_CLI_JAR:-${TMPDIR:-/tmp}/gkeyll-jenkins-cli.jar}"
CURL_CONFIG=''
QUEUE_ID=''
RESOLVED_BUILD_NUMBER=''

source "$SCRIPT_DIR/jenkins-client-common.sh"
source "$SCRIPT_DIR/jenkins-client-artifacts.sh"

die() { echo "ERROR: $*" >&2; exit 1; }
usage() {
    cat <<'EOF'
Usage: jenkins-personal.sh <command> [flags]

Commands:
  run      Queue a pull-request or comparison build.
  active   List this job's queued and running work.
  recent   List retained builds.
  info     Show detailed information for a retained build.
  artifact List or download retained build artifacts.
  status   Show a queued or known build's state.
  follow   Wait for and stream the console of a queued or known build.
  abort    Cancel a queued or running build.

The run command returns after Jenkins accepts the request. --follow streams
the build console and returns its final Jenkins result. Press Ctrl-C to stop
following without aborting the Jenkins build.

Jenkins must already be running. Credentials default to
~/.config/gkeyll/jenkins/personal.auth, a protected file containing one line:
jenkins-user:api-token. JENKINS_CLI_AUTH_FILE, JENKINS_URL, and JENKINS_JOB
override their defaults.
EOF
}
command_usage() {
    case "$1" in
        run) cat <<'EOF'
Usage: jenkins-personal.sh run (--pr NUMBER | --candidate-ref REF --baseline-ref REF) [--ci-ref REF] [--allow-behind-candidate] [--follow]

Flags:
  --pr NUMBER           Build GitHub pull request NUMBER.
  --candidate-ref REF   Candidate branch or commit; requires --baseline-ref.
  --baseline-ref REF    Baseline branch or commit; requires --candidate-ref.
  --ci-ref REF          Jenkinsfile and CI helpers branch or full commit SHA.
  --allow-behind-candidate  Permit a candidate that does not contain the baseline.
  --follow              Stream the build console after Jenkins queues it.
EOF
        ;;
        active) cat <<'EOF'
Usage: jenkins-personal.sh active

This command takes no options.
EOF
        ;;
        recent) cat <<'EOF'
Usage: jenkins-personal.sh recent [--limit NUMBER]

Flags:
  --limit NUMBER        Number of retained builds to list (default: 10).
EOF
        ;;
        status) cat <<'EOF'
Usage: jenkins-personal.sh status (--queue ID | --build NUMBER)

Flags:
  --queue ID            Inspect a Jenkins queue item.
  --build NUMBER        Inspect a known Jenkins build.
EOF
        ;;
        info) cat <<'EOF'
Usage: jenkins-personal.sh info --build NUMBER

Flags:
  --build NUMBER        Show detailed information for a retained Jenkins build.
EOF
        ;;
        artifact) cat <<'EOF'
Usage: jenkins-personal.sh artifact --build NUMBER [--list | --fetch [--only PATH[,PATH...]] [--output-dir DIR]]

Flags:
  --build NUMBER        Select a retained Jenkins build.
  --list                List artifacts (the default).
  --fetch               Download artifacts.
  --only PATHS          Comma-separated artifact-relative paths to download.
  --output-dir DIR      New directory for downloaded artifacts.
EOF
        ;;
        follow) cat <<'EOF'
Usage: jenkins-personal.sh follow (--queue ID | --build NUMBER)

Flags:
  --queue ID            Wait for this queue item, then stream its build.
  --build NUMBER        Stream this known Jenkins build.
EOF
        ;;
        abort) cat <<'EOF'
Usage: jenkins-personal.sh abort (--queue ID | --build NUMBER)

Flags:
  --queue ID            Cancel this Jenkins queue item.
  --build NUMBER        Abort this running Jenkins build.
EOF
        ;;
    esac
}
java_command() {
    local candidate
    if [[ -n "${JAVA_HOME:-}" && -x "$JAVA_HOME/bin/java" ]]; then
        printf '%s' "$JAVA_HOME/bin/java"
    elif candidate="$(/usr/libexec/java_home 2>/dev/null)" && [[ -x "$candidate/bin/java" ]]; then
        printf '%s' "$candidate/bin/java"
    elif candidate="$(ls -d /opt/homebrew/opt/openjdk*/bin/java /usr/local/opt/openjdk*/bin/java 2>/dev/null | sort -V | tail -n 1)" && [[ -x "$candidate" ]]; then
        # Homebrew's JDK (the one Jenkins itself runs on); macOS's /usr/bin/java is only a stub.
        printf '%s' "$candidate"
    else
        command -v java || die 'Set JAVA_HOME or put Java 21 or newer on PATH to stream the Jenkins console'
    fi
}
submit() {
    local pr="$1" candidate="$2" baseline="$3" allow_behind="$4" ci_ref="$5" headers queue
    headers="$(mktemp "${TMPDIR:-/tmp}/gkeyll-jenkins-headers.XXXXXX")"
    curl_auth --dump-header "$headers" --output /dev/null --request POST \
        --data-urlencode "CANDIDATE_PR=$pr" --data-urlencode "CANDIDATE_REF=$candidate" --data-urlencode "BASELINE_REF=$baseline" --data-urlencode "ALLOW_BEHIND_CANDIDATE=$allow_behind" --data-urlencode "CI_REF=$ci_ref" \
        "${JENKINS_URL}$(job_path)/buildWithParameters" || { rm -f "$headers"; die 'Jenkins rejected the build'; }
    queue="$(awk 'BEGIN{IGNORECASE=1} /^Location:/{sub(/^[^:]*: /,""); sub(/\r$/,""); print; exit}' "$headers")"; rm -f "$headers"
    [[ "$queue" =~ /queue/item/([0-9]+)/ ]] || die 'Jenkins accepted the build but returned no queue ID'
    QUEUE_ID="${BASH_REMATCH[1]}"
    echo "Queued $JENKINS_JOB as queue item $QUEUE_ID"
}
run() {
    local pr='' candidate='' baseline='' ci_ref='' allow_behind=false want_follow=false
    while (($#)); do case "$1" in --pr) (($#>=2))||die '--pr requires a number'; pr="$2"; shift 2;; --candidate-ref) (($#>=2))||die '--candidate-ref requires a ref'; candidate="$2"; shift 2;; --baseline-ref) (($#>=2))||die '--baseline-ref requires a ref'; baseline="$2"; shift 2;; --ci-ref) (($#>=2)) && [[ -n "$2" ]] || die '--ci-ref requires a ref'; ci_ref="$2"; shift 2;; --allow-behind-candidate) allow_behind=true; shift;; --follow) want_follow=true; shift;; *) die "unknown run option: $1";; esac; done
    if [[ -n "$pr" ]]; then positive '--pr' "$pr"; [[ -z "$candidate$baseline" ]] || die '--pr cannot be combined with refs'; else [[ -n "$candidate" && -n "$baseline" ]] || die 'provide --pr, or both --candidate-ref and --baseline-ref'; fi
    if [[ -n "$ci_ref" ]]; then
        local definitions
        definitions="$(curl_auth "${JENKINS_URL}$(job_path)/api/json?tree=property[parameterDefinitions[name]]")"
        python3 -c 'import json,sys; sys.exit(not any(p.get("name") == "CI_REF" for prop in json.load(sys.stdin).get("property", []) for p in prop.get("parameterDefinitions", [])))' <<< "$definitions" \
            || die 'This job does not expose CI_REF yet; run the updated personal Jenkinsfile once to register it.'
    fi
    submit "$pr" "$candidate" "$baseline" "$allow_behind" "$ci_ref"
    if [[ "$want_follow" == true ]]; then
        wait_for_build_number "$QUEUE_ID"
        follow "$RESOLVED_BUILD_NUMBER"
    fi
}
format_builds() {
    local limit="$1" active_only="$2" payload="$3"
    python3 -c '
import datetime,json,sys
limit=int(sys.argv[1]); active_only=sys.argv[2] == "true"; count=0
def selectors(actions):
    values={}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            if parameter.get("name") in ("CANDIDATE_PR","CANDIDATE_REF","BASELINE_REF","CI_REF"):
                values[parameter["name"]]=parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2} ci={3}".format(values.get("CANDIDATE_PR","-"),values.get("CANDIDATE_REF","-"),values.get("BASELINE_REF","-"),values.get("CI_REF","default"))
for build in json.load(sys.stdin).get("builds",[]):
    if active_only and not build.get("building",False): continue
    if count >= limit: break
    timestamp=build.get("timestamp")
    when="-" if timestamp is None else datetime.datetime.fromtimestamp(timestamp/1000.0,datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    state="RUNNING" if build.get("building",False) else (build.get("result") or "UNKNOWN")
    print("BUILD #{number} {state} queue={queue} {when} {selectors} url={url}".format(number=build.get("number","?"),state=state,queue=build.get("queueId") if build.get("queueId") is not None else "-",when=when,selectors=selectors(build.get("actions")),url=build.get("url","-")))
    count+=1
' "$limit" "$active_only" <<< "$payload"
}
format_queue_items() {
    local payload="$1"
    python3 -c '
import json,sys
job=sys.argv[1]
def selectors(actions):
    values={}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            if parameter.get("name") in ("CANDIDATE_PR","CANDIDATE_REF","BASELINE_REF","CI_REF"):
                values[parameter["name"]]=parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2} ci={3}".format(values.get("CANDIDATE_PR","-"),values.get("CANDIDATE_REF","-"),values.get("BASELINE_REF","-"),values.get("CI_REF","default"))
for item in json.load(sys.stdin).get("items",[]):
    if (item.get("task") or {}).get("name") != job or item.get("executable"): continue
    state="CANCELLED" if item.get("cancelled",False) else "QUEUED"
    print("QUEUE {id} {state} reason={why} {selectors}".format(id=item.get("id","?"),state=state,why=" ".join((item.get("why") or "waiting").split()),selectors=selectors(item.get("actions"))))
' "$JENKINS_JOB" <<< "$payload"
}
main() {
    (($#)) || { usage; exit 2; }
    case "$1" in -h|--help|help) usage; return;; esac
    if (($# == 2)) && [[ "$2" == -h || "$2" == --help ]]; then
        case "$1" in run|active|recent|info|artifact|status|follow|abort) command_usage "$1"; return;; esac
    fi
    prepare_auth
    case "$1" in run) shift; run "$@";; follow) shift; follow_command "$@";; status) shift; status_command "$@";; abort) shift; abort "$@";; active) shift; [[ $# == 0 ]] || die 'usage: active'; active;; recent) shift; [[ $# == 0 || ( $# == 2 && $1 == --limit ) ]] || die 'usage: recent [--limit NUMBER]'; recent "${2:-10}";; info) shift; ci_info_command "$@";; artifact) shift; ci_artifact_command "$@";; -h|--help|help) usage;; *) usage >&2; die "unknown command: $1";; esac
}
main "$@"
