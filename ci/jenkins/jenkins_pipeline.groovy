// Shared trusted Pipeline operations. Loaded from the trusted implementation
// revision, never from a candidate checkout. Machine-specific resource
// budgets and stage orchestration stay in the entrypoint Jenkinsfiles.
import groovy.transform.Field
@Field def reporting

def configure(def reporter) { reporting = reporter }
def isCommitSha(String value) { value ==~ /[0-9a-fA-F]{40}/ }

def usePrebuiltConfig(String prefix) {
    def source = env.GKEYLL_CI_PREBUILT_CONFIG?.trim()
    if (!source) return false
    if (!source.startsWith('/') || !fileExists(source)) {
        error('GKEYLL_CI_PREBUILT_CONFIG must name an existing absolute config.mak path on the agent.')
    }
    if (!env.CI_RUN_DIR?.trim()) error('CI_RUN_DIR must be prepared before copying dependencies.')
    def launcher = env.PERSONAL_MPIEXEC?.trim()
        ?: (env.PERSONAL_MPI_HOME?.trim() ? "${env.PERSONAL_MPI_HOME.trim()}/bin/mpiexec" : '')
    // Keep the owner's launcher across candidate/baseline environment changes.
    if (!env.CI_PREBUILT_SOURCE_MPIEXEC) env.CI_PREBUILT_SOURCE_MPIEXEC = launcher ?: 'auto'
    launcher = env.CI_PREBUILT_SOURCE_MPIEXEC == 'auto' ? '' : env.CI_PREBUILT_SOURCE_MPIEXEC
    def dependencies = env.CI_BASELINE_DEPENDENCIES ?: "${env.CI_RUN_DIR}/dependencies"
    def runtime = [:]
    withEnv(["CI_PREBUILT_SOURCE=${source}", "CI_PREBUILT_PREFIX=${prefix}",
             "CI_PREBUILT_DEPENDENCIES=${dependencies}",
             "CI_PREBUILT_SCRIPT=${reporting.settings.prebuiltScript}",
             "CI_PREBUILT_MPIEXEC=${launcher}"]) {
        // Literal KEY=value lines avoid JSON APIs that need sandbox approval.
        def output = sh(returnStdout: true, script: '''#!/bin/bash
            set -euo pipefail
            if [[ -n "${CI_BASELINE_DEPENDENCIES:-}" ]]; then
                cmp "$CI_PREBUILT_SOURCE" "$CI_RUN_DIR/dependencies/source-config.mak" || {
                    echo 'Prebuilt config changed during this run; start a new CI run.' >&2
                    exit 1
                }
            fi
            python3 -I "$CI_PREBUILT_SCRIPT" --config "$CI_PREBUILT_SOURCE" \
                --prefix "$CI_PREBUILT_PREFIX" --dependencies "$CI_PREBUILT_DEPENDENCIES" \
                --mpiexec "${CI_PREBUILT_MPIEXEC:-}" --format env
        ''')
        output.readLines().each { line ->
            def fields = line.split('=', 2)
            runtime[fields[0]] = fields[1]
        }
    }
    env.LD_LIBRARY_PATH = runtime.LD_LIBRARY_PATH
    env.GKEYLL_CI_DEPENDENCY_ENV = runtime.GKEYLL_CI_DEPENDENCY_ENV
    if (runtime.MPI_HOME) {
        env.MPI_HOME = runtime.MPI_HOME
        env.OPAL_PREFIX = runtime.OPAL_PREFIX
        env.PERSONAL_MPI_HOME = runtime.MPI_HOME
    }
    if (runtime.MPIEXEC) env.PERSONAL_MPIEXEC = runtime.MPIEXEC
    writeFile file: 'alltargets.mak', text: '''.PHONY: all install
all: gkeyll
install: gkeyll
	$(MAKE) gkeyll-install
'''
    echo "Using dependency copy at ${dependencies}; installing Gkeyll in ${prefix}"
    return true
}

def checkoutPullRequest(String prNumber) {
    checkout([
        $class: 'GitSCM',
        branches: [[name: "refs/remotes/origin/candidate"]],
        doGenerateSubmoduleConfigurations: false,
        extensions: [
            [$class: 'CleanBeforeCheckout'],
            [$class: 'CloneOption', depth: 1, honorRefspec: true, noTags: true, shallow: true]
        ],
        userRemoteConfigs: [[
            refspec: "+refs/pull/${prNumber}/head:refs/remotes/origin/candidate",
            url: 'https://github.com/gkeyllorg/gkeyll.git'
        ]]
    ])
}

def checkoutRef(String selector, String localRef) {
    def refspec = isCommitSha(selector)
        ? "+${selector}:refs/remotes/origin/${localRef}"
        : "+refs/heads/${selector}:refs/remotes/origin/${localRef}"
    checkout([
        $class: 'GitSCM',
        branches: [[name: "refs/remotes/origin/${localRef}"]],
        doGenerateSubmoduleConfigurations: false,
        extensions: [
            [$class: 'CleanBeforeCheckout'],
            [$class: 'CloneOption', depth: 1, honorRefspec: true, noTags: true, shallow: true]
        ],
        userRemoteConfigs: [[
            refspec: refspec,
            url: 'https://github.com/gkeyllorg/gkeyll.git'
        ]]
    ])
}

def resolveRefCommit(String selector) {
    if (isCommitSha(selector)) return selector.toLowerCase()
    def commit = sh(returnStdout: true, script: "git ls-remote https://github.com/gkeyllorg/gkeyll.git refs/heads/${selector} | cut -f1 | head -n 1").trim()
    if (!(commit ==~ /[0-9a-f]{40}/)) error("Could not resolve baseline ref '${selector}' to a commit SHA.")
    return commit
}

def configureAndCompileCRegressions(String label, String driver, String sourceDir, String prefix, String moduleFile, String jobs, String runOptions = '', String configureOptions = '') {
    withEnv([
        "CI_REGRESSION_LABEL=${label}",
        "CI_REGRESSION_COMPILE_JOBS=${jobs}",
        "CI_REGRESSION_DRIVER=${driver}",
        "CI_REGRESSION_SOURCE_DIR=${sourceDir}",
        "CI_REGRESSION_PREFIX=${prefix}",
        "CI_REGRESSION_MODULE_FILE=${moduleFile}",
        "CI_REGRESSION_RUN_OPTIONS=${runOptions}",
        "CI_REGRESSION_CONFIGURE_OPTIONS=${configureOptions}",
        "CI_WORKSPACE=${env.WORKSPACE}"
    ]) {
        reporting.loggedSh("${label}-c-compile", '''#!/bin/bash -l
            set -euo pipefail
            . "$CI_REGRESSION_MODULE_FILE"
            if [ -n "${GKEYLL_CI_DEPENDENCY_ENV:-}" ]; then . "$GKEYLL_CI_DEPENDENCY_ENV"; fi
            "$CI_REGRESSION_DRIVER" runregression configure \
                --source-dir "$CI_REGRESSION_SOURCE_DIR" \
                --prefix "$CI_REGRESSION_PREFIX" ${CI_REGRESSION_CONFIGURE_OPTIONS:-}
            started="$(date +%s)"
            "$CI_REGRESSION_DRIVER" runregression run -c ${CI_REGRESSION_RUN_OPTIONS:-} --jobs "$CI_REGRESSION_COMPILE_JOBS" compile
            elapsed="$(( $(date +%s) - started ))"
            printf '%s\n' "$elapsed" > "$CI_WORKSPACE/$CI_REGRESSION_LABEL-c-compile-seconds.txt"
            echo "$CI_REGRESSION_LABEL C-regression compile time: $elapsed seconds"
        ''')
    }
}

def finalizeUnfinishedSlurmJob() {
    // Jenkins may terminate a durable-task process tree before its shell trap
    // completes. Run this from the Pipeline finally block as a second chance
    // to cancel a submitted job and preserve an abort-status artifact.
    sh '''#!/bin/bash -l
        set -euo pipefail

        for id_file in slurm-*-job-id.txt; do
            [[ -e "$id_file" ]] || continue
            status_file="${id_file%-id.txt}-status.txt"
            [[ -e "$status_file" ]] && continue

            job_id="$(<"$id_file")"
            if [[ ! "$job_id" =~ ^[0-9]+(_[0-9]+)?$ ]]; then
                printf 'job_id=%s\nstate=UNKNOWN\nexit_code=\n' "$job_id" > "$status_file"
                echo "WARNING: Cannot cancel malformed Slurm job ID '$job_id'" >&2
                continue
            fi

            echo "Pipeline finalizer: cancelling unfinished Slurm job $job_id"
            scancel "$job_id" || true
            while :; do
                state="$(squeue --noheader --jobs="$job_id" --format='%T' | head -n 1 || true)"
                if [[ -z "$state" ]]; then
                    break
                fi
                echo "Pipeline finalizer: waiting for job $job_id to stop (state: $state)"
                sleep 5
            done
            printf 'job_id=%s\nstate=CANCELLED_BY_JENKINS\nexit_code=\n' "$job_id" > "$status_file"
        done
    '''
}

return this
