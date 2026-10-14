#!/usr/bin/env bash
# Machine defaults for the shared cluster controller and Jenkins client.
set -euo pipefail
readonly SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly SCRIPT_DIR="${SCRIPT_PATH%/*}"
CI_ROOT="${GKEYLL_CI_ROOT:-}"
JAVA_HOME="${JAVA_HOME:-}"
JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-perlmutter_gpu}"
source "$SCRIPT_DIR/jenkins-cluster.sh"
main "$@"
