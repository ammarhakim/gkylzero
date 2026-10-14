#!/usr/bin/env bash
# Machine defaults for the shared cluster controller and Jenkins client.
set -euo pipefail
readonly SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly SCRIPT_DIR="${SCRIPT_PATH%/*}"
CI_ROOT="${GKEYLL_CI_ROOT:-/scratch/gpfs/${USER:?USER must be set}/gkeyll_ci}"
JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64}"
JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-stellar_cpu}"
source "$SCRIPT_DIR/jenkins-cluster.sh"
main "$@"
