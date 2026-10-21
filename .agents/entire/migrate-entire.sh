#!/usr/bin/env bash
#
# Migrate an existing Entire setup to a different checkpoint repository.
#
# Existing checkpoint data is COPIED and verified. Nothing is deleted
# from the source repository.
#
# Usage:
#   ./.agents/entire/migrate-entire.sh [options]
#
# Options:
#   --source REPO     Repository currently storing Entire data.
#                     Default: gkeyll
#
#   --dest REPO       Repository to migrate Entire data to.
#                     Default: gkeyll-agvcs
#
#   -h, --help        Show this help.
#
# Examples:
#   ./.agents/entire/migrate-entire.sh
#
#   ./.agents/entire/migrate-entire.sh \
#       --source gkeyll \
#       --dest gkeyll-agvcs
#

set -euo pipefail

ORG="gkeyllorg"
SOURCE="gkeyll"
DEST="gkeyll-agvcs"

usage()
{
  cat <<EOF
Usage:
  $0 [options]

Options:
  --source REPO     Repository currently storing Entire checkpoint data.
                    Default: gkeyll

  --dest REPO       Repository to migrate Entire checkpoint data to.
                    Default: gkeyll-agvcs

  -h, --help        Show this help.

Examples:
  $0
  $0 --source gkeyll --dest gkeyll-agvcs

The migration copies and verifies checkpoint data, then configures
future Entire checkpoints to use the destination repository.

It does NOT delete checkpoint data from the source repository.
EOF
}

fail()
{
  echo "ERROR: $*" >&2
  exit 1
}

# ----------------------------------------------------------------------
# Parse arguments
# ----------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source)
      [[ $# -ge 2 ]] || fail "--source requires an argument."
      SOURCE="$2"
      shift 2
      ;;
    --dest)
      [[ $# -ge 2 ]] || fail "--dest requires an argument."
      DEST="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      fail "Unknown option: $1. Use --help for usage."
      ;;
  esac
done

[[ "$SOURCE" == */* ]] || SOURCE="$ORG/$SOURCE"
[[ "$DEST"   == */* ]] || DEST="$ORG/$DEST"

[[ "$SOURCE" != "$DEST" ]] \
  || fail "Source and destination repositories are the same."

# ----------------------------------------------------------------------
# Basic checks
# ----------------------------------------------------------------------

command -v git >/dev/null 2>&1 \
  || fail "git is not installed."

command -v entire >/dev/null 2>&1 \
  || fail "Entire is not installed or is not in PATH."

ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" \
  || fail "Run this script from inside a Gkeyll checkout."

cd "$ROOT"

# ----------------------------------------------------------------------
# Construct repository URLs using the same authentication/transport
# style as origin.
#
# Examples supported:
#
#   https://github.com/user/gkeyll.git
#   https://janedoe@github.com/user/gkeyll.git
#   git@github.com:user/gkeyll.git
#   git@github-work:user/gkeyll.git
#   ssh://git@github.com/user/gkeyll.git
# ----------------------------------------------------------------------

ORIGIN="$(git remote get-url origin 2>/dev/null || true)"
[[ -n "$ORIGIN" ]] || fail "No Git remote named 'origin' was found."

case "$ORIGIN" in
  git@*:*)
    SSH_HOST="${ORIGIN%%:*}"

    repo_url()
    {
      printf '%s:%s.git\n' "$SSH_HOST" "$1"
    }
    ;;

  ssh://*)
    REST="${ORIGIN#ssh://}"
    SSH_AUTHORITY="${REST%%/*}"

    repo_url()
    {
      printf 'ssh://%s/%s.git\n' "$SSH_AUTHORITY" "$1"
    }
    ;;

  http://*github.com/*|https://*github.com/*)
    if [[ "$ORIGIN" =~ ^(https?://)([^/]+@)?github\.com/ ]]; then
      HTTPS_BASE="${BASH_REMATCH[1]}${BASH_REMATCH[2]:-}github.com"

      repo_url()
      {
        printf '%s/%s.git\n' "$HTTPS_BASE" "$1"
      }
    else
      fail "Could not parse GitHub HTTPS origin: $ORIGIN"
    fi
    ;;

  *)
    fail "Could not determine how to access GitHub from origin: $ORIGIN"
    ;;
esac

SOURCE_URL="$(repo_url "$SOURCE")"
DEST_URL="$(repo_url "$DEST")"

# ----------------------------------------------------------------------
# Check access.
# ----------------------------------------------------------------------

git ls-remote "$SOURCE_URL" >/dev/null 2>&1 \
  || fail "Cannot read $SOURCE."

git ls-remote "$DEST_URL" >/dev/null 2>&1 \
  || fail "Cannot access $DEST. Request access from its administrators."

git push --dry-run --no-verify \
  "$DEST_URL" HEAD:refs/entire/access-check/migrate \
  >/dev/null 2>&1 \
  || fail "Write access to $DEST is required."

# Make migration into public Gkeyll an explicit choice.
if [[ "$DEST" == "$ORG/gkeyll" ]]; then
  echo
  echo "WARNING: Entire data migrated to $DEST will be PUBLIC."
  echo

  read -r -p "Continue? [y/N] " reply
  [[ "$reply" =~ ^[Yy]$ ]] || exit 1
fi

echo
echo "Migrating Entire data:"
echo
echo "  $SOURCE"
echo "      ->"
echo "  $DEST"
echo

TMP="$(mktemp -d "${TMPDIR:-/tmp}/entire-migrate.XXXXXX")"

cleanup()
{
  git for-each-ref --format='delete %(refname)' refs/migrate-entire/ \
    | git update-ref --stdin 2>/dev/null || true
  rm -rf "$TMP"
}
trap cleanup EXIT

# ----------------------------------------------------------------------
# Copy modern refs from SOURCE using a temporary local namespace.
# ----------------------------------------------------------------------

git ls-remote "$SOURCE_URL" 'refs/entire/checkpoints/*' \
  | sort -u > "$TMP/source"

if [[ -s "$TMP/source" ]]; then
  echo "Copying checkpoint refs..."

  git fetch "$SOURCE_URL" \
    '+refs/entire/checkpoints/*:refs/migrate-entire/checkpoints/*'

  git push --no-verify "$DEST_URL" \
    'refs/migrate-entire/checkpoints/*:refs/entire/checkpoints/*'
fi

# ----------------------------------------------------------------------
# Handle legacy Entire storage.
# ----------------------------------------------------------------------

LEGACY_SHA="$(
  git ls-remote "$SOURCE_URL" refs/heads/entire/checkpoints/v1 \
    | awk 'NR == 1 {print $1}'
)"

if [[ -n "$LEGACY_SHA" ]]; then
  echo "Legacy entire/checkpoints/v1 storage detected."

  LOCAL_LEGACY="$(
    git rev-parse --verify refs/heads/entire/checkpoints/v1 \
      2>/dev/null || true
  )"

  if [[ -n "$LOCAL_LEGACY" &&
        "$LOCAL_LEGACY" != "$LEGACY_SHA" ]]; then
    fail "Local entire/checkpoints/v1 differs from $SOURCE."
  fi

  if [[ -z "$LOCAL_LEGACY" ]]; then
    git fetch "$SOURCE_URL" \
      refs/heads/entire/checkpoints/v1:refs/heads/entire/checkpoints/v1
  fi

  echo "Converting legacy checkpoints to git-refs..."
  entire doctor migrate-checkpoints
fi

# ----------------------------------------------------------------------
# Include any checkpoints currently present only in this checkout.
# ----------------------------------------------------------------------

if git for-each-ref --format='%(refname)' \
     refs/entire/checkpoints/ | grep -q .; then

  git push --no-verify "$DEST_URL" \
    'refs/entire/checkpoints/*:refs/entire/checkpoints/*'
fi

# ----------------------------------------------------------------------
# Verify all expected checkpoint refs and object IDs.
# ----------------------------------------------------------------------

cp "$TMP/source" "$TMP/expected"

git for-each-ref \
  --format='%(objectname) %(refname)' \
  refs/entire/checkpoints/ |
while read -r oid ref; do
  printf '%s\t%s\n' "$oid" "$ref"
done >> "$TMP/expected"

sort -u -o "$TMP/expected" "$TMP/expected"

git ls-remote "$DEST_URL" 'refs/entire/checkpoints/*' \
  | sort -u > "$TMP/destination"

comm -23 "$TMP/expected" "$TMP/destination" > "$TMP/missing"

if [[ -s "$TMP/missing" ]]; then
  echo
  echo "The following checkpoints were not verified at $DEST:" >&2
  cat "$TMP/missing" >&2
  fail "Migration verification failed."
fi

# ----------------------------------------------------------------------
# Point future Entire checkpoints at DEST.
# ----------------------------------------------------------------------

entire configure \
  --checkpoint-remote "github:$DEST"

echo
echo "Migration complete."
echo
echo "Future Entire checkpoint data will use:"
echo
echo "  $DEST"
echo
echo "No checkpoint data was deleted from:"
echo
echo "  $SOURCE"
echo
echo "Verify the migrated history before removing old source refs."
echo
entire status
