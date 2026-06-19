#!/usr/bin/env bash
# Copies the bundled sample data out of the scNexus image onto the host, so it
# can be uploaded through the browser-based file pickers in Modules II–IV.
#
# Usage: ./extract-test-data.sh [destination-dir]
#   destination-dir  where to copy the data (default: ~/Downloads)
#
# Override the image name with IMAGE=... if you tagged the build differently.
set -euo pipefail

IMAGE="${IMAGE:-scnexus-demo:latest}"
DEST="${1:-$HOME/Downloads}"
SRC_PATH="/srv/app/test-data"

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "Image '$IMAGE' not found. Build it first (docker compose build) or set IMAGE=..." >&2
  exit 1
fi

mkdir -p "$DEST"

# docker cp needs a container (not an image), so create a stopped one, copy, remove.
container_id="$(docker create "$IMAGE")"
trap 'docker rm -f "$container_id" >/dev/null 2>&1 || true' EXIT

docker cp "$container_id:$SRC_PATH" "$DEST/scNexus-test-data"

echo "Test data copied to: $DEST/scNexus-test-data"
