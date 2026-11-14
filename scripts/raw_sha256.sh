#!/bin/sh
# Reproducible source-identity audit. Runs after Stata source reconstruction.
# macOS /usr/bin/shasum is bundled with the system.
set -eu
raw="$1"
out="$2"
if [ ! -d "$raw" ]; then
    printf 'RAW SOURCE DIRECTORY MISSING: %s\n' "$raw" > "$out"
    exit 1
fi
if [ -x /usr/bin/shasum ]; then
    hash=/usr/bin/shasum
elif command -v shasum >/dev/null 2>&1; then
    hash=shasum
else
    printf 'SHA256 NOT PRODUCED: shasum unavailable\n' > "$out"
    exit 1
fi
printf 'SOURCE FILE SHA-256 / source date and collection can change across vintages\n' > "$out"
find "$raw" -type f -exec "$hash" -a 256 '{}' ';' >> "$out"
