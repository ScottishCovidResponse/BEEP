#!/bin/bash

set -e
set -u

outfile=$1
tmpfile=${outfile}.tmp

version=$(git describe --dirty --always)
cat >"$tmpfile" <<EOF
#pragma once
#define GIT_VERSION "$version"
EOF

if [ -r $outfile ]; then
    if diff -q "$outfile" "$tmpfile" >/dev/null; then
        rm -f "$tmpfile"
        # No change to version
        exit 0
    fi
fi

mv "$tmpfile" "$outfile"
