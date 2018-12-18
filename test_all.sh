#!/bin/bash

# exit on first error, see: http://stackoverflow.com/a/185900/432509
error() {
    local parent_lineno="$1"
    local message="$2"
    local code="${3:-1}"
    if [[ -n "$message" ]] ; then
        echo "Error on or near line ${parent_lineno}: ${message}; exiting with status ${code}"
    else
        echo "Error on or near line ${parent_lineno}; exiting with status ${code}"
    fi
    exit "${code}"
}
trap 'error ${LINENO}' ERR
# done with trap

# support cargo command override
if [[ -z $CARGO_BIN ]]; then
    CARGO_BIN=cargo
fi

# toplevel git repo
ROOT=$(git rev-parse --show-toplevel)

for cargo_dir in $(find "$ROOT" -name Cargo.toml -printf '%h\n'); do
    echo "Running tests in: $cargo_dir"
    pushd "$cargo_dir"
    RUST_BACKTRACE=0 $CARGO_BIN test --no-run
    RUST_BACKTRACE=1 $CARGO_BIN test -- --nocapture
    popd
done
