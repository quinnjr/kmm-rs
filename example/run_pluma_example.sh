#!/usr/bin/env bash
# Run the KMM example through a real pluma binary.
#
# Steps:
#   1. Build libkmm_rs.so in release mode.
#   2. Stage it under a temporary PLUMA_PLUGIN_PATH.
#   3. Invoke pluma against example/config.txt from the repo root.
#   4. Diff the produced output against example/expected_output.tsv.
#
# Usage:
#     ./run_pluma_example.sh                # uses `pluma` from $PATH
#     ./run_pluma_example.sh /path/to/pluma # explicit binary
#     PLUMA=/path/to/pluma ./run_pluma_example.sh
#
# Requirements:
#   * `pluma` binary on $PATH or via $1/$PLUMA. Must be built with
#     `--with-rust` and the HAVE_RUST fix from FIUBioRG/PluMA#13.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PLUMA_BIN="${1:-${PLUMA:-pluma}}"

if ! command -v "$PLUMA_BIN" >/dev/null 2>&1 && ! [ -x "$PLUMA_BIN" ]; then
    echo "error: pluma binary not found ('$PLUMA_BIN'). Pass a path or set \$PLUMA." >&2
    exit 1
fi

cd "$ROOT"
echo "→ cargo build --release"
cargo build --release --quiet

PLUGIN_DIR="$(mktemp -d)"
trap 'rm -rf "$PLUGIN_DIR"' EXIT
mkdir -p "$PLUGIN_DIR/KMM"
ln -sfn "$ROOT/target/release/libkmm_rs.so" "$PLUGIN_DIR/KMM/libKMMPlugin.so"
ln -sfn "$ROOT/Cargo.toml"                  "$PLUGIN_DIR/KMM/Cargo.toml"

mkdir -p example/out
rm -f example/out/results.tsv

echo "→ pluma example/config.txt"
PLUMA_PLUGIN_PATH="$PLUGIN_DIR" "$PLUMA_BIN" example/config.txt

echo
echo "→ diff example/out/results.tsv example/expected_output.tsv"
if diff -q example/out/results.tsv example/expected_output.tsv >/dev/null; then
    echo "✓ output matches expected_output.tsv"
else
    echo "✗ output differs from expected_output.tsv"
    diff example/out/results.tsv example/expected_output.tsv | head -20
    exit 1
fi
