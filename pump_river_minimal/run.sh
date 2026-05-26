#!/usr/bin/env bash
# 最小泵+河位算例 — 地表(SWE)+排水耦合, 无地下水
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SERGHEI_BIN="${SERGHEI_BIN:-$(cd "$CASE_DIR/../.." && pwd)/bin/serghei}"
OUT_DIR="${1:-$CASE_DIR/output_with_pump}"

if [[ -d "$OUT_DIR" ]]; then
  echo "[run] removing existing output: $OUT_DIR"
  rm -rf "$OUT_DIR"
fi

echo "[run] input : $CASE_DIR/"
echo "[run] output: $OUT_DIR/"
echo "[run] binary: $SERGHEI_BIN"

exec "$SERGHEI_BIN" "$CASE_DIR/" "$OUT_DIR/" 1
