#!/usr/bin/env bash
# 无泵站算例 — 与上级目录共享 dem/parameters/rainfall 等, 仅 Drainage.inp 不同
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(cd "$CASE_DIR/.." && pwd)"
SERGHEI_BIN="${SERGHEI_BIN:-$(cd "$PARENT_DIR/../.." && pwd)/bin/serghei}"
OUT_DIR="${1:-$CASE_DIR/output}"

if [[ -d "$OUT_DIR" ]]; then
  echo "[run] removing existing output: $OUT_DIR"
  rm -rf "$OUT_DIR"
fi

# serghei 从单一输入目录读文件: 用上级目录的地表输入 + 本目录排水拓扑
LINK_DIR="$CASE_DIR/_run_input"
rm -rf "$LINK_DIR"
mkdir -p "$LINK_DIR"
for f in dem.input parameters.input sw.input rainfall.input Drainage-parameter.input; do
  ln -sf "$PARENT_DIR/$f" "$LINK_DIR/$f"
done
cp "$CASE_DIR/Drainage.inp" "$LINK_DIR/Drainage.inp"

echo "[run] input : $LINK_DIR/ (Drainage from no_pump/, surface from parent)"
echo "[run] output: $OUT_DIR/"
echo "[run] binary: $SERGHEI_BIN"

exec "$SERGHEI_BIN" "$LINK_DIR/" "$OUT_DIR/" 1
