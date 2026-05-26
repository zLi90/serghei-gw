#!/usr/bin/env bash
# 依次运行: 有泵(时变河位) / 无泵(时变河位), 输出到 output_with_pump 与 output_no_pump
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SERGHEI_BIN="${SERGHEI_BIN:-$(cd "$CASE_DIR/../.." && pwd)/bin/serghei}"

if [[ ! -x "$SERGHEI_BIN" ]]; then
  echo "[run_compare] 请先编译: SERGHEI_BIN=$SERGHEI_BIN 不存在或不可执行"
  exit 1
fi

export SERGHEI_BIN

echo "========== 有泵 + 时变河位 =========="
"$CASE_DIR/run.sh" "$CASE_DIR/output_with_pump"

echo ""
echo "========== 无泵 + 时变河位 =========="
"$CASE_DIR/no_pump/run.sh" "$CASE_DIR/output_no_pump"

echo ""
echo "[run_compare] 完成."
echo "  有泵: $CASE_DIR/output_with_pump/DrainageTimeSeries.out"
echo "  无泵: $CASE_DIR/output_no_pump/DrainageTimeSeries.out"
