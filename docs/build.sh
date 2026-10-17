#!/usr/bin/env bash
# 把 docs/ 下所有 .typ 编译成 PDF，输出到 docs/build/。
#
#   ./docs/build.sh              # 编译全部
#   ./docs/build.sh theory       # 只编译 docs/theory/ 下的
#
# 注意 --root：Typst 默认把“输入文件所在目录”当作项目根，而各篇文档都用
# `../template.typ` 引用上级目录，所以必须把根显式设成仓库根目录。
# 图片路径统一写成相对 docs/ 的形式（`assets/...`），因为 fig()/fig-gallery() 定义在
# docs/template.typ 里，Typst 按“路径字符串所在文件”解析相对路径。

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

TYPST="${TYPST:-typst}"
OUT="$REPO_ROOT/docs/build"
mkdir -p "$OUT"

shopt -s nullglob
targets=()
if [ $# -gt 0 ]; then
  for sub in "$@"; do
    targets+=("docs/$sub"/*.typ)
  done
else
  targets+=(docs/*.typ docs/*/*.typ)
fi

fail=0
for f in "${targets[@]}"; do
  # template.typ 只是排版模板，不是文档，跳过
  [ "$(basename "$f")" = "template.typ" ] && continue
  rel="${f#docs/}"
  out="$OUT/${rel%.typ}.pdf"
  mkdir -p "$(dirname "$out")"
  if "$TYPST" compile --root . "$f" "$out"; then
    echo "  ok  $f"
  else
    echo "  FAIL $f"
    fail=1
  fi
done

exit "$fail"
