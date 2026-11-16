#!/usr/bin/env bash
# Build ymf's sdist and wheel, check them, install the wheel into a fresh venv
# and run it. The same script runs in release.yml and by hand:
#
#   tools/check_dist.sh [OUTDIR]        # default: a new temporary directory
#
# Needs `python -m build` and `twine` (pip install build twine). Nothing is
# uploaded. The sdist must stay small (limit below) and carry no logs, HDF5
# data, compiled files or CI files; a file that matters to the tests is
# tracked in git, so a failure here is a stray file or a missing prune in
# MANIFEST.in.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$(mktemp -d)}"
LIMIT_KB="${YMF_SDIST_LIMIT_KB:-1024}"
mkdir -p "$OUT"
PY="${PYTHON:-python}"

"$PY" -m build --outdir "$OUT/dist" "$ROOT"
"$PY" -m twine check --strict "$OUT"/dist/*

SDIST="$(ls "$OUT"/dist/ymf-*.tar.gz)"
WHEEL="$(ls "$OUT"/dist/ymf-*.whl)"
KB=$(( $(wc -c < "$SDIST") / 1024 ))
echo "sdist: $SDIST (${KB} KB, limit ${LIMIT_KB} KB)"
echo "wheel: $WHEEL ($(( $(wc -c < "$WHEEL") / 1024 )) KB)"
[ "$KB" -le "$LIMIT_KB" ] || { echo "sdist too large"; exit 1; }

# No accidental payload.
if tar tzf "$SDIST" | grep -E '\.(log|h5|so|pyc|o)$|\.coverage|/\.github/|\.DS_Store'; then
  echo "unexpected files in the sdist (listed above)"; exit 1
fi

# The wheel carries the package and nothing else but its metadata.
if unzip -Z1 "$WHEEL" | grep -vE '^ymf/|^ymf-[^/]*\.dist-info/'; then
  echo "unexpected files in the wheel (listed above)"; exit 1
fi

# Fresh venv, wheel only, run from a directory without the checkout on the path.
"$PY" -m venv "$OUT/venv"
"$OUT/venv/bin/python" -m pip install --quiet "${WHEEL}[frontend,examples]"
mkdir -p "$OUT/run"
cd "$OUT/run"
"$OUT/venv/bin/python" -I - "$ROOT" <<'PYEOF'
import sys, pathlib, ymf
root = pathlib.Path(sys.argv[1])
print("import ymf", ymf.__version__, "from", ymf.__file__)
assert "site-packages" in ymf.__file__, "not the installed wheel"
assert ymf.__version__ not in ("0.0.0", "0.0.0+unknown"), "no version"
doc = ymf.load_ymf(root / "examples" / "poisson.ymf")
print("loaded examples/poisson.ymf:", type(doc).__name__)
PYEOF
"$OUT/venv/bin/ymf2xmf" --help > /dev/null
"$OUT/venv/bin/python" "$ROOT/examples/write_archive.py" > /dev/null
"$OUT/venv/bin/ymf2xmf" ymf-example-output/heat.ymf -o heat-check.xmf
echo "ok: build, twine check, wheel install, import, spec load, ymf2xmf"
echo "$OUT"
