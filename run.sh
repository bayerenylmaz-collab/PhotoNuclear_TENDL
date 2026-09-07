#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

echo
echo "============================================"
echo " Photonuclear Katalog - kolay çalıştırıcı"
echo "============================================"
echo

if ! command -v python3 >/dev/null 2>&1; then
  echo "[HATA] python3 bulunamadı. Python 3.10+ kurun."
  exit 1
fi

if [[ ! -x .venv/bin/python ]]; then
  echo "Sanal ortam oluşturuluyor..."
  python3 -m venv .venv
fi

# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install -q --upgrade pip
python -m pip install -q -e .
echo
python -m photonuclear_catalog --interactive
