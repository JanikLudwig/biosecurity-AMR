#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${project_dir}/.tools/amrfinder"
environment_file="${project_dir}/environments/amrfinder.yml"

if [[ ! -x "${prefix}/bin/amrfinder" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    micromamba create --yes --prefix "${prefix}" --file "${environment_file}"
  elif command -v mamba >/dev/null 2>&1; then
    mamba env create --yes --prefix "${prefix}" --file "${environment_file}"
  elif command -v conda >/dev/null 2>&1; then
    conda env create --yes --prefix "${prefix}" --file "${environment_file}"
  else
    echo "No micromamba, mamba, or conda executable was found." >&2
    echo "Install micromamba (recommended), then rerun this script." >&2
    exit 1
  fi
fi

CONDA_PREFIX="${prefix}" "${prefix}/bin/amrfinder" --version

database_dir="${prefix}/share/amrfinderplus/data/latest"
if [[ ! -e "${database_dir}/version.txt" ]]; then
  echo "Downloading the AMRFinderPlus reference database..."
  CONDA_PREFIX="${prefix}" "${prefix}/bin/amrfinder" -u
fi

echo "Database version: $(sed -n '1p' "${database_dir}/version.txt")"
