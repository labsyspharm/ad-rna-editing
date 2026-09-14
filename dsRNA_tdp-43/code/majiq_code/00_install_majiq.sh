#!/usr/bin/env bash
# ============================================================================
# 00_install_majiq.sh -- install MAJIQ + VOILA (academic) into a conda env.
#
# MAJIQ is a compiled C++/Cython tool linked against HTSlib. It is NOT on
# PyPI/conda directly; it is pip-installed from source with HTSLIB_* pointing
# at an htslib you provide. We get htslib from bioconda so nothing is compiled
# by hand.
#
# Needs a C/C++ toolchain (provided here via conda compilers) and Python 3.8-3.11.
# Run once. Requires conda or mamba on PATH.
# ============================================================================
set -euo pipefail
ENV_NAME="${1:-majiq}"
CONDA="$(command -v mamba || command -v conda)"; [ -n "$CONDA" ] || { echo "need conda/mamba"; exit 1; }

# 1) env with python 3.11, htslib, and a compiler toolchain
"$CONDA" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    "python=3.11" "htslib>=1.15" pysam numpy cython pip \
    c-compiler cxx-compiler zlib

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# 2) point the MAJIQ build at the conda htslib
export HTSLIB_LIBRARY_DIR="${CONDA_PREFIX}/lib"
export HTSLIB_INCLUDE_DIR="${CONDA_PREFIX}/include"
python -m pip install --upgrade pip setuptools wheel

# 3) install MAJIQ academic. Primary: maintained GitHub fork; fallback: biociphers.
pip install "git+https://github.com/OncoHarmony-Network/majiq_academic.git@v2.5.7" \
  || pip install "git+https://bitbucket.org/biociphers/majiq_academic.git@v2.5.7"

# 4) sanity: these should print help (build/deltapsi need the LICENSE at runtime)
majiq --version || true
voila --help >/dev/null 2>&1 && echo "voila OK" || echo "voila not on PATH?"

cat <<EOF

============================================================================
Installed MAJIQ into conda env: ${ENV_NAME}
Before running the pipeline, in EVERY shell:

    conda activate ${ENV_NAME}
    export HTSLIB_LIBRARY_DIR="${CONDA_PREFIX}/lib"     # runtime linking
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"
    export MAJIQ_LICENSE_FILE=/path/to/your_academic.lic   # REQUIRED

Get the academic license file from https://majiq.biociphers.org/ (I cannot).
============================================================================
EOF
