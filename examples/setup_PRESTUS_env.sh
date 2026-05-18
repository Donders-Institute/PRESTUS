#!/usr/bin/env bash
# setup_PRESTUS_env.sh
#
# Creates a conda environment (PRESTUS_env_4.6.0) containing SimNIBS 4.6.0
# and all optional dependencies required for PlanTUS transducer placement.
#
# Usage
# -----
#   bash setup_PRESTUS_env.sh [--no-plantus] [--no-workbench]
#
#   --no-plantus    Skip nilearn / vtk / h5py (PlanTUS python deps)
#   --no-workbench  Skip Connectome Workbench (wb_command)
#                   Use this if wb_command is already available as an HPC module
#
# After installation, point PRESTUS at the new environment:
#   parameters.simnibs_bin_path      = '<HOME>/.conda/envs/PRESTUS_env_4.6.0/bin'
#   parameters.placement.plantus.env_path = '<HOME>/.conda/envs/PRESTUS_env_4.6.0/bin'
#
# Requirements: conda (Anaconda / Miniconda), internet access

set -euo pipefail

ENV_NAME="PRESTUS_env_4.6.0"
SIMNIBS_VERSION="4.6.0"
SIMNIBS_ENV_URL="https://github.com/simnibs/simnibs/releases/download/v${SIMNIBS_VERSION}/environment_linux.yml"
SIMNIBS_WHL_URL="https://github.com/simnibs/simnibs/releases/download/v${SIMNIBS_VERSION}/simnibs-${SIMNIBS_VERSION}-cp311-cp311-linux_x86_64.whl"

INSTALL_PLANTUS=true
INSTALL_WORKBENCH=true

for arg in "$@"; do
    case $arg in
        --no-plantus)   INSTALL_PLANTUS=false ;;
        --no-workbench) INSTALL_WORKBENCH=false ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

echo "============================================================"
echo " PRESTUS environment setup: ${ENV_NAME}"
echo "============================================================"

# ── 1. Download SimNIBS conda environment spec ────────────────────────────────
TMP_YML=$(mktemp /tmp/environment_simnibs_XXXXXX.yml)
echo "[1/5] Downloading SimNIBS ${SIMNIBS_VERSION} environment spec..."
wget -q "${SIMNIBS_ENV_URL}" -O "${TMP_YML}"

# ── 2. Create conda environment ───────────────────────────────────────────────
echo "[2/5] Creating conda environment '${ENV_NAME}'..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "  Environment '${ENV_NAME}' already exists — skipping creation."
    echo "  To recreate it, first run: conda env remove -n ${ENV_NAME}"
else
    conda env create -f "${TMP_YML}" -n "${ENV_NAME}"
fi
rm -f "${TMP_YML}"

# ── 3. Install SimNIBS ────────────────────────────────────────────────────────
echo "[3/5] Installing SimNIBS ${SIMNIBS_VERSION}..."
conda run -n "${ENV_NAME}" pip install "${SIMNIBS_WHL_URL}"

# ── 4. Install PlanTUS Python dependencies (optional) ────────────────────────
if [ "${INSTALL_PLANTUS}" = true ]; then
    echo "[4/5] Installing PlanTUS dependencies (nilearn, vtk, h5py)..."
    conda run -n "${ENV_NAME}" pip install nilearn vtk h5py
else
    echo "[4/5] Skipping PlanTUS dependencies (--no-plantus)."
fi

# ── 5. Install Connectome Workbench (optional) ────────────────────────────────
if [ "${INSTALL_WORKBENCH}" = true ]; then
    echo "[5/5] Installing Connectome Workbench (wb_command)..."
    conda run -n "${ENV_NAME}" conda install -y -c conda-forge connectome-workbench
else
    echo "[5/5] Skipping Connectome Workbench (--no-workbench)."
    echo "  If wb_command is available as an HPC module, set connectome_wb_path in your config:"
    echo "    parameters.placement.plantus.connectome_wb_path = '/path/to/workbench/bin_linux64';"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
PYTHON_BIN=$(conda run -n "${ENV_NAME}" which python)
ENV_BIN=$(dirname "${PYTHON_BIN}")

echo ""
echo "============================================================"
echo " Setup complete: ${ENV_NAME}"
echo "============================================================"
echo ""
echo "Add the following to your PRESTUS config or MATLAB script:"
echo ""
echo "  parameters.simnibs_bin_path = '${ENV_BIN}';"
echo "  parameters.placement.plantus.env_path = '${ENV_BIN}';"
echo ""
echo "To verify the installation:"
echo "  conda activate ${ENV_NAME}"
echo "  python -c \"import simnibs; import nilearn; print('OK')\""
if [ "${INSTALL_WORKBENCH}" = true ]; then
    echo "  which wb_command"
fi
