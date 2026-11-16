#!/usr/bin/env bash
# Install the PRESTUS Nexus (macOS, Linux) in one step:
#   bash nexus/install.sh            # from the PRESTUS folder, or from anywhere
#   bash nexus/install.sh --yes      # install Miniforge without asking, if conda is missing
#
# Uses conda/mamba if it is installed; otherwise offers to install Miniforge
# (conda-forge's conda) into ~/miniforge3. Then creates the environment
# 'prestus-nexus' from environment.yml, or updates it if it exists. Run it
# again after a PRESTUS update to update the environment.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env_name="prestus-nexus"
assume_yes=false
for arg in "$@"; do
    case "$arg" in
        -y|--yes) assume_yes=true ;;
        -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done

find_conda() {
    for c in mamba conda; do
        if command -v "$c" >/dev/null 2>&1; then command -v "$c"; return; fi
    done
    for d in miniforge3 mambaforge miniconda3 anaconda3; do
        if [ -x "$HOME/$d/bin/conda" ]; then echo "$HOME/$d/bin/conda"; return; fi
    done
}

conda_bin="$(find_conda || true)"
if [ -z "$conda_bin" ]; then
    echo "conda was not found. The Nexus needs it once, to install its Python environment."
    if [ "$assume_yes" != true ]; then
        read -r -p "Install Miniforge into $HOME/miniforge3 now? [y/N] " answer
        case "$answer" in [yY]*) ;; *) echo "Stopped. Install conda yourself, then run this again."; exit 1 ;; esac
    fi
    installer="$(mktemp -t miniforge.XXXXXX).sh"
    url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    echo "Downloading $url"
    curl -fsSL "$url" -o "$installer"
    bash "$installer" -b -p "$HOME/miniforge3"
    rm -f "$installer"
    conda_bin="$HOME/miniforge3/bin/conda"
fi
echo "Using $conda_bin"

cd "$here"
if "$conda_bin" env list | awk '{print $1}' | grep -qx "$env_name"; then
    echo "Updating the environment '$env_name'..."
    "$conda_bin" env update -n "$env_name" -f environment.yml
else
    echo "Creating the environment '$env_name'..."
    "$conda_bin" env create -n "$env_name" -f environment.yml
fi

env_dir="$("$conda_bin" env list | awk -v n="$env_name" '$1 == n {print $NF}')"
python="$env_dir/bin/python"
if "$python" -c "import yaml, PyQt6.QtWidgets, prestus_nexus" 2>/dev/null; then
    echo
    echo "PRESTUS Nexus installed."
    echo "  Start it in MATLAB:   prestus_nexus"
    echo "  or in a terminal:     $env_dir/bin/prestus-nexus"
else
    echo "The environment was created, but $python cannot import the Nexus." >&2
    echo "See 'Installation' in doc/doc_nexus.md." >&2
    exit 1
fi
