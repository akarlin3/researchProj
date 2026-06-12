#!/usr/bin/env bash
# Install the Proteus LOCAL toolchain so the test suite runs instead of skipping.
# Shared by the Claude-Code-on-the-web SessionStart hook (.claude/hooks/session-start.sh)
# and CI (.github/workflows/ci.yml). Idempotent + non-interactive.
#
# Installs: Python deps (numpy/pyyaml/biotite/biopython/pytest/ruff/vina), MMseqs2 +
# Foldseek static binaries, fpocket (built from source), Open Babel, and the control
# structures. Set PROTEUS_FETCH_WEIGHTS=1 to also fetch the ~2.4 GB ProstT5 weights
# (needed by the S1/S2/pipeline tests) — skipped by default so CI stays fast.
#
# Tools land in $PROTEUS_TOOLS_DIR (default ~/.proteus-tools); add $PROTEUS_TOOLS_DIR/bin
# to PATH and export PROTEUS_PROSTT5_MODEL=$PROTEUS_TOOLS_DIR/prostt5 to use them.
set -uo pipefail

TOOLS_DIR="${PROTEUS_TOOLS_DIR:-$HOME/.proteus-tools}"
BIN="$TOOLS_DIR/bin"
REPO_DIR="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
mkdir -p "$BIN"

log() { echo "[setup_tools] $*"; }
have() { command -v "$1" >/dev/null 2>&1; }

# --- Python dependencies (pure-Python pipeline + test deps) ------------------ #
log "installing Python dependencies …"
python3 -m pip install --quiet \
    numpy pyyaml biotite biopython pytest ruff vina \
    || log "WARN: some Python deps failed to install"

# --- MMseqs2 + Foldseek (static linux-avx2 binaries) ------------------------- #
fetch_tar() {  # url dest_subdir bin_relpath linkname
    local url="$1" sub="$2" rel="$3" name="$4"
    if have "$name" || [ -x "$BIN/$name" ]; then return 0; fi
    log "installing $name …"
    if wget -qO "$TOOLS_DIR/$name.tar.gz" "$url"; then
        tar xzf "$TOOLS_DIR/$name.tar.gz" -C "$TOOLS_DIR" && ln -sf "$TOOLS_DIR/$rel" "$BIN/$name"
        rm -f "$TOOLS_DIR/$name.tar.gz"
    else
        log "WARN: could not download $name"
    fi
}
fetch_tar "https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz" mmseqs "mmseqs/bin/mmseqs" mmseqs
fetch_tar "https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz" foldseek "foldseek/bin/foldseek" foldseek

# --- fpocket (build from source) -------------------------------------------- #
if ! have fpocket && [ ! -x "$BIN/fpocket" ]; then
    log "building fpocket from source …"
    rm -rf "$TOOLS_DIR/fpocket"
    # NB: fpocket's Makefile is not parallel-safe — build serially (no -j).
    if git clone --depth 1 https://github.com/Discngine/fpocket.git "$TOOLS_DIR/fpocket" 2>/dev/null \
        && make -C "$TOOLS_DIR/fpocket" >/dev/null 2>&1; then
        ln -sf "$TOOLS_DIR/fpocket/bin/fpocket" "$BIN/fpocket"
    else
        log "WARN: fpocket build failed (need gcc/make) — S4/S5 tests will skip"
    fi
fi

# --- Open Babel (receptor/ligand prep for docking) -------------------------- #
if ! have obabel; then
    log "installing Open Babel …"
    SUDO=""; [ "$(id -u)" -ne 0 ] && have sudo && SUDO="sudo"
    ($SUDO apt-get update -qq && $SUDO apt-get install -y -qq openbabel) >/dev/null 2>&1 \
        || log "WARN: could not apt-install openbabel — live docking test will skip"
fi

export PATH="$BIN:$PATH"

# --- Control structures (S4/S5/calibration/recovery/docking-live) ----------- #
if [ -f "$REPO_DIR/controls/fetch_controls.py" ]; then
    log "fetching control structures …"
    python3 "$REPO_DIR/controls/fetch_controls.py" --out "$REPO_DIR/structures" \
        >/dev/null 2>&1 || log "WARN: control fetch failed (network?) — structure tests will skip"
fi

# --- ProstT5 weights (optional; large) -------------------------------------- #
WEIGHTS_DIR="$TOOLS_DIR/prostt5"
if [ "${PROTEUS_FETCH_WEIGHTS:-0}" = "1" ]; then
    if [ -e "$WEIGHTS_DIR/prostt5-f16.gguf" ]; then
        log "ProstT5 weights already present"
    elif have foldseek || [ -x "$BIN/foldseek" ]; then
        log "downloading ProstT5 weights (~2.4 GB, first run only) …"
        "$BIN/foldseek" databases ProstT5 "$WEIGHTS_DIR" "$TOOLS_DIR/pt5_tmp" >/dev/null 2>&1 \
            && rm -rf "$TOOLS_DIR/pt5_tmp" \
            || log "WARN: ProstT5 weights download failed — S1/S2 tests will skip"
    fi
fi

log "done. tools in $BIN"
have mmseqs   && log "  mmseqs:   $(command -v mmseqs)"   || log "  mmseqs:   (missing)"
have foldseek && log "  foldseek: $(command -v foldseek)" || log "  foldseek: (missing)"
have fpocket  && log "  fpocket:  $(command -v fpocket)"  || log "  fpocket:  (missing)"
have obabel   && log "  obabel:   $(command -v obabel)"   || log "  obabel:   (missing)"
[ -e "$WEIGHTS_DIR/prostt5-f16.gguf" ] && log "  prostt5:  $WEIGHTS_DIR" || log "  prostt5:  (not fetched)"
exit 0
