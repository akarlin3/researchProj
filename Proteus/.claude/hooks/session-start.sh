#!/usr/bin/env bash
# Proteus SessionStart hook — install the local toolchain so the test suite runs
# in Claude Code on the web (instead of skipping the tool-gated tests). Synchronous:
# the session waits until deps are ready (no race conditions). Web sessions only.
set -uo pipefail

# Only run in the remote (Claude Code on the web) environment.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

REPO="${CLAUDE_PROJECT_DIR:-$(pwd)}"
export PROTEUS_TOOLS_DIR="${PROTEUS_TOOLS_DIR:-$HOME/.proteus-tools}"
export PROTEUS_FETCH_WEIGHTS=1   # fetch ProstT5 weights so the S1/S2/pipeline tests run

bash "$REPO/scripts/setup_tools.sh"

# Persist tool PATH + env for the whole session.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  {
    echo "export PATH=\"$PROTEUS_TOOLS_DIR/bin:\$PATH\""
    echo "export PYTHONPATH=\"$REPO/src:\${PYTHONPATH:-}\""
    [ -e "$PROTEUS_TOOLS_DIR/prostt5/prostt5-f16.gguf" ] \
      && echo "export PROTEUS_PROSTT5_MODEL=\"$PROTEUS_TOOLS_DIR/prostt5\""
  } >> "$CLAUDE_ENV_FILE"
fi
