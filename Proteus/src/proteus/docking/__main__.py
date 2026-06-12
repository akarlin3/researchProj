"""Package entrypoint so `python -m proteus.docking` works.

`proteus.docking` is a package (the docking helpers live in __init__.py), and
`python -m <package>` executes the package's __main__ submodule — not __init__.
Without this file the documented invocation fails with "cannot be directly
executed". The CLI itself is `proteus.docking.main`."""
from __future__ import annotations

from proteus.docking import main

if __name__ == "__main__":
    raise SystemExit(main())
