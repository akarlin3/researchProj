import os
import sys

# Make the sextant package importable without installation.
_CORE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CORE not in sys.path:
    sys.path.insert(0, _CORE)

_EXTRACT = os.path.join(os.path.dirname(_CORE), "data", "osipi", "extracted")


def osipi_available() -> bool:
    return os.path.exists(os.path.join(_EXTRACT, "abdomen.nii.gz"))


def extract_dir() -> str:
    return _EXTRACT
