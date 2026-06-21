"""Shared test setup."""
import os

# macOS libomp duplicate-initialization guard (torch + numpy/scipy both ship libomp).
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
