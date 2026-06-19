"""Make the in-tree ``vernier`` package importable when running pytest from the
Vernier root without an editable install. pytest inserts the rootdir (the
directory holding this conftest) onto sys.path in its default prepend import
mode, so ``import vernier`` resolves.
"""
