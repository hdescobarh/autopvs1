# Bundled dependencies

These private copies ship with AutoPVS1 and are imported through
`autopvs1._vendor`. They are not separate installed packages.

## maxentpy

Upstream: https://github.com/kepbod/maxentpy

Provides MaxEntScan splice-site scoring. AutoPVS1 uses `maxent.py` and its
text scoring matrices. The optional accelerator (`maxent_fast.py`, `_hashseq.c`,
`_hashseq.pyx`, and binary matrices) is retained as source material, but is not
built or supported by this package installation. Neither Cython nor msgpack is
required for normal use. The accelerator's import was made package-relative
when this copy moved into `_vendor`.

## pyhgvs

Upstream: https://github.com/counsyl/hgvs

The original AutoPVS1 README describes this copy as modified for Python 3 and
additional functionality. Preserve those changes; this is not an unmodified
upstream release.

The exact upstream revisions and complete prior patch histories of both copies
are unknown. Existing source notices and attribution are retained. This layout
change does not update their algorithms or replace either copy from PyPI.
