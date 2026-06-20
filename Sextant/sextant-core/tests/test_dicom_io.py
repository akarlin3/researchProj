"""The robust b-value extractor handles both the standard DICOM tag and the
Siemens SequenceName encoding used by TCGA-LIHC."""
from pydicom.dataset import Dataset

from sextant.dicom_io import bvalue


def test_bvalue_from_sequence_name():
    ds = Dataset()
    ds.SequenceName = "*ep_b500t"
    assert bvalue(ds) == 500.0
    ds.SequenceName = "*ep_b0"
    assert bvalue(ds) == 0.0


def test_bvalue_from_standard_tag_takes_priority():
    ds = Dataset()
    ds.add_new((0x0018, 0x9087), "FD", 600.0)
    assert bvalue(ds) == 600.0


def test_bvalue_absent_returns_none():
    assert bvalue(Dataset()) is None
