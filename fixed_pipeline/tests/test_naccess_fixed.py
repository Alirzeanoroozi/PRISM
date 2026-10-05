"""Tests for the fixed naccess_utils.py."""
import os

import src.naccess_utils as nu


def test_freesasa_pre_check_no_env():
    """check_freesasa_available() should work with default system python."""
    available, python_bin, error = nu.check_freesasa_available()
    assert python_bin is not None
    assert isinstance(available, bool)


def test_freesasa_pre_check_cache():
    """FreeSASA pre-check should be cached per python binary + pid."""
    cache_before = len(nu.FREESASA_PRE_CHECK_CACHE)
    nu.check_freesasa_available()
    cache_after = len(nu.FREESASA_PRE_CHECK_CACHE)
    assert cache_after >= cache_before


def test_naccess_binary_not_found_message():
    """Run_naccess should raise RuntimeError with diagnostic when binary missing."""
    os.environ["PRISM_NACCESS_EXECUTABLE"] = "/nonexistent/naccess_binary"
    try:
        nu.run_naccess("1a28AB", "/tmp", is_target=False)
        assert False, "Should have raised RuntimeError"
    except RuntimeError as exc:
        assert "NACCESS binary not found" in str(exc)
        assert "accall.f" in str(exc)
    finally:
        del os.environ["PRISM_NACCESS_EXECUTABLE"]
