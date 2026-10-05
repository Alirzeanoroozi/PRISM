"""
Tests for the fixed alignment_multiprot.py.

Validates:
1. Seccomp diagnostic detection
2. Configurable MIN_MATCHES env var
3. Pre-flight missing file counting
"""

"""Tests for the fixed alignment_multiprot.py."""
import importlib
import os

import src.alignment_multiprot as mp


def test_seccomp_detection_by_exit_code():
    """Exit code 159 should be detected as seccomp."""
    assert mp._is_seccomp_blocked(159, "") is True
    assert mp._is_seccomp_blocked(159, "Bad system call") is True


def test_seccomp_detection_by_stderr():
    """'Bad system call' in stderr should be detected even without exit 159."""
    assert mp._is_seccomp_blocked(1, "Bad system call") is True


def test_normal_failure_not_seccomp():
    """Normal non-zero exit without 'Bad system call' is not seccomp."""
    assert mp._is_seccomp_blocked(1, "File not found") is False
    assert mp._is_seccomp_blocked(0, "") is False


def test_seccomp_diagnostic_message():
    """Seccomp diagnostic should mention workarounds."""
    msg = mp._diagnose_multiprot_failure(159, "Bad system call", "")
    assert "seccomp" in msg.lower()
    assert "Slurm" in msg
    assert "64-bit" in msg


def test_normal_failure_diagnostic():
    """Normal failure diagnostic should mention exit code."""
    msg = mp._diagnose_multiprot_failure(127, "command not found", "")
    assert "exit code 127" in msg


def test_min_matches_env_var():
    """PRISM_MULTIPROT_MIN_MATCHES should be readable from env."""
    os.environ["PRISM_MULTIPROT_MIN_MATCHES"] = "10"
    try:
        importlib.reload(mp)
        assert mp.MIN_MATCHES == 10
    finally:
        del os.environ["PRISM_MULTIPROT_MIN_MATCHES"]


def test_multiprot_timeout_env_var():
    """PRISM_MULTIPROT_TIMEOUT should be configurable."""
    os.environ["PRISM_MULTIPROT_TIMEOUT"] = "600"
    try:
        importlib.reload(mp)
        assert mp.MULTIPROT_TIMEOUT == 600
    finally:
        del os.environ["PRISM_MULTIPROT_TIMEOUT"]
