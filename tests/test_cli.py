"""Tests for nearscaff.cli._build_config wiring."""
from argparse import Namespace

from nearscaff.cli import _build_config


def test_chr_purity_defaults_when_flags_absent():
    cfg = _build_config(Namespace(threads=4, output="o"))
    assert cfg.scaffold.enforce_chr_purity is True
    assert cfg.scaffold.min_chr_share == 0.20
    assert cfg.scaffold.min_chr_len == 1_000_000


def test_no_chr_purity_flag_disables_and_overrides_thresholds():
    cfg = _build_config(Namespace(
        threads=4, output="o",
        no_chr_purity=True, min_chr_share=0.30, min_chr_len=2_000_000,
    ))
    assert cfg.scaffold.enforce_chr_purity is False
    assert cfg.scaffold.min_chr_share == 0.30
    assert cfg.scaffold.min_chr_len == 2_000_000
