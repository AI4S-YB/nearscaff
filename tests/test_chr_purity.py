"""Unit tests for nearscaff chromosome-purity enforcement."""
from nearscaff.config import NearscaffConfig


def test_scaffold_config_chr_purity_defaults():
    cfg = NearscaffConfig()
    assert cfg.scaffold.enforce_chr_purity is True
    assert cfg.scaffold.min_chr_share == 0.20
    assert cfg.scaffold.min_chr_len == 1_000_000
