import importlib

from importlib.metadata import PackageNotFoundError

import ymf


def test_version_fallback_when_not_installed(monkeypatch):
    def raise_not_found(name):
        raise PackageNotFoundError(name)

    monkeypatch.setattr("importlib.metadata.version", raise_not_found)
    importlib.reload(ymf)
    assert ymf.__version__ == "0.0.0+unknown"

    importlib.reload(ymf)
