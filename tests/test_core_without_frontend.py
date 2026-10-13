"""The core imports and works with no front end installed.

``dependencies`` is pyyaml alone; strictyaml, pint and sympy live behind the
``spec``, ``units``, ``symbolic`` and ``frontend`` extras. A solver installs
the core: Proteus writes its archive of record through ``ymf.archive`` on
every run and has no use for schema validation.

Nothing in ymf's own test environment enforces that, because ``dev`` installs
all three extras -- so a core module could grow a front-end import and no test
here would notice. One did. ``_is_archive_text`` answered "archive or
specification?", the question that decides *which loader to use*, and it lived
in ``ymf.compose``, which imports strictyaml at module level. Proteus's
``ymf_run`` imported it and died on "No module named 'strictyaml'"
(cekees/proteus#46) -- in a core-only environment, asking the question cost
the whole validating stack.

These tests run in a subprocess with the three blocked, which is the only way
to see the boundary from inside an environment that has them.
"""
import subprocess
import sys

import pytest

#: Importable with pyyaml and nothing else. Proteus imports from every one.
CORE = ["ymf", "ymf.archive", "ymf.xdmf", "ymf.cli", "ymf.closure",
        "ymf.domain", "ymf.references"]

#: Front end, and honestly so. Asserted as well as CORE: if the blocker below
#: silently stopped blocking, every CORE assertion would pass for the wrong
#: reason, and this is what catches that.
FRONTEND = {"ymf.schema": "strictyaml", "ymf.compose": "strictyaml",
            "ymf.units": "pint", "ymf.symbolic": "sympy"}

_BLOCKER = '''
import sys
BLOCKED = {"strictyaml", "pint", "sympy"}
class Blocker:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in BLOCKED:
            raise ModuleNotFoundError("No module named %r" % fullname,
                                      name=fullname)
        return None
sys.meta_path.insert(0, Blocker())
'''


def _without_frontend(code):
    """Run code in a subprocess where strictyaml, pint and sympy are absent."""
    return subprocess.run([sys.executable, "-c", _BLOCKER + code],
                          capture_output=True, text=True)


@pytest.mark.parametrize("module", CORE)
def test_a_core_module_needs_no_front_end(module):
    out = _without_frontend("import %s" % module)
    assert out.returncode == 0, out.stderr


@pytest.mark.parametrize("module,missing", sorted(FRONTEND.items()))
def test_a_front_end_module_says_what_it_needs(module, missing):
    out = _without_frontend("import %s" % module)
    assert out.returncode != 0
    assert missing in out.stderr


def test_which_loader_a_file_needs_is_a_core_question():
    """An archive goes to the fast loader, a specification to the strict one.

    Asking must not itself require the strict one.
    """
    out = _without_frontend('''
from ymf.archive import is_archive_text
assert is_archive_text("ymf: 0.2.0\\noutputs: {}\\n") is True
assert is_archive_text("outputs:\\n  a: {}\\n") is True
assert is_archive_text("model:\\n  name: x\\n") is False
print("ok")
''')
    assert out.returncode == 0, out.stderr
    assert "ok" in out.stdout


def test_compose_still_offers_the_old_private_name():
    """Proteus imports ymf.compose._is_archive_text; it must keep working."""
    from ymf.archive import is_archive_text
    from ymf.compose import _is_archive_text
    assert _is_archive_text is is_archive_text
