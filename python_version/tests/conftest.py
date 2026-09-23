import importlib.util
import sys
from pathlib import Path

import pytest

PYTHON_VERSION_DIR = Path(__file__).resolve().parents[1]

# prc1_state.py / prc1.py / gillespie.py / run_gillespie.py import each other by
# bare module name (e.g. "from prc1 import Prc1"), so python_version/ must be on
# sys.path for those imports to resolve regardless of where pytest is invoked from.
if str(PYTHON_VERSION_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_VERSION_DIR))


def _load_smc_abc():
    """
    smc-abc.py has a hyphen in its filename, so it isn't importable as a normal
    module (`import smc-abc` is a syntax error). Load it directly from its path.
    """
    spec = importlib.util.spec_from_file_location("smc_abc", PYTHON_VERSION_DIR / "smc-abc.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def smc_abc():
    return _load_smc_abc()
