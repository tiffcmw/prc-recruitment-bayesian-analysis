"""
Pins the production lattice used by run_gillespie.py so that any change is deliberate.

site_spacing = 0.2 (nm) differs from the 8 nm tubulin dimer spacing used in the CyLaKS
papers (Steckhahn et al. 2025 Table 1; Fiorenza et al. 2021). It is unresolved and logged
in ISSUES.md. If you change these values on purpose, update this test and ISSUES.md together,
because every fitted per-site rate depends on the number of sites.
"""
import pytest

import run_gillespie as rg


@pytest.fixture
def captured_states(monkeypatch):
    created = []

    class SpyState(rg.State):
        def __init__(self, *args, **kwargs):
            created.append((args, self))
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(rg, "State", SpyState)
    return created


def test_production_lattice_is_pinned(captured_states):
    rg.run_gillespie_prc1(4e-5, 1.0, 10.0, 0.1, end_time=0.01, max_steps=1)
    (args, state), = captured_states
    microtubule_length, site_spacing, microtubule_offset = args[0], args[1], args[2]

    assert microtubule_length == 5000.0   # nm (5 um, matches the papers)
    assert site_spacing == 0.2            # nm (papers use 8 nm; see ISSUES.md issue 2)
    assert microtubule_offset == 2000     # nm
    assert state.num_sites == 25001       # sites per lane
