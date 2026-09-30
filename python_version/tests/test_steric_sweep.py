"""
Tests for the sweep tooling (python_version/steric_sweep.py): the driver matches run_gillespie_prc1,
saved files have the documented fields, and results reload unchanged.
"""
import numpy as np
import pytest

import run_gillespie as rg
import steric_sweep as ss
from prc1_state import State


def tiny_config(**overrides):
    """A short microtubule so a whole sweep takes a fraction of a second."""
    cfg = ss.SweepConfig(
        ranges_nm=(0.0, 4.0), modes=("graded", "hard"), seeds=(0, 1), end_time=4.0, max_steps=2000,
        max_distance_nm=20.0,
    )
    cfg.lattice = dict(ss.PRODUCTION_LATTICE, microtubule_length=300.0, microtubule_offset=100)
    cfg.rates = dict(cfg.rates, initial_binding_rate_per_site=2e-3)
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def test_sweep_driver_matches_run_gillespie_prc1_on_the_production_lattice():
    # guards against the sweep's copy of the rate function and constants drifting from run_gillespie.py
    np.random.seed(5)
    expected, expected_times = rg.run_gillespie_prc1(4e-5, 1.0, 10.0, 0.1, end_time=6, max_steps=300)
    np.random.seed(5)
    state = State(**ss.PRODUCTION_LATTICE, **ss.DEFAULT_RATES)
    got, got_times = ss.simulate(state, end_time=6, max_steps=300)
    assert list(got) == list(expected)
    assert np.array_equal(got_times, expected_times)


def test_production_constants_match_run_gillespie(monkeypatch):
    created = []

    class Spy(rg.State):
        def __init__(self, *args, **kwargs):
            created.append(args)
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(rg, "State", Spy)
    rg.run_gillespie_prc1(4e-5, 1.0, 10.0, 0.1, end_time=0.01, max_steps=1)
    (length, spacing, offset, spring, rest, kbt, sep, *_rest), = created
    lat = ss.PRODUCTION_LATTICE
    assert (length, spacing, offset, spring, rest, kbt, sep) == (
        lat["microtubule_length"], lat["site_spacing"], lat["microtubule_offset"], lat["spring_constant"],
        lat["rest_length"], lat["k_B_T"], lat["microtubule_separation"])


def test_settings_have_a_single_off_baseline():
    cfg = ss.SweepConfig(ranges_nm=(0.0, 2.0, 8.0), modes=("graded", "hard"))
    assert cfg.settings() == [(0.0, "off"), (2.0, "graded"), (2.0, "hard"), (8.0, "graded"), (8.0, "hard")]


def test_run_one_is_reproducible_and_has_expected_shapes():
    cfg = tiny_config()
    a = ss.run_one(cfg, 4.0, "graded", seed=3)
    b = ss.run_one(cfg, 4.0, "graded", seed=3)
    for key in ("pair_counts", "pair_expected", "attach_nn_counts"):
        assert np.array_equal(a[key], b[key])
    assert a["n_final"] == b["n_final"]
    assert a["pair_counts"].shape == (cfg.max_distance_sites + 1,)
    assert a["pair_expected"][0] == 0.0
    assert a["attach_nn_counts"].sum() <= a["attach_events"]


def test_hard_exclusion_leaves_no_first_head_inside_the_footprint():
    cfg = tiny_config(seeds=(0, 1, 2), end_time=8.0)
    range_sites = int(round(4.0 / cfg.lattice["site_spacing"]))
    total_inside = 0
    for seed in cfg.seeds:
        out = ss.run_one(cfg, 4.0, "hard", seed)
        total_inside += out["attach_nn_counts"][1:range_sites + 1].sum()
    assert total_inside == 0


def test_binding_profile_shapes_for_each_mode():
    cfg = tiny_config()
    d, w, sites = ss.binding_profile(cfg, 4.0, "graded")
    k = sites
    assert sites == 20 and w[0] == 0.0
    assert np.all(np.diff(w[1:k + 1]) > 0)          # strictly increasing inside the range
    assert np.all(w[k + 1:] == 1.0)                  # no effect beyond it
    d, w, sites = ss.binding_profile(cfg, 4.0, "hard")
    assert np.all(w[:sites + 1] == 0.0) and np.all(w[sites + 1:] == 1.0)
    d, w, sites = ss.binding_profile(cfg, 0.0, "off")
    assert sites == 0 and w[0] == 0.0 and np.all(w[1:] == 1.0)


def test_saved_files_have_the_documented_fields_and_reload(tmp_path):
    cfg = tiny_config()
    paths = ss.run_sweep(cfg, tmp_path, verbose=False)
    assert len(paths) == len(cfg.settings()) == 3
    assert (tmp_path / "summary.csv").exists() and (tmp_path / "manifest.json").exists()

    results = ss.load_results(tmp_path)
    assert set(results) == {(0.0, "off"), (4.0, "graded"), (4.0, "hard")}
    for (range_nm, mode), item in results.items():
        for key in ("seed", "range_nm", "mode", "runtime_s", "n_final", "pair_counts", "pair_expected",
                    "attach_nn_counts", "profile_distance_nm", "profile_weight", "site_spacing"):
            assert key in item
        assert list(item["seed"]) == [0, 1]
        assert item["pair_counts"].shape[0] == 2
        assert item["mode"] == mode and item["range_nm"] == range_nm
        assert np.all(item["runtime_s"] > 0)


def test_load_results_reports_a_missing_sweep(tmp_path):
    with pytest.raises(FileNotFoundError):
        ss.load_results(tmp_path)


# ---------------------------------------------------------------- aggregation and plotting

def synthetic_item(seeds=3, max_sites=20, spacing=1.0):
    counts = np.zeros((seeds, max_sites + 1))
    expected = np.zeros((seeds, max_sites + 1))
    counts[:, 1:] = 4.0
    expected[:, 1:] = 2.0            # g = 2 everywhere
    attach = np.zeros((seeds, max_sites + 1), dtype=np.int64)
    attach[:, 1:6] = 1               # 5 of 10 attachments have a neighbor within 5 sites
    return dict(site_spacing=spacing, pair_counts=counts, pair_expected=expected, attach_nn_counts=attach,
                attach_events=np.full(seeds, 10), seed=np.arange(seeds), n_final=np.array([10, 12, 14]),
                runtime_s=np.array([1.0, 1.0, 1.0]))


def test_pair_correlation_bins_and_averages_correctly():
    x, mean, sem = ss.pair_correlation(synthetic_item(), bin_nm=4.0, max_nm=20.0)
    assert np.allclose(mean, 2.0) and np.allclose(sem, 0.0)
    assert np.allclose(x, [2.0, 6.0, 10.0, 14.0, 18.0])   # centers of 4 nm bins starting at 0


def test_attach_fraction_is_cumulative_and_bounded():
    d, mean, sem = ss.attach_fraction_within(synthetic_item(), max_nm=10.0)
    assert mean[0] == 0.0 and mean[5] == pytest.approx(0.5) and mean[-1] == pytest.approx(0.5)
    assert np.all(np.diff(mean) >= 0) and mean.max() <= 1.0


def test_error_bands_use_the_spread_across_seeds():
    item = synthetic_item()
    item["pair_counts"][:, 1:] = np.array([[2.0], [4.0], [6.0]])   # g = 1, 2, 3 across seeds
    x, mean, sem = ss.pair_correlation(item, bin_nm=4.0, max_nm=20.0)
    assert np.allclose(mean, 2.0) and np.allclose(sem, 1.0 / np.sqrt(3))


def test_range_colors_darken_with_range():
    ranges = [2.0, 4.0, 8.0, 16.0, 32.0]
    lum = [sum(ss._hex_to_rgb(ss.range_color(r, ranges))) for r in ranges]
    assert lum == sorted(lum, reverse=True) and len(set(lum)) == 5
    assert ss.range_color(2.0, ranges) == ss.RANGE_RAMP[0] and ss.range_color(32.0, ranges) == ss.RANGE_RAMP[-1]


def test_summarize_reports_one_row_per_setting(tmp_path):
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)
    rows = ss.summarize(ss.load_results(tmp_path))
    assert [(r["range_nm"], r["mode"]) for r in rows] == [(0.0, "off"), (4.0, "graded"), (4.0, "hard")]
    assert all(r["seeds"] == 2 and r["n_final_mean"] >= 0 for r in rows)


def test_figures_regenerate_from_saved_results_without_simulating(tmp_path, monkeypatch):
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)

    def no_simulation(*args, **kwargs):
        raise AssertionError("plotting must not run simulations")

    monkeypatch.setattr(ss, "run_one", no_simulation)
    monkeypatch.setattr(ss, "simulate", no_simulation)
    paths = ss.make_figures(tmp_path, tmp_path / "figs")
    assert [p.name for p in paths] == ["binding_profile.png", "pair_profile.png", "attach_profile.png", "summary.png"]
    assert all(p.exists() and p.stat().st_size > 2000 for p in paths)


def test_summary_panels_do_not_share_a_y_axis(tmp_path):
    # counts and seconds have different scales; sharing the axis flattened the runtime panel
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)
    fig = ss.plot_summary(ss.load_results(tmp_path))
    left, right = fig.axes[:2]
    assert not left.get_shared_y_axes().joined(left, right)


def test_command_line_plot_prints_a_summary(tmp_path, capsys):
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)
    ss.main(["plot", "--results", str(tmp_path), "--figures", str(tmp_path / "figs")])
    out = capsys.readouterr().out
    assert "range_nm" in out and "wrote" in out
