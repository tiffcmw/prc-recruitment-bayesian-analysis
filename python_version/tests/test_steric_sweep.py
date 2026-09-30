"""
Tests for the sweep tooling (python_version/steric_sweep.py): the driver matches run_gillespie_prc1,
saved files have the documented fields, and results reload unchanged.
"""
import warnings

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
                attach_events=np.full(seeds, 10), attach_with_ref=np.full(seeds, 10),
                attach_expected=np.zeros((seeds, max_sites + 1)),
                seed=np.arange(seeds), n_final=np.array([10, 12, 14]),
                runtime_s=np.array([1.0, 1.0, 1.0]))


def test_pair_correlation_bins_and_averages_correctly():
    x, mean, sem = ss.pair_correlation(synthetic_item(), bin_nm=4.0, max_nm=20.0)
    assert np.allclose(mean, 2.0) and np.allclose(sem, 0.0)
    # bins cover sites 1-4, 5-8, ...; the center of sites 1..4 is 2.5 sites (a bin from 0 would start at 0)
    assert np.allclose(x, [2.5, 6.5, 10.5, 14.5, 18.5])


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


# ---------------------------------------------------------------- fixes from the fresh-context review

class FakeState:
    """Just enough of State for the pair accumulator: sorted taken sites per lane."""
    def __init__(self, num_sites, bottom, top=()):
        self.num_sites = num_sites
        self.bottom_taken_sites = sorted(bottom)
        self.top_taken_sites = sorted(top)


def test_overlap_windows_follow_the_offset():
    lat = dict(ss.PRODUCTION_LATTICE)   # 5000 nm, 0.2 nm sites, offset 2000 nm
    assert ss.overlap_windows(lat) == {"bottom": (10000, 25000), "top": (0, 15000)}
    assert ss.overlap_windows(dict(lat, microtubule_offset=-2000)) == {"bottom": (0, 15000), "top": (10000, 25000)}
    whole = {"bottom": (0, 25000), "top": (0, 25000)}
    assert ss.overlap_windows(dict(lat, microtubule_offset=9000)) == whole   # no overlap: whole lane


def test_pair_expectation_is_one_for_heads_uniform_inside_the_overlap():
    # 40 heads uniform inside a 60% window: the whole-lane normalization gives about 1 / 0.6 (the
    # artifact the review found); the window normalization gives about 1.
    rng = np.random.default_rng(0)
    n_sites, lo, hi, max_d = 25001, 10000, 25000, 100
    windows = {"bottom": (lo, hi), "top": (0, hi - lo)}
    counts_w = np.zeros(max_d + 1); expected_w = np.zeros(max_d + 1)
    counts_l = np.zeros(max_d + 1); expected_l = np.zeros(max_d + 1)
    for _ in range(1500):
        heads = rng.choice(np.arange(lo, hi + 1), size=60, replace=False)
        state = FakeState(n_sites, heads)
        c, e, *_ = ss._accumulate_pairs(state, max_d, windows)
        counts_w += c; expected_w += e
        c, e, *_ = ss._accumulate_pairs(state, max_d, None)
        counts_l += c; expected_l += e
    assert counts_w[1:].sum() / expected_w[1:].sum() == pytest.approx(1.0, rel=0.05)
    assert counts_l[1:].sum() / expected_l[1:].sum() == pytest.approx(1 / 0.6, rel=0.05)


def test_pair_accumulator_reports_the_fraction_of_heads_in_the_window():
    state = FakeState(1000, [5, 500, 600, 700], [10, 20])
    windows = {"bottom": (400, 999), "top": (0, 599)}
    *_, in_window, total = ss._accumulate_pairs(state, 50, windows)
    assert (in_window, total) == (5, 6)


def test_uniform_binding_cdf_matches_brute_force():
    heads, n_sites, max_d = [3, 10], 15, 6
    free = [s for s in range(n_sites) if s not in heads]
    nearest = np.array([min(abs(s - h) for h in heads) for s in free])
    brute = np.array([(nearest <= d).sum() / len(free) for d in range(max_d + 1)])
    assert np.allclose(ss._uniform_binding_cdf(heads, n_sites, max_d), brute)


def test_attach_enrichment_is_one_for_uniform_binding_and_zero_inside_a_hard_footprint():
    rng = np.random.default_rng(3)
    n_sites, max_d, draws = 4000, 60, 60000
    heads = np.sort(rng.choice(n_sites, size=40, replace=False))
    free = np.setdiff1d(np.arange(n_sites), heads)
    nearest = np.abs(free[:, None] - heads[None, :]).min(axis=1)
    cdf = ss._uniform_binding_cdf(list(heads), n_sites, max_d)

    def item_from(sites):
        near = np.abs(sites[:, None] - heads[None, :]).min(axis=1)
        hist = np.bincount(near[near <= max_d], minlength=max_d + 1)[None, :]
        return dict(site_spacing=1.0, attach_nn_counts=hist, attach_with_ref=np.array([len(sites)]),
                    attach_expected=(len(sites) * cdf)[None, :])

    uniform = item_from(rng.choice(free, size=draws))
    d, mean, sem = ss.attach_enrichment(uniform, max_nm=40.0)
    assert np.allclose(mean[10:], 1.0, atol=0.06)
    hard = item_from(rng.choice(free[nearest > 8], size=draws))
    d, mean, sem = ss.attach_enrichment(hard, max_nm=40.0)
    assert np.all(mean[1:9] == 0.0) and mean[40] < 1.0


def test_mean_sem_of_a_single_seed_is_nan_not_zero():
    mean, sem = ss._mean_sem(np.array([[1.0, 2.0]]))
    assert np.allclose(mean, [1.0, 2.0]) and np.all(np.isnan(sem))


def test_mean_sem_with_sparse_columns_does_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean, sem = ss._mean_sem(np.array([[1.0, np.nan, np.nan], [3.0, np.nan, np.nan], [5.0, 2.0, np.nan]]))
    assert sem[0] == pytest.approx(2.0 / np.sqrt(3)) and np.isnan(sem[1]) and np.isnan(sem[2])


def test_attach_denominator_excludes_attachments_onto_an_empty_lane():
    cfg = tiny_config()
    np.random.seed(1)
    state = ss.RecordingState(**cfg.lattice, **cfg.rates, max_distance_sites=cfg.max_distance_sites)
    for _ in range(30):
        state.single_attach_prc1()
    empty_lane_attaches = state.attach_events - state.attach_with_ref
    assert state.attach_events == 30 and empty_lane_attaches in (1, 2)   # at most one per lane
    assert state.attach_nn_counts.sum() <= state.attach_with_ref


def test_a_fraction_of_one_is_reachable():
    item = synthetic_item(seeds=2, max_sites=10)
    item["attach_nn_counts"][:] = 0
    item["attach_nn_counts"][:, 2] = 10          # every attachment lands 2 sites from a bound head
    d, mean, sem = ss.attach_fraction_within(item, max_nm=10.0)
    assert mean[-1] == pytest.approx(1.0) and mean[1] == 0.0


def test_stale_files_from_an_earlier_sweep_do_not_mix_into_a_new_one(tmp_path):
    ss.run_sweep(tiny_config(ranges_nm=(0.0, 4.0)), tmp_path, verbose=False)
    ss.run_sweep(tiny_config(ranges_nm=(0.0, 6.0)), tmp_path, verbose=False)
    assert set(ss.load_results(tmp_path)) == {(0.0, "off"), (6.0, "graded"), (6.0, "hard")}
    assert not list(tmp_path.glob("range_4*"))


def test_files_not_listed_in_the_manifest_are_ignored(tmp_path):
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)
    stray = next(tmp_path.glob("range_*_graded.npz"))
    with np.load(stray) as data:
        fields = {k: data[k] for k in data.files}
    fields["range_nm"] = np.array(99.0)   # a file from some other sweep that reports a different range
    np.savez(tmp_path / "range_99nm_graded.npz", **fields)
    results = ss.load_results(tmp_path)
    assert (99.0, "graded") not in results and len(results) == 3


def test_small_ranges_get_distinct_file_names():
    assert ss.setting_name(0.001, "graded") != ss.setting_name(0.004, "graded")
    assert ss.setting_name(8.0, "hard") == "range_8nm_hard"


def test_a_truncated_run_is_flagged_and_gives_no_pair_statistics():
    cfg = tiny_config(max_steps=20, end_time=50.0)
    with pytest.warns(UserWarning, match="truncated"):
        out = ss.run_one(cfg, 4.0, "graded", seed=0)
    assert out["truncated"] and out["end_time_reached"] < cfg.end_time
    assert out["pair_counts"].sum() == 0 and out["pair_expected"].sum() == 0


def test_a_normal_run_is_not_truncated():
    out = ss.run_one(tiny_config(), 4.0, "graded", seed=0)
    assert not out["truncated"] and 0.0 <= out["window_fraction"] <= 1.0


def test_config_rejects_a_range_below_one_site_instead_of_running_it_as_steric_off():
    with pytest.raises(ValueError, match="below one site"):
        ss.SweepConfig(ranges_nm=(0.0, 0.1)).validate()
    ss.SweepConfig(ranges_nm=(0.0, 0.2)).validate()   # exactly one site is fine


def test_command_line_rejects_empty_or_unresolvable_settings(tmp_path):
    with pytest.raises(ValueError, match="below one site"):
        ss.main(["run", "--out", str(tmp_path), "--ranges", "0", "0.1", "--seeds", "1"])
    with pytest.raises(ValueError, match="seeds"):
        ss.main(["run", "--out", str(tmp_path), "--ranges", "0", "4", "--seeds", "0"])


def test_summary_includes_the_density_controlled_enrichment(tmp_path):
    ss.run_sweep(tiny_config(), tmp_path, verbose=False)
    rows = ss.summarize(ss.load_results(tmp_path))
    assert all("attach_enrichment_8nm_mean" in r for r in rows)


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
