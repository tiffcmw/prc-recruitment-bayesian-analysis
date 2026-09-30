"""
Steric blocking sweep: run simulations, save one .npz per setting, load and plot them.

Run (from the repo root):
    python python_version/steric_sweep.py run  [--out results/steric_sweep] [options]
    python python_version/steric_sweep.py plot [--results results/steric_sweep]

`run` writes results/steric_sweep/<setting>.npz plus summary.csv and manifest.json.
`plot` only reads those files, so figures can be regenerated without rerunning simulations.
The results folder is gitignored.

What is measured per seed. The same seed list is used for every setting, but these are not common random
numbers: the steric sampler draws different random numbers from the default one, so runs diverge at the
first attachment and settings are effectively independent samples.
  pair_counts / pair_expected  head pairs on the same lane at distance d sites, and the count expected if
                               heads were placed uniformly at random inside the overlap of the two
                               microtubules (where almost all heads are; see overlap_windows). g(d) =
                               counts / expected, so g = 1 means no correlation. Pooled over late snapshots.
  attach_nn_counts             distance (sites) from each newly attached first head to the nearest
                               same-lane head that was already bound.
  attach_with_ref              number of first-head attachments that had a bound head on the same lane
                               (the denominator; attachments onto an empty lane have no reference head).
  attach_expected              the same distances if each attachment had picked a free site uniformly
                               (no steric effect), given the occupancy at that moment. Observed over
                               expected therefore removes the dependence on density (attach_enrichment).
  n_final, runtime_s           PRC1 count at the end of the run and wall time of the simulation.
  truncated                    the run hit max_steps before end_time; its pair statistics are left empty.
Per setting (deterministic, from the model):
  profile_distance_nm / profile_weight   relative binding rate of a free site at distance d from one
                               occupied site (the binding-rate profile).

Steric blocking currently covers first-head binding only, so pair statistics also include second heads,
which are not steric aware (see ISSUES.md, issue 3).
"""
import argparse
import bisect
import csv
import json
import math
import os
import subprocess
import time
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from gillespie import run_gillespie
from prc1 import Prc1
from prc1_state import State, range_nm_to_sites

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "results" / "steric_sweep"

# Same constants as run_gillespie.run_gillespie_prc1 (guarded by a test that compares the two drivers).
PRODUCTION_LATTICE = dict(
    microtubule_length=5000.0, site_spacing=0.2, microtubule_offset=2000, spring_constant=2,
    rest_length=32, k_B_T=4.1, microtubule_separation=32, base_hopping_rate=1,
)
DEFAULT_RATES = dict(
    initial_binding_rate_per_site=4e-5, singly_bound_detachment_rate=1.0,
    base_double_attachment_rate=10.0, base_double_detachment_rate=0.1,
)


@dataclass
class SweepConfig:
    ranges_nm: tuple = (0.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    modes: tuple = ("graded", "hard")
    seeds: tuple = tuple(range(30))
    contact_factor: float = 0.0
    end_time: float = 30.0
    max_steps: int = 40000
    max_distance_nm: float = 64.0
    snapshot_fractions: tuple = (0.6, 0.7, 0.8, 0.9, 1.0)
    lattice: dict = field(default_factory=lambda: dict(PRODUCTION_LATTICE))
    # 2e-4 rather than the 4e-5 used in the tests: about 175 PRC1 per run instead of about 30, so pair
    # statistics have enough counts. A density choice for the sweep, not a fitted value.
    rates: dict = field(default_factory=lambda: dict(DEFAULT_RATES, initial_binding_rate_per_site=2e-4))

    @property
    def max_distance_sites(self):
        return int(round(self.max_distance_nm / self.lattice["site_spacing"]))

    def settings(self):
        """(range_nm, mode) pairs; range 0 is a single 'off' baseline."""
        out = [(0.0, "off")] if 0.0 in self.ranges_nm else []
        out += [(r, m) for r in self.ranges_nm if r > 0 for m in self.modes]
        return out

    def validate(self):
        """Refuse configurations that would silently run something other than what was asked for."""
        if not self.seeds:
            raise ValueError("seeds must not be empty")
        if not self.settings():
            raise ValueError("no settings to run: give at least one range")
        spacing = self.lattice["site_spacing"]
        for range_nm in self.ranges_nm:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sites = range_nm_to_sites(range_nm, spacing)   # also rejects negative and non-finite
            if range_nm > 0 and sites < 1:
                raise ValueError(f"range {range_nm} nm is below one site ({spacing} nm); it would run as steric "
                                 "off under a steric label. Use 0 for off or a range of at least one site.")


def setting_name(range_nm, mode):
    return f"range_{range_nm:g}nm_{mode}"


def overlap_windows(lattice):
    """
    Inclusive site ranges, in each lane's own index, where the two microtubules overlap. Nearly all heads
    sit there (doubly attached PRC1 need a partner on the other lane), so pair statistics are normalized
    to this window instead of the whole lane. Without an overlap the whole lane is used.
    """
    spacing = lattice["site_spacing"]
    n = int(lattice["microtubule_length"] / spacing) + 1
    offset = int(round(lattice["microtubule_offset"] / spacing))   # top lane origin, in bottom lane sites
    if abs(offset) >= n:
        return {"bottom": (0, n - 1), "top": (0, n - 1)}
    if offset >= 0:
        return {"bottom": (offset, n - 1), "top": (0, n - 1 - offset)}
    return {"bottom": (0, n - 1 + offset), "top": (-offset, n - 1)}


def _uniform_binding_cdf(heads, num_sites, max_d):
    """
    P(nearest bound head is within d sites) for d = 0..max_d, if a free site of this lane were chosen
    uniformly (no steric effect). Exact, from the gaps between the sorted heads. It is the null the
    attachment profile is compared with, and it depends on the occupancy at the moment of attachment.
    """
    h = np.asarray(heads)
    d = np.arange(max_d + 1)
    within = np.minimum(h[0], d) + np.minimum(num_sites - 1 - h[-1], d)          # the two ends of the lane
    if len(h) > 1:
        within = within + np.minimum((np.diff(h) - 1)[:, None], 2 * d[None, :]).sum(axis=0)   # gaps between heads
    return within / (num_sites - len(h))


# ------------------------------------------------------------------ simulation

class RecordingState(State):
    """
    State that records, for every first-head attachment onto a lane that already has a bound head, the
    distance to the nearest bound head, and what that distance would be if the free site were chosen uniformly.
    """

    def __init__(self, *args, max_distance_sites=0, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_distance_sites = max_distance_sites
        self.attach_nn_counts = np.zeros(max_distance_sites + 1, dtype=np.int64)
        self.attach_expected = np.zeros(max_distance_sites + 1)
        self.attach_events = 0
        self.attach_with_ref = 0

    def single_attach_prc1(self):
        before = (list(self.bottom_taken_sites), list(self.top_taken_sites))
        super().single_attach_prc1()
        self.attach_events += 1
        for lane_before, taken in zip(before, (self.bottom_taken_sites, self.top_taken_sites)):
            if len(taken) == len(lane_before) or not lane_before:
                continue   # the other lane, or the first head on an empty lane (no reference head)
            self.attach_with_ref += 1
            self.attach_expected += _uniform_binding_cdf(lane_before, self.num_sites, self.max_distance_sites)
            (new_site,) = set(taken).difference(lane_before)
            i = bisect.bisect_left(lane_before, new_site)
            nearest = min(abs(new_site - lane_before[j]) for j in (i - 1, i) if 0 <= j < len(lane_before))
            if nearest <= self.max_distance_sites:
                self.attach_nn_counts[nearest] += 1


def simulate(state, end_time, max_steps, statistic_function=None):
    """
    Single attach, double attach and detach (no hopping) on `state`, mirroring
    run_gillespie.run_gillespie_prc1 with enable_hopping=False. Returns (statistics, times); the default
    statistic is the number of PRC1 at each step.
    """
    def rate_function(s):
        initial_binding_rate = s.initial_binding_rate
        if len(s) == 0:
            return [initial_binding_rate, 0, 0, 0, 0]
        rates = [initial_binding_rate] + [0] * (len(s) - 1)
        for prc1 in s:
            rates.append(prc1.double_attachment_rate)
        for prc1 in s:
            rates.append(prc1.detachment_rate)
        return rates

    reactions = [
        lambda s, index: s.single_attach_prc1(),
        lambda s, index: s.double_attach_prc1(index),
        lambda s, index: s.detach_prc1(index),
    ]

    def default_statistic(s, stat_list, time_):
        stat_list.append(len(s))

    return run_gillespie(state, end_time, rate_function, reactions,
                         statistic_function or default_statistic, lambda s, t, dt: t + dt,
                         np.inf, max_steps, lambda s, t, dt: None)


def _accumulate_pairs(state, max_d, windows=None):
    """
    Same-lane head pairs at distance 1..max_d sites among heads inside `windows` (per-lane inclusive site
    ranges; the whole lane if None), and the pairs expected if those heads were uniform inside the window.
    Returns (counts, expected, heads_in_window, heads_total).
    """
    counts = np.zeros(max_d + 1)
    expected = np.zeros(max_d + 1)
    distances = np.arange(max_d + 1)
    in_window = total = 0
    for lane, taken in (("bottom", state.bottom_taken_sites), ("top", state.top_taken_sites)):
        lo, hi = windows[lane] if windows else (0, state.num_sites - 1)
        sites = [s for s in taken if lo <= s <= hi]
        total += len(taken)
        in_window += len(sites)
        n, width = len(sites), hi - lo + 1
        if n < 2 or width < 2:
            continue
        for i, site in enumerate(sites):
            j = i + 1
            while j < n and sites[j] - site <= max_d:
                counts[sites[j] - site] += 1
                j += 1
        expected += (n * (n - 1) / 2) * np.maximum(width - distances, 0) / (width * (width - 1) / 2)
    return counts, expected, in_window, total


def run_one(cfg, range_nm, mode, seed):
    """One seeded simulation; returns the measurements listed in the module docstring."""
    steric_mode = "graded" if mode == "off" else mode
    max_d = min(cfg.max_distance_sites, int(cfg.lattice["microtubule_length"] / cfg.lattice["site_spacing"]))
    state = RecordingState(
        **cfg.lattice, **cfg.rates, max_distance_sites=max_d,
        steric_range_nm=range_nm, steric_mode=steric_mode, steric_contact_factor=cfg.contact_factor,
    )
    windows = overlap_windows(cfg.lattice)
    counts = np.zeros(max_d + 1)
    expected = np.zeros(max_d + 1)
    heads_in_window = heads_total = 0
    snapshot_times = [f * cfg.end_time for f in cfg.snapshot_fractions]
    taken_snapshots = 0

    def statistic(s, stat_list, time_):
        nonlocal taken_snapshots, heads_in_window, heads_total
        stat_list.append(len(s))
        while taken_snapshots < len(snapshot_times) and time_ >= snapshot_times[taken_snapshots]:
            c, e, in_window, total = _accumulate_pairs(s, max_d, windows)
            counts[:] += c
            expected[:] += e
            heads_in_window += in_window
            heads_total += total
            taken_snapshots += 1

    np.random.seed(seed)
    t0 = time.perf_counter()
    statistics, times = simulate(state, cfg.end_time, cfg.max_steps, statistic)
    runtime = time.perf_counter() - t0
    truncated = bool(times[-1] < cfg.end_time)
    if truncated:
        # the run stopped in the transient, before the snapshot times: report it and keep no pair statistics
        warnings.warn(f"run truncated at max_steps={cfg.max_steps} (reached t={times[-1]:.2f} of {cfg.end_time}); "
                      "pair statistics are left empty for this run", UserWarning, stacklevel=2)
        counts[:] = 0.0
        expected[:] = 0.0
        heads_in_window = heads_total = 0
    expected[0] = 0.0
    return dict(
        pair_counts=counts, pair_expected=expected, attach_nn_counts=state.attach_nn_counts.copy(),
        attach_expected=state.attach_expected.copy(), attach_with_ref=state.attach_with_ref,
        n_final=len(state), attach_events=state.attach_events, runtime_s=runtime, steps=len(times),
        end_time_reached=float(times[-1]), truncated=truncated,
        window_fraction=heads_in_window / heads_total if heads_total else float("nan"),
    )


def binding_profile(cfg, range_nm, mode):
    """
    Relative binding rate of a free site at distance d from one occupied site (1 = unaffected), from
    the model itself. Returns (distance_nm, weight, range_sites).
    """
    spacing = cfg.lattice["site_spacing"]
    span_sites = int(round((max(range_nm, 8.0) + 8.0) / spacing))
    if mode == "off":
        weights = np.ones(span_sites + 1)
        weights[0] = 0.0
        return np.arange(span_sites + 1) * spacing, weights, 0
    length = (2 * span_sites + 4) * spacing
    state = State(**{**cfg.lattice, "microtubule_length": length, "microtubule_offset": 0},
                  **cfg.rates, steric_range_nm=range_nm, steric_mode=mode,
                  steric_contact_factor=cfg.contact_factor)
    center = span_sites + 2
    head = Prc1(state)
    head.binding_site_bottom = center
    head.set_closest_neighbors()
    weights = np.array([state.steric_weight("bottom", center + k) for k in range(span_sites + 1)])
    return np.arange(span_sites + 1) * spacing, weights, state.steric_range_sites


def run_setting(cfg, range_nm, mode, out_dir, verbose=True):
    runs = [run_one(cfg, range_nm, mode, seed) for seed in cfg.seeds]
    distance_nm, weight, range_sites = binding_profile(cfg, range_nm, mode)
    path = Path(out_dir) / f"{setting_name(range_nm, mode)}.npz"
    partial = path.with_name(path.stem + ".partial.npz")   # written fully, then renamed, so no half-written file
    windows = overlap_windows(cfg.lattice)
    np.savez_compressed(
        partial,
        seed=np.array(cfg.seeds), range_nm=range_nm, mode=mode, range_sites=range_sites,
        truncated=np.array([r["truncated"] for r in runs]),
        window_fraction=np.array([r["window_fraction"] for r in runs]),
        pair_window_bottom=np.array(windows["bottom"]), pair_window_top=np.array(windows["top"]),
        attach_with_ref=np.array([r["attach_with_ref"] for r in runs]),
        attach_expected=np.array([r["attach_expected"] for r in runs]),
        contact_factor=cfg.contact_factor, site_spacing=cfg.lattice["site_spacing"],
        rate_per_site=cfg.rates["initial_binding_rate_per_site"],
        runtime_s=np.array([r["runtime_s"] for r in runs]),
        n_final=np.array([r["n_final"] for r in runs]),
        steps=np.array([r["steps"] for r in runs]),
        attach_events=np.array([r["attach_events"] for r in runs]),
        end_time_reached=np.array([r["end_time_reached"] for r in runs]),
        pair_counts=np.array([r["pair_counts"] for r in runs]),
        pair_expected=np.array([r["pair_expected"] for r in runs]),
        attach_nn_counts=np.array([r["attach_nn_counts"] for r in runs]),
        profile_distance_nm=distance_nm, profile_weight=weight,
        config_json=json.dumps(asdict(cfg), default=list),
    )
    os.replace(partial, path)
    if verbose:
        mean_n = np.mean([r["n_final"] for r in runs])
        total = sum(r["runtime_s"] for r in runs)
        print(f"  {path.name}: {len(runs)} seeds, mean final n = {mean_n:.1f}, {total:.1f} s", flush=True)
    return path, runs


def run_sweep(cfg, out_dir=DEFAULT_OUT, verbose=True):
    """Run every setting and save it. Files from an earlier sweep in `out_dir` are removed first."""
    cfg.validate()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in [*out_dir.glob("range_*.npz"), out_dir / "summary.csv", out_dir / "manifest.json"]:
        stale.unlink(missing_ok=True)
    rows, paths = [], []
    for range_nm, mode in cfg.settings():
        path, runs = run_setting(cfg, range_nm, mode, out_dir, verbose)
        paths.append(path)
        rows += [dict(setting=path.stem, range_nm=range_nm, mode=mode, seed=seed, n_final=r["n_final"],
                      attach_events=r["attach_events"], steps=r["steps"], runtime_s=round(r["runtime_s"], 4))
                 for seed, r in zip(cfg.seeds, runs)]
    with open(out_dir / "summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                cwd=REPO_ROOT, check=True).stdout.strip()
    except Exception:
        commit = None
    (out_dir / "manifest.json").write_text(json.dumps(
        dict(config=asdict(cfg), git_commit=commit, files=[p.name for p in paths]), indent=2, default=list))
    return paths


# ------------------------------------------------------------------ loading

def load_results(results_dir=DEFAULT_OUT):
    """
    Load the settings listed in manifest.json (written last, so an interrupted sweep has none), and only
    those, so files left over from another sweep are never mixed in.
    Returns {(range_nm, mode): dict of arrays and scalars}.
    """
    results_dir = Path(results_dir)
    manifest = results_dir / "manifest.json"
    if not manifest.exists():
        raise FileNotFoundError(f"no manifest.json in {results_dir}; run the sweep first "
                                "(an interrupted sweep leaves no manifest)")
    results = {}
    for name in json.loads(manifest.read_text())["files"]:
        with np.load(results_dir / name, allow_pickle=False) as data:
            item = {k: data[k] for k in data.files}
        for key in ("range_nm", "site_spacing", "rate_per_site", "contact_factor"):
            item[key] = float(item[key])
        item["mode"] = str(item["mode"])
        item["range_sites"] = int(item["range_sites"])
        results[(item["range_nm"], item["mode"])] = item
    return results


# ------------------------------------------------------------------ aggregation (reads loaded results only)

def _mean_sem(values):
    """
    Mean and standard error across seeds (axis 0), ignoring NaN. The standard error is NaN, not 0,
    where fewer than two seeds have a value: one seed says nothing about the spread.
    """
    values = np.asarray(values, dtype=float)
    n = np.sum(~np.isnan(values), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # empty or single-value columns are handled below
        mean = np.nanmean(values, axis=0)
        std = np.nanstd(values, axis=0, ddof=1)
    sem = np.where(n >= 2, std / np.sqrt(np.maximum(n, 1)), np.nan)
    return mean, sem


def pair_correlation(item, bin_nm=4.0, max_nm=64.0):
    """
    g(d) binned in nm: same-lane head pairs divided by the pairs expected under uniform placement inside
    the overlap window, computed per seed, then mean and standard error across seeds (seeds with no
    expected pairs in a bin, such as truncated runs, are left out). Returns (centers_nm, mean, sem).
    """
    spacing = item["site_spacing"]
    bin_sites = max(1, int(round(bin_nm / spacing)))
    max_sites = min(int(round(max_nm / spacing)), item["pair_counts"].shape[1] - 1)
    starts = np.arange(1, max_sites + 1, bin_sites)
    ends = np.minimum(starts + bin_sites - 1, max_sites)
    per_seed = []
    for counts, expected in zip(item["pair_counts"], item["pair_expected"]):
        c = np.array([counts[s:e + 1].sum() for s, e in zip(starts, ends)])
        ex = np.array([expected[s:e + 1].sum() for s, e in zip(starts, ends)])
        per_seed.append(np.where(ex > 0, c / np.where(ex > 0, ex, 1), np.nan))
    mean, sem = _mean_sem(per_seed)
    return (starts + ends) / 2.0 * spacing, mean, sem


def attach_fraction_within(item, max_nm=40.0):
    """
    Fraction of first-head attachments, among those with a bound head on the same lane, that land within d nm
    of the nearest bound head, per seed, then mean and standard error. It depends on how crowded the lane
    is, so compare settings through attach_enrichment. Returns (distance_nm, mean, sem).
    """
    spacing = item["site_spacing"]
    max_sites = min(int(round(max_nm / spacing)), item["attach_nn_counts"].shape[1] - 1)
    cumulative = np.cumsum(item["attach_nn_counts"][:, :max_sites + 1], axis=1)
    fractions = cumulative / np.maximum(item["attach_with_ref"][:, None], 1)
    mean, sem = _mean_sem(fractions)
    return np.arange(max_sites + 1) * spacing, mean, sem


def attach_enrichment(item, max_nm=40.0):
    """
    Observed over expected number of first-head attachments within d nm of a bound head, where the
    expectation assumes each attachment picked a free site uniformly (no steric effect) given the
    occupancy at that moment. 1 means no steric effect and 0 means none land that close, whatever the
    density. Per seed, then mean and standard error; a seed is left out at a distance where fewer than one
    attachment was expected. Returns (distance_nm, mean, sem).
    """
    spacing = item["site_spacing"]
    max_sites = min(int(round(max_nm / spacing)), item["attach_nn_counts"].shape[1] - 1)
    observed = np.cumsum(item["attach_nn_counts"][:, :max_sites + 1], axis=1)
    expected = item["attach_expected"][:, :max_sites + 1]
    ratio = np.where(expected >= 1.0, observed / np.where(expected >= 1.0, expected, 1.0), np.nan)
    mean, sem = _mean_sem(ratio)
    return np.arange(max_sites + 1) * spacing, mean, sem


def summarize(results, near_nm=(8.0, 16.0)):
    """One row per setting: final count, runtime, attachments, and how first heads land near a bound head."""
    rows = []
    for (range_nm, mode), item in sorted(results.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        row = dict(range_nm=range_nm, mode=mode, seeds=len(item["seed"]))
        for key in ("n_final", "runtime_s", "attach_events"):
            mean, sem = _mean_sem(item[key])
            row[f"{key}_mean"], row[f"{key}_sem"] = float(mean), float(sem)
        for nm in near_nm:
            d, mean, sem = attach_fraction_within(item, max_nm=nm)
            row[f"attach_within_{nm:g}nm_mean"], row[f"attach_within_{nm:g}nm_sem"] = float(mean[-1]), float(sem[-1])
            d, mean, sem = attach_enrichment(item, max_nm=nm)
            row[f"attach_enrichment_{nm:g}nm_mean"], row[f"attach_enrichment_{nm:g}nm_sem"] = float(mean[-1]), float(sem[-1])
        rows.append(row)
    return rows


# ------------------------------------------------------------------ plotting (reads loaded results only)
# Palette from the dataviz reference: range is an ordered magnitude, so it takes the blue sequential ramp
# (steps 250 to 650, light = short range, dark = long range); the mode is a secondary encoding (solid =
# graded, dashed = hard). The off baseline is neutral gray ink. Validated with validate_palette.js --ordinal.

SURFACE, INK, INK_2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
RANGE_RAMP = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
BASELINE_COLOR = INK_2
MODE_LINESTYLE = {"graded": "-", "hard": (0, (5, 2))}
MODE_TITLE = {"graded": "Graded penalty", "hard": "Hard exclusion"}


def _hex_to_rgb(color):
    return np.array([int(color[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)


def range_color(range_nm, ranges):
    """Color of a range on the ordinal ramp; interpolates when there are more than five ranges."""
    ranges = sorted(ranges)
    t = 0.5 if len(ranges) == 1 else ranges.index(range_nm) / (len(ranges) - 1)
    ramp = np.array([_hex_to_rgb(c) for c in RANGE_RAMP])
    channels = [np.interp(t, np.linspace(0, 1, len(ramp)), ramp[:, k]) for k in range(3)]
    return "#%02x%02x%02x" % tuple(int(round(c)) for c in channels)


def _positive_ranges(results):
    return sorted({r for (r, m) in results if r > 0})


def _new_figure(ncols=2, width=9.0, height=3.9, sharey=True):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(width, height), facecolor=SURFACE, layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, ncols, sharey=sharey))
    return fig, axes


def _style(ax, xlabel, ylabel=None, title=None):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.set_xlabel(xlabel, color=INK_2, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK_2, fontsize=10)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left")


def _legend(fig, results, baseline_label="Off (0 nm)"):
    from matplotlib.lines import Line2D
    ranges = _positive_ranges(results)
    handles = []
    if (0.0, "off") in results:
        handles.append(Line2D([], [], color=BASELINE_COLOR, linewidth=1.8, label=baseline_label))
    handles += [Line2D([], [], color=range_color(r, ranges), linewidth=2.0, label=f"{r:g} nm") for r in ranges]
    fig.legend(handles=handles, loc="outside lower center", ncols=len(handles), frameon=False,
               fontsize=9, labelcolor=INK_2)


def _band(ax, x, mean, sem, color, linestyle="-", linewidth=1.8, alpha=0.16, zorder=2):
    ax.plot(x, mean, color=color, linestyle=linestyle, linewidth=linewidth, zorder=zorder + 1)
    ax.fill_between(x, mean - sem, mean + sem, color=color, alpha=alpha, linewidth=0, zorder=zorder)


def plot_binding_profiles(results, max_nm=40.0):
    """Relative binding rate of a free site at distance d from one occupied site (from the model)."""
    fig, axes = _new_figure()
    ranges = _positive_ranges(results)
    for ax, mode in zip(axes, ("graded", "hard")):
        ax.axhline(1.0, color=BASELINE_COLOR, linewidth=1.2)
        for r in ranges:
            item = results.get((r, mode))
            if item is None:
                continue
            d, w = item["profile_distance_nm"], item["profile_weight"]
            keep = (d > 0) & (d <= max_nm)
            ax.plot(d[keep], w[keep], color=range_color(r, ranges), linestyle=MODE_LINESTYLE[mode], linewidth=2.0)
        _style(ax, "Distance from an occupied site (nm)", "Relative binding rate of a free site" if mode == "graded" else None,
               MODE_TITLE[mode])
        ax.set_xlim(0, max_nm)
        ax.set_ylim(-0.03, 1.08)
    _legend(fig, results, baseline_label="Off (rate 1)")
    return fig


def plot_pair_profiles(results, bin_nm=4.0, max_nm=64.0):
    """Occupancy near an occupied head: same-lane pair density relative to uniform placement, mean and SEM."""
    fig, axes = _new_figure()
    ranges = _positive_ranges(results)
    base = results.get((0.0, "off"))
    for ax, mode in zip(axes, ("graded", "hard")):
        if base is not None:
            x, mean, sem = pair_correlation(base, bin_nm, max_nm)
            _band(ax, x, mean, sem, BASELINE_COLOR, alpha=0.14, zorder=1)
        for r in ranges:
            item = results.get((r, mode))
            if item is not None:
                x, mean, sem = pair_correlation(item, bin_nm, max_nm)
                _band(ax, x, mean, sem, range_color(r, ranges), MODE_LINESTYLE[mode])
        _style(ax, "Distance between two heads on the same lane (nm)",
               "Pair density vs uniform in the overlap, g(d)" if mode == "graded" else None, MODE_TITLE[mode])
        ax.set_xlim(0, max_nm)
        ax.set_ylim(bottom=0)
    _legend(fig, results)
    return fig


def plot_attach_profiles(results, max_nm=40.0):
    """
    Empirical binding profile: observed over expected number of first-head attachments within d of a bound
    head, where expected assumes uniform binding at the same occupancy. 1 = no steric effect, whatever the
    density, so the steric-off line should sit at 1.
    """
    fig, axes = _new_figure()
    ranges = _positive_ranges(results)
    base = results.get((0.0, "off"))
    for ax, mode in zip(axes, ("graded", "hard")):
        ax.axhline(1.0, color=GRID, linewidth=1.4, zorder=0)
        if base is not None:
            x, mean, sem = attach_enrichment(base, max_nm)
            _band(ax, x, mean, sem, BASELINE_COLOR, alpha=0.14, zorder=1)
        for r in ranges:
            item = results.get((r, mode))
            if item is not None:
                x, mean, sem = attach_enrichment(item, max_nm)
                _band(ax, x, mean, sem, range_color(r, ranges), MODE_LINESTYLE[mode])
        _style(ax, "Distance to the nearest bound head (nm)",
               "New first heads within d, observed / expected" if mode == "graded" else None, MODE_TITLE[mode])
        ax.set_xlim(0, max_nm)
        ax.set_ylim(bottom=0)
    _legend(fig, results)
    return fig


def plot_summary(results):
    """Final PRC1 count and simulation time against range (mean and SEM across seeds)."""
    fig, axes = _new_figure(width=9.0, height=3.5, sharey=False)   # counts and seconds need their own scales
    rows = summarize(results)
    ranges = _positive_ranges(results)
    for ax, key, label in ((axes[0], "n_final", "PRC1 bound at the end of the run"),
                           (axes[1], "runtime_s", "Simulation time per run (s)")):
        base = next((r for r in rows if r["mode"] == "off"), None)
        if base is not None:
            ax.axhline(base[f"{key}_mean"], color=BASELINE_COLOR, linewidth=1.2, label="Off (0 nm)")
            ax.axhspan(base[f"{key}_mean"] - base[f"{key}_sem"], base[f"{key}_mean"] + base[f"{key}_sem"],
                       color=BASELINE_COLOR, alpha=0.12, linewidth=0)
        for mode, marker, filled in (("graded", "o", True), ("hard", "s", False)):
            pts = [(r["range_nm"], r[f"{key}_mean"], r[f"{key}_sem"]) for r in rows if r["mode"] == mode]
            if pts:
                x, y, e = zip(*pts)
                ax.errorbar(x, y, yerr=e, color=INK, marker=marker, markersize=7, linewidth=1.4, capsize=3,
                            markerfacecolor=INK if filled else SURFACE, markeredgewidth=1.6,
                            label=MODE_TITLE[mode])
        ax.set_xscale("log", base=2)
        ax.set_xticks(ranges)
        ax.set_xticklabels([f"{r:g}" for r in ranges])
        ax.minorticks_off()
        _style(ax, "Steric range (nm)", label)
        ax.set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK_2, loc="lower left")
    return fig


def make_figures(results_dir=DEFAULT_OUT, out_dir=None):
    """Regenerate every figure from the saved results (no simulation). Returns the written paths."""
    results_dir = Path(results_dir)
    out_dir = Path(out_dir) if out_dir else results_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = load_results(results_dir)
    paths = []
    for name, fig in (("binding_profile", plot_binding_profiles(results)),
                      ("pair_profile", plot_pair_profiles(results)),
                      ("attach_profile", plot_attach_profiles(results)),
                      ("summary", plot_summary(results))):
        path = out_dir / f"{name}.png"
        fig.savefig(path, dpi=150)
        paths.append(path)
    return paths


# ------------------------------------------------------------------ command line

def print_summary(results):
    print(f"{'range_nm':>8} {'mode':>7} {'seeds':>5} {'n_final':>14} {'runtime_s':>12} "
          f"{'attach<=8nm':>13} {'obs/exp 8nm':>13} {'obs/exp 16nm':>13}")
    for r in summarize(results):
        print(f"{r['range_nm']:8.1f} {r['mode']:>7} {r['seeds']:5d} "
              f"{r['n_final_mean']:7.1f}+-{r['n_final_sem']:4.1f} {r['runtime_s_mean']:6.2f}+-{r['runtime_s_sem']:4.2f} "
              f"{r['attach_within_8nm_mean']:7.3f}+-{r['attach_within_8nm_sem']:5.3f} "
              f"{r['attach_enrichment_8nm_mean']:7.2f}+-{r['attach_enrichment_8nm_sem']:4.2f} "
              f"{r['attach_enrichment_16nm_mean']:7.2f}+-{r['attach_enrichment_16nm_sem']:4.2f}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1], formatter_class=argparse.RawTextHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run the sweep and save results")
    run.add_argument("--out", default=str(DEFAULT_OUT))
    run.add_argument("--ranges", type=float, nargs="+", help="steric ranges in nm (include 0 for the baseline)")
    run.add_argument("--modes", nargs="+", choices=["graded", "hard"])
    run.add_argument("--seeds", type=int, help="number of seeds (0..n-1)")
    run.add_argument("--end-time", type=float)
    run.add_argument("--max-steps", type=int)
    run.add_argument("--binding-rate", type=float, help="initial_binding_rate_per_site")
    run.add_argument("--length", type=float, help="microtubule length in nm")
    run.add_argument("--offset", type=float, help="microtubule offset in nm")
    plot = sub.add_parser("plot", help="regenerate figures and print the summary from saved results")
    plot.add_argument("--results", default=str(DEFAULT_OUT))
    plot.add_argument("--figures", default=None)
    args = parser.parse_args(argv)

    if args.command == "run":
        cfg = SweepConfig()
        if args.ranges is not None: cfg.ranges_nm = tuple(args.ranges)
        if args.modes is not None: cfg.modes = tuple(args.modes)
        if args.seeds is not None: cfg.seeds = tuple(range(args.seeds))
        if args.end_time is not None: cfg.end_time = args.end_time
        if args.max_steps is not None: cfg.max_steps = args.max_steps
        if args.binding_rate is not None: cfg.rates = dict(cfg.rates, initial_binding_rate_per_site=args.binding_rate)
        if args.length is not None: cfg.lattice = dict(cfg.lattice, microtubule_length=args.length)
        if args.offset is not None: cfg.lattice = dict(cfg.lattice, microtubule_offset=args.offset)
        cfg.validate()
        print(f"running {len(cfg.settings())} settings x {len(cfg.seeds)} seeds")
        run_sweep(cfg, args.out)
    else:
        results = load_results(args.results)
        print_summary(results)
        for path in make_figures(args.results, args.figures):
            print("wrote", path)


if __name__ == "__main__":
    main()
