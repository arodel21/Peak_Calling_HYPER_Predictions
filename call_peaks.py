#!/usr/bin/env python3
## Call peaks from a bigWig or bedGraph against a log-normal background model.

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import brentq, fsolve
from scipy.stats import lognorm

#: The three log-normal model options.
LOGNORMAL_MODELS = ("hybrid", "bulk", "full")

BEDGRAPH_COLUMNS = ["chrom", "start", "end", "value"]
BED6_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


class Stages:
    """Elapsed time per stage, so a slow run says which part was slow."""

    def __init__(self) -> None:
        self.laps: list[tuple[str, float]] = []
        self._last = self._start = time.perf_counter()

    def lap(self, name: str) -> None:
        now = time.perf_counter()
        self.laps.append((name, now - self._last))
        self._last = now

    def report(self) -> str:
        total = time.perf_counter() - self._start
        parts = "  ".join(f"{n} {s:.1f}s" for n, s in self.laps if s >= 0.05)
        return f"total {total:.1f}s  ({parts})"


# --------------------------------------------------------------------------- #
# Arguments
# --------------------------------------------------------------------------- #

def parse_loc(text: str):
    """--loc takes a number, or one of the keywords."""
    lowered = str(text).strip().lower()
    if lowered in {"auto", "free", "zero"}:
        return lowered
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"expected a number or auto/free/zero, got {text!r}"
        ) from None
    if value < 0:
        raise argparse.ArgumentTypeError(
            f"a log-normal cannot start below zero, got {value}"
        )
    return value


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("signal", type=Path,
                        help="input bigWig (.bw) or bedGraph (4 columns); the "
                             "format is detected from the file itself")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="output BED of called peaks")

    reading = parser.add_argument_group("reading the signal")
    reading.add_argument(
        "--chroms", nargs="*", default=None, metavar="CHROM",
        help="only these chromosomes (default: all)",
    )
    reading.add_argument(
        "--bin-size", type=int, default=None, metavar="BP",
        help="re-bin a bigWig into fixed windows of this width. Leave unset to "
             "use the intervals stored in the file, which for a binned track "
             "are already the bins.",
    )

    fitting = parser.add_argument_group("background and fitting")
    fitting.add_argument(
        "--background", default="auto", metavar="VALUE",
        help="drop bins at or below this signal value before fitting. 'auto' "
             "fits the model above a dozen candidate cuts and keeps the "
             "smallest one that fits the kept data best (see --plot to check "
             "it). (default: %(default)s)",
    )
    fitting.add_argument(
        "--model", choices=list(LOGNORMAL_MODELS), default="hybrid",
        help="which background model a bin is judged against. 'bulk' is the "
             "log-normal fitted to the bottom --bulk-quantile of bins; 'full' "
             "is the one fitted to every kept bin, whose tail the peaks "
             "themselves make heavier; 'hybrid' is bulk below their crossing "
             "point and full above it (default: %(default)s)",
    )
    fitting.add_argument(
        "--loc", type=parse_loc, default="auto", metavar="VALUE",
        help="where the log-normal starts: a number to pin it there, or "
             "'zero', 'free' (let it float, hybrid model only) or 'auto' (try "
             "both and keep the better fit, hybrid model only). A number must "
             "be at least 0 and below the smallest kept value. Under --model "
             "hybrid this sets the bulk fit, and the tail fit stays pinned at "
             "zero (default: %(default)s)",
    )
    fitting.add_argument(
        "--bulk-quantile", default="auto", metavar="Q",
        help="the bulk fit uses values up to this quantile. 'auto' picks the "
             "one whose fit matches the core best; too high and the signal "
             "tail widens the fit and flattens its peak below the data "
             "(default: %(default)s)",
    )

    calling = parser.add_argument_group("peak calling")
    calling.add_argument(
        "--threshold-mode", choices=["density", "pvalue"], default="density",
        help="'density' calls a bin when the hybrid PDF at its value falls "
             "below --density-threshold; 'pvalue' uses the tail survival "
             "function, which is a real p-value (default: %(default)s)",
    )
    calling.add_argument("--density-threshold", type=float, default=1e-3,
                         metavar="D", help="for --threshold-mode density "
                                           "(default: %(default)s)")
    calling.add_argument("--pvalue", type=float, default=1e-4, metavar="P",
                         help="for --threshold-mode pvalue (default: %(default)s)")
    calling.add_argument("--fold-change", type=float, default=5.0, metavar="FC",
                         help="minimum value / background mean (default: %(default)s)")

    other = parser.add_argument_group("other")
    other.add_argument("--merge-adjacent", action="store_true",
                       help="merge touching called bins into single peaks")
    other.add_argument("--name-prefix", default="peak",
                       help="peak name prefix (default: %(default)s)")
    other.add_argument("--plot", type=Path, default=None,
                       help="save a diagnostic plot of the fits and thresholds")
    other.add_argument("--plot-scale", choices=["linear", "log"], default="linear",
                       help="axes for the fit panel. Linear shows the log-normal "
                            "shape of binned coverage; log suits a track spanning "
                            "orders of magnitude (default: %(default)s)")
    other.add_argument("--bed6", action="store_true",
                       help="write 6 BED columns instead of 4")
    return parser.parse_args(argv)


# --------------------------------------------------------------------------- #
# Reading the signal
# --------------------------------------------------------------------------- #

BIGWIG_MAGIC = (b"\x26\xfc\x8f\x88", b"\x88\x8f\xfc\x26")

def is_bigwig(path: Path) -> bool:
    """Detect by magic number rather than extension, which can be anything."""
    with open(path, "rb") as handle:
        return handle.read(4) in BIGWIG_MAGIC


def read_bedgraph(path: Path, chroms: list[str] | None) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", header=None, comment="#",
                        names=BEDGRAPH_COLUMNS, usecols=[0, 1, 2, 3])
    if chroms:
        table = table[table["chrom"].isin(chroms)].reset_index(drop=True)
    return table


def read_bigwig(path: Path, chroms: list[str] | None,
                bin_size: int | None) -> pd.DataFrame:
    #Read a bigWig as the interval table it is.
    import pyBigWig

    handle = pyBigWig.open(str(path))
    try:
        available = handle.chroms()
        wanted = chroms if chroms else list(available)
        missing = [c for c in wanted if c not in available]
        if missing:
            raise KeyError(f"not in {path}: {missing}")

        frames = []
        for chrom in wanted:
            length = int(available[chrom])
            if bin_size:
                starts, ends, values = _binned(handle, chrom, length, bin_size)
            else:
                starts, ends, values = _native(handle, chrom)
            if starts.size == 0:
                continue
            keep = np.isfinite(values)
            frames.append(pd.DataFrame({
                "chrom": chrom, "start": starts[keep],
                "end": ends[keep], "value": values[keep],
            }))
    finally:
        handle.close()

    if not frames:
        raise ValueError(f"no usable intervals in {path}")
    return pd.concat(frames, ignore_index=True)


def _native(handle, chrom: str):
    intervals = handle.intervals(chrom)
    if not intervals:
        return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, float)
    block = np.asarray(intervals, dtype=np.float64)
    return (block[:, 0].astype(np.int64), block[:, 1].astype(np.int64),
            block[:, 2].copy())


def _binned(handle, chrom: str, length: int, bin_size: int, chunk: int = 200_000):
    """Fixed-width bins, a chunk at a time so a whole chromosome is never held."""
    starts_all, values_all = [], []
    for begin in range(0, length, chunk * bin_size):
        stop = min(begin + chunk * bin_size, length)
        n_bins = (stop - begin) // bin_size
        if n_bins == 0:
            continue
        end = begin + n_bins * bin_size
        stats = handle.stats(chrom, begin, end, nBins=n_bins, type="mean")
        values_all.append(np.array([np.nan if v is None else v for v in stats],
                                   dtype=float))
        starts_all.append(begin + np.arange(n_bins, dtype=np.int64) * bin_size)
    if not starts_all:
        return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, float)
    starts = np.concatenate(starts_all)
    return starts, starts + bin_size, np.concatenate(values_all)


def read_signal(path: Path, chroms: list[str] | None,
                bin_size: int | None) -> tuple[pd.DataFrame, str]:
    """Read either format into one chrom/start/end/value table."""
    if is_bigwig(path):
        table = read_bigwig(path, chroms, bin_size)
        kind = f"bigWig{f' re-binned at {bin_size} bp' if bin_size else ''}"
    else:
        if bin_size:
            raise ValueError("--bin-size only applies to a bigWig input")
        table = read_bedgraph(path, chroms)
        kind = "bedGraph"

    widths = np.unique(table["end"].to_numpy() - table["start"].to_numpy())
    if widths.size > 1:
        print(
            f"WARNING: intervals are not all the same width ({widths.size} "
            f"distinct widths, {widths.min()}-{widths.max()} bp). Each interval "
            "counts once however long it is, so the fitted distribution is "
            "skewed. Pass --bin-size to re-bin.",
            file=sys.stderr,
        )
    return table, kind


# --------------------------------------------------------------------------- #
# Background level
# --------------------------------------------------------------------------- #


#: Fractions of the track the automatic scan is willing to discard as background.
BACKGROUND_CANDIDATES = (0.0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05,
                         0.075, 0.10, 0.15, 0.20, 0.25)

#: How much worse than the best a smaller cut may score and still be preferred.
#: Throwing bins away is not free, so a cut has to earn its keep.
BACKGROUND_TOLERANCE = 1.05


def estimate_background(values: np.ndarray, loc: float = 0.0,
                        candidates=BACKGROUND_CANDIDATES) -> tuple[float, str]:
    """Choose the cut by fitting the real model above each candidate.

    A predicted track often holds two populations: a low one the model put
    everywhere it saw nothing, and the real signal. Where to cut between them
    used to be guessed from the shape of a histogram, which says nothing about
    whether the resulting fit is any good -- on a track whose two populations
    overlap it finds one mode and cuts nothing.

    So the cut is scored directly instead: fit above each candidate and measure
    how far the fitted distribution is from the kept data (the
    Kolmogorov-Smirnov distance). Cutting too low leaves the low population in,
    which widens the fit; cutting too high removes the left half of the
    background and makes ordinary bins look extreme. Both show up as a worse
    KS, so the scan has a minimum rather than improving forever, and the
    smallest cut within ``BACKGROUND_TOLERANCE`` of the best one wins.

    What is fitted at each candidate is a single log-normal over the kept bins,
    which is exactly the assumption the background cut exists to make true.
    The KS it reports is worth reading: a large one at every cut says no cut
    makes the kept bins log-normal, whatever value you pick.

    A 50k sample per candidate is plenty for a KS distance and keeps the whole
    scan to a second or two.
    """
    from scipy.stats import kstest

    positive = values[values > 0]
    if positive.size == 0:
        return 0.0, "no positive values; nothing removed"

    rng = np.random.default_rng(0)
    scored = []
    for share in candidates:
        cut = 0.0 if share == 0 else float(np.quantile(positive, share))
        kept = positive[positive > cut]
        if kept.size < 10_000:
            continue
        check = (kept if kept.size <= 50_000
                 else rng.choice(kept, 50_000, replace=False))
        try:
            shape, start, scale = lognorm.fit(
                check, floc=loc if isinstance(loc, float) else 0.0)
        except (ValueError, FloatingPointError):
            continue
        survival = lambda x: lognorm.sf(x, s=shape, loc=start, scale=scale)
        # Conditional on being above the cut, which is the only thing the kept
        # data can speak to.
        above = float(survival(max(cut, kept.min() * 0.999)))
        if not np.isfinite(above) or above <= 0:
            continue
        stat = kstest(check, lambda x: 1.0 - survival(x) / above).statistic
        scored.append((cut, float(stat), share))

    if not scored:
        return 0.0, "no candidate cut could be fitted; nothing removed"

    best = min(s for _, s, _ in scored)
    cut, stat, share = next((c, s, f) for c, s, f in scored
                            if s <= best * BACKGROUND_TOLERANCE)
    if cut == 0.0:
        return 0.0, (f"scan of {len(scored)} cuts: cutting nothing fits best "
                     f"(KS {stat:.4f}), so the track has one population")
    return cut, (f"scan of {len(scored)} cuts: best KS {stat:.4f} at the "
                 f"{share:.1%} quantile; removes {np.mean(values <= cut):.1%} "
                 f"of bins")


# --------------------------------------------------------------------------- #
# Fitting
# --------------------------------------------------------------------------- #


def fit_lognormals(data: np.ndarray, bulk_quantile: float,
                   mode: str = "auto") -> dict[str, float]:
    """Fit the bulk and the full distribution"""
    bulk = data[data <= np.quantile(data, bulk_quantile)]
    if bulk.size < 10:
        raise ValueError(
            f"only {bulk.size} values below the {bulk_quantile} quantile; "
            "the background cut is probably too high"
        )

    from scipy.stats import kstest

    def fit_fixed(sample):
        return lognorm.fit(sample, floc=0)

    def fit_free(sample):
        shape, loc, scale = lognorm.fit(sample)
        # Coverage cannot be negative, and the density must be positive
        if not np.isfinite([shape, loc, scale]).all() or not 0 <= loc < sample.min():
            return None
        return shape, loc, scale

    def fit(sample: np.ndarray) -> tuple[float, float, float]:
        if isinstance(mode, float):
            return tuple(float(x) for x in lognorm.fit(sample, floc=mode))
        if mode == "zero":
            return tuple(float(x) for x in fit_fixed(sample))

        free = fit_free(sample)
        if mode == "free":
            return tuple(float(x) for x in (free or fit_fixed(sample)))

        # auto: take whichever describes the sample better.
        fixed = fit_fixed(sample)
        if free is None:
            return tuple(float(x) for x in fixed)
        sub = sample if sample.size <= 50_000 else np.random.default_rng(0).choice(
            sample, 50_000, replace=False)
        scores = [
            (kstest(sub, "lognorm", args=(c[0],), kwds=None).statistic
             if False else kstest(sub, lambda x, c=c: lognorm.cdf(
                 x, s=c[0], loc=c[1], scale=c[2])).statistic, c)
            for c in (fixed, free)
        ]
        return tuple(float(x) for x in min(scores, key=lambda s: s[0])[1])

    shape_bulk, loc_bulk, scale_bulk = fit(bulk)
    shape_full, loc_full, scale_full = (float(x) for x in fit_fixed(data))
    return {
        "shape_bulk": shape_bulk, "loc_bulk": loc_bulk, "scale_bulk": scale_bulk,
        "shape_full": shape_full, "loc_full": loc_full, "scale_full": scale_full,
    }


def bulk_args(fit: dict[str, float]) -> dict[str, float]:
    return {"s": fit["shape_bulk"], "loc": fit["loc_bulk"], "scale": fit["scale_bulk"]}


def full_args(fit: dict[str, float]) -> dict[str, float]:
    return {"s": fit["shape_full"], "loc": fit["loc_full"], "scale": fit["scale_full"]}


def choose_bulk_quantile(data: np.ndarray, candidates=None,
                         mode="auto") -> tuple[float, str]:
    """Pick the quantile whose bulk fit matches the core best"""
    from scipy.stats import kstest

    if candidates is None:
        candidates = np.round(np.arange(0.70, 0.995, 0.01), 3)

    rng = np.random.default_rng(0)
    sample = data if data.size <= 100_000 else rng.choice(data, 100_000, replace=False)
    ordered = np.sort(data)
    sample = np.sort(sample)

    best_q, best_stat = None, np.inf
    for q in candidates:
        cut = float(ordered[min(int(q * (ordered.size - 1)), ordered.size - 1)])
        core = sample[:np.searchsorted(sample, cut, side="right")]
        if core.size < 100:
            continue
        # Score against the same model the run will use, or the quantile is
        # chosen for a fit that is never made.
        try:
            if isinstance(mode, float):
                shape, loc, scale = lognorm.fit(core, floc=mode)
            elif mode == "zero":
                shape, loc, scale = lognorm.fit(core, floc=0)
            else:
                shape, loc, scale = lognorm.fit(core)
                if not np.isfinite([shape, loc, scale]).all() or loc >= core.min():
                    shape, loc, scale = lognorm.fit(core, floc=0)
        except Exception:
            continue
        total = lognorm.cdf(cut, s=shape, loc=loc, scale=scale)
        if not np.isfinite(total) or total <= 0:
            continue
        stat = kstest(
            core, lambda x: lognorm.cdf(x, s=shape, loc=loc, scale=scale) / total
        ).statistic
        if stat < best_stat:
            best_q, best_stat = float(q), float(stat)

    if best_q is None:
        return 0.98, "could not be chosen; fell back to 0.98"
    return best_q, f"chosen by KS distance ({best_stat:.4f})"


def density_difference(x, fit: dict[str, float]):
    """Zero where the two fitted densities cross."""
    return lognorm.pdf(x, **bulk_args(fit)) - lognorm.pdf(x, **full_args(fit))


def find_split(data: np.ndarray, fit: dict[str, float]) -> float:
    """The signal value where the bulk and full densities cross"""
    low, high = np.quantile(data, 0.5), np.quantile(data, 0.999)
    try:
        return float(brentq(density_difference, low, high, args=(fit,)))
    except ValueError:
        guess = np.quantile(data, 0.85)
        return float(fsolve(density_difference, guess, args=(fit,), xtol=1e-10)[0])


def hybrid_density(x, split: float, fit: dict[str, float]) -> np.ndarray:
    """Bulk density below the split, full density above it."""
    x = np.asarray(x, dtype=float)
    return np.where(x < split,
                    lognorm.pdf(x, **bulk_args(fit)),
                    lognorm.pdf(x, **full_args(fit)))


def tail_pvalue(x, fit: dict[str, float]) -> np.ndarray:
    """Survival function of the full fit: the probability of a value this high."""
    return lognorm.sf(x, **full_args(fit))


def model_density(x, model: str, split: float, fit: dict[str, float]) -> np.ndarray:
    """The density the chosen model puts at each value."""
    if model == "bulk":
        return lognorm.pdf(np.asarray(x, dtype=float), **bulk_args(fit))
    if model == "full":
        return lognorm.pdf(np.asarray(x, dtype=float), **full_args(fit))
    return hybrid_density(x, split, fit)


def model_pvalue(x, model: str, fit: dict[str, float]) -> np.ndarray:
    """The upper-tail probability the chosen model puts at each value"""
    args = bulk_args(fit) if model == "bulk" else full_args(fit)
    return lognorm.sf(np.asarray(x, dtype=float), **args)


# --------------------------------------------------------------------------- #
# Calling
# --------------------------------------------------------------------------- #


def call_peaks(data: np.ndarray, split: float, fit: dict[str, float],
               args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Boolean mask of called bins, plus the score for each bin."""
    background_mean = float(np.mean(data[data <= split]))
    fold_change = data / background_mean

    if args.threshold_mode == "pvalue":
        statistic = model_pvalue(data, args.model, fit)
        extreme = statistic < args.pvalue
    else:
        # Interpolating the density keeps this identical to evaluating it
        # directly, but is far cheaper over millions of bins.
        grid = np.linspace(max(data.min(), 1e-9), data.max() * 1.1, 4000)
        spline = interp1d(grid, model_density(grid, args.model, split, fit),
                          kind="cubic", bounds_error=False, fill_value=0.0)
        statistic = spline(data)
        extreme = statistic < args.density_threshold

    mask = extreme & (fold_change > args.fold_change)
    with np.errstate(divide="ignore"):
        score = -np.log10(np.clip(statistic, 1e-300, None))
    return mask, score


def merge_adjacent(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse runs of touching bins into one peak, keeping the best score."""
    frame = frame.sort_values(["chrom", "start"]).reset_index(drop=True)
    new_group = (
        (frame["chrom"] != frame["chrom"].shift())
        | (frame["start"] > frame["end"].shift())
    )
    frame["_group"] = new_group.cumsum()
    merged = frame.groupby("_group").agg(
        chrom=("chrom", "first"), start=("start", "min"),
        end=("end", "max"), score=("score", "max"),
    )
    return merged.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #


def save_plot(path: Path, raw: np.ndarray, data: np.ndarray, background: float,
              split: float | None, fit: dict[str, float] | None,
              scale: str = "linear", model: str = "hybrid") -> None:
    """The fitted distribution over the kept bins."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Linear by default: binned coverage is log-normal and that shape reads
    # directly on linear axes.
    if scale == "log":
        grid = np.logspace(np.log10(data.min()), np.log10(data.max()), 1000)
        bins = np.logspace(np.log10(data.min()), np.log10(data.max()), 200)
    else:
        upper = float(np.quantile(data, 0.99))
        if split is not None:
            upper = max(upper, split * 2.5)
        grid = np.linspace(max(data.min() * 0.5, 1e-9), upper, 1000)
        bins = np.linspace(0, upper, 120)

    ax.hist(data, bins=bins, density=True, color="#8fb8d8",
            edgecolor="none", label="data")
    # All three are drawn whichever one is in use
    if model != "bulk":
        ax.plot(grid, lognorm.pdf(grid, **bulk_args(fit)),
                "g--", lw=1.6, label="bulk fit (peak)")
    if model != "full":
        ax.plot(grid, lognorm.pdf(grid, **full_args(fit)),
                "r--", lw=1.6, label="full fit (tail)")
    if model != "hybrid":
        ax.plot(grid, hybrid_density(grid, split, fit), color="grey",
                lw=1.2, ls="-.", label="hybrid fit")
    ax.plot(grid, model_density(grid, model, split, fit), color="purple",
            lw=3, label=f"{model.upper()} fit (used for calling)")
    ax.axvline(split, color="navy", lw=1.5, ls=":", label=f"split {split:.4g}")
    ax.set(xlabel="binned coverage", ylabel="density")
    if scale == "log":
        ax.set(xscale="log", yscale="log")
        floor = float(lognorm.pdf(data.max(), **full_args(fit)))
        ax.set_ylim(bottom=max(1e-8, floor / 10))
    else:
        ax.set_xlim(0, grid[-1])
    ax.legend(frameon=False)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    stages = Stages()

    if not args.signal.is_file():
        print(f"ERROR: no such file: {args.signal}", file=sys.stderr)
        return 1

    try:
        table, kind = read_signal(args.signal, args.chroms, args.bin_size)
    except (KeyError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    raw = table["value"].to_numpy(dtype=float)
    stages.lap("read")
    print(f"read {len(table):,} bins from {args.signal} ({kind})")

    if args.background == "auto":
        background, how = estimate_background(
            raw, loc=args.loc if isinstance(args.loc, float) else 0.0)
        stages.lap("background scan")
        print(f"background cut: {background:.6g}  ({how})")
    else:
        background = float(args.background)
        print(f"background cut: {background:.6g}  (given)")

    kept = table[table["value"] > background].reset_index(drop=True)
    data = kept["value"].to_numpy(dtype=float)
    if data.size < 100:
        print(f"ERROR: only {data.size} bins above the background cut", file=sys.stderr)
        return 1
    print(f"kept {len(kept):,} bins ({len(kept)/len(table):.1%})")

    if isinstance(args.loc, float) and args.loc >= data.min():
        print(
            f"ERROR: --loc {args.loc:g} is not below the smallest kept value "
            f"({data.min():.6g}); the density would be zero where there is data",
            file=sys.stderr,
        )
        return 1

    if args.bulk_quantile == "auto":
        bulk_quantile, why = choose_bulk_quantile(data, mode=args.loc)
        print(f"bulk quantile: {bulk_quantile:.3f}  ({why})")
    else:
        bulk_quantile = float(args.bulk_quantile)
        print(f"bulk quantile: {bulk_quantile:.3f}  (given)")

    try:
        fit = fit_lognormals(data, bulk_quantile, mode=args.loc)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    split = find_split(data, fit)
    print(f"model: {args.model}"
          + ("" if args.model == "hybrid" else
             f" (the {args.model} fit alone; the split and the fold change "
             "are unchanged)"))
    print(f"bulk  fit: shape {fit['shape_bulk']:.4f}  loc {fit['loc_bulk']:.4f}  "
          f"scale {fit['scale_bulk']:.4f}")
    print(f"full  fit: shape {fit['shape_full']:.4f}  loc {fit['loc_full']:.4f}  "
          f"scale {fit['scale_full']:.4f}")
    print(f"split at {split:.4g} "
          f"(quantile {np.mean(data <= split):.4f} of the kept bins)")
    stages.lap("fit")
    mask, score = call_peaks(data, split, fit, args)
    stages.lap("call")

    print(f"{int(mask.sum()):,} bins called "
          f"({args.threshold_mode} threshold, fold change > {args.fold_change})")

    peaks = kept.loc[mask, ["chrom", "start", "end"]].copy()
    peaks["score"] = score[mask]
    if args.merge_adjacent:
        before = len(peaks)
        peaks = merge_adjacent(peaks)
        print(f"merged {before:,} bins into {len(peaks):,} peaks")

    peaks = peaks.sort_values("score", ascending=False).reset_index(drop=True)
    peaks["name"] = [f"{args.name_prefix}_{i + 1}" for i in range(len(peaks))]
    peaks["strand"] = "."

    columns = BED6_COLUMNS if args.bed6 else ["chrom", "start", "end", "name"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    peaks[columns].to_csv(args.output, sep="\t", header=False, index=False)
    print(f"wrote {len(peaks):,} peaks -> {args.output}")

    if args.plot:
        save_plot(args.plot, raw, data, background, split, fit,
                  args.plot_scale, args.model)
        stages.lap("plot")
        print(f"plot -> {args.plot}")
    print(stages.report())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
