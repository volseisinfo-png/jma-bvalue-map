#!/usr/bin/env python3
"""Make configurable gridded b-value maps from JMA hypocenter records."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sys
import zipfile
import calendar
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator

import numpy as np


MAG_MIN_TENTH = -30
MAG_MAX_TENTH = 99
MAG_CENTERS = np.arange(MAG_MIN_TENTH, MAG_MAX_TENTH + 1, dtype=np.int16) / 10.0


def decode_jma_magnitude(raw: bytes) -> float | None:
    """Decode columns 53-54 (Magnitude 1), including JMA negative notation."""
    text = raw.decode("ascii", errors="strict")
    if text == "  ":
        return None
    if len(text) != 2:
        return None
    if text[0] == "-" and text[1].isdigit():       # -1 ... -9 => -0.1 ... -0.9
        return -int(text[1]) / 10.0
    if "A" <= text[0] <= "I" and text[1].isdigit():
        # A0=-1.0, A9=-1.9, B0=-2.0, ... (future-proof through I9).
        return -((ord(text[0]) - ord("A") + 1) + int(text[1]) / 10.0)
    if text.isdigit():
        return int(text) / 10.0
    return None


def _ascii_number(raw: bytes) -> float | None:
    try:
        text = raw.decode("ascii").strip()
        return float(text) if text else None
    except (UnicodeDecodeError, ValueError):
        return None


def _implied_decimal(raw: bytes, decimals: int) -> float | None:
    """Parse a Fortran fixed-width field such as F4.2 (`3025` => 30.25).

    JMA may leave the fractional positions blank for a fixed hypocenter.
    """
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError:
        return None
    if not text.strip():
        return None
    if "." in text:
        try:
            return float(text)
        except ValueError:
            return None
    fractional = text[-decimals:]
    if fractional.isspace():
        try:
            return float(text[:-decimals].strip())
        except ValueError:
            return None
    compact = text.strip()
    try:
        return int(compact) / (10 ** decimals)
    except ValueError:
        return None


def parse_hypocenter_record(line: bytes) -> tuple[float, float, float] | None:
    """Return (latitude, longitude, magnitude1) for a valid J record.

    Python slices are zero-based; these correspond exactly to JMA columns
    01, 22-24, 25-28, 33-36, 37-40, and 53-54.
    """
    line = line.rstrip(b"\r\n")
    if len(line) < 96 or line[0:1] != b"J":
        return None
    lat_deg = _ascii_number(line[21:24])
    lat_min = _implied_decimal(line[24:28], 2)
    lon_deg = _ascii_number(line[32:36])
    lon_min = _implied_decimal(line[36:40], 2)
    mag = decode_jma_magnitude(line[52:54])
    if None in (lat_deg, lat_min, lon_deg, lon_min, mag):
        return None
    if not (0.0 <= lat_min < 60.0 and 0.0 <= lon_min < 60.0):
        return None
    return lat_deg + lat_min / 60.0, lon_deg + lon_min / 60.0, mag


def parse_hypocenter_date(line: bytes) -> date | None:
    line = line.rstrip(b"\r\n")
    if len(line) < 13 or line[0:1] != b"J":
        return None
    try:
        return date(int(line[1:5]), int(line[5:7]), int(line[7:9]))
    except ValueError:
        return None


def iter_catalog_stream(stream: BinaryIO) -> Iterator[bytes]:
    yield from stream


def iter_input_records(paths: Iterable[Path]) -> Iterator[bytes]:
    for path in paths:
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                for name in archive.namelist():
                    if name.endswith("/"):
                        continue
                    with archive.open(name) as stream:
                        yield from iter_catalog_stream(stream)
        else:
            with path.open("rb") as stream:
                yield from iter_catalog_stream(stream)


def iter_provisional_events(path: Path) -> Iterator[tuple[date, float, float, float]]:
    """Read the normalized CSV produced by update_jma_provisional.py."""
    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                event_date = datetime.fromisoformat(row["datetime_jst"]).date()
                yield event_date, float(row["latitude"]), float(row["longitude"]), float(row["magnitude"])
            except (KeyError, TypeError, ValueError):
                continue


def iter_all_events(args: argparse.Namespace) -> Iterator[tuple[str, tuple[date, float, float, float] | None]]:
    for raw in iter_input_records(args.files):
        event = parse_hypocenter_record(raw)
        event_date = parse_hypocenter_date(raw)
        yield "fixed", ((event_date, *event) if event is not None and event_date is not None else None)
    if args.provisional_file is not None:
        for event in iter_provisional_events(args.provisional_file):
            yield "provisional", event


def find_year_files(input_dir: Path, first_year: int, last_year: int) -> list[Path]:
    found: list[Path] = []
    missing: list[int] = []
    for year in range(first_year, last_year + 1):
        matches = sorted(p for p in input_dir.glob(f"h{year}*") if p.is_file())
        if not matches:
            missing.append(year)
        found.extend(matches)
    if missing:
        raise FileNotFoundError("Missing catalog year(s): " + ", ".join(map(str, missing)))
    return found


def choose_valid(mc: np.ndarray, rule: str, fixed_mc: float) -> np.ndarray:
    finite = np.isfinite(mc)
    if rule == "ignore":
        return finite
    if rule == "below":
        return finite & (mc < fixed_mc)
    if rule == "at-or-below":
        return finite & (mc <= fixed_mc)
    if rule == "at-or-above":
        return finite & (mc >= fixed_mc)
    raise ValueError(rule)


def grid_label(grid: float) -> str:
    """Return a filesystem-friendly grid label: 0.2 -> 0p2deg, 1.0 -> 1deg."""
    text = f"{grid:.10g}".replace(".", "p")
    return f"{text}deg"


def analyze(args: argparse.Namespace) -> dict[str, object]:
    lat_edges = np.arange(args.lat_min, args.lat_max + args.grid / 2, args.grid)
    lon_edges = np.arange(args.lon_min, args.lon_max + args.grid / 2, args.grid)
    nlat, nlon = len(lat_edges) - 1, len(lon_edges) - 1
    shape = (nlat, nlon)
    fmd = np.zeros((nlat, nlon, len(MAG_CENTERS)), dtype=np.uint32)
    n_above = np.zeros(shape, dtype=np.uint32)
    sum_above = np.zeros(shape, dtype=np.float64)

    stats = {"records": 0, "parsed_j": 0, "provisional": 0,
             "in_region": 0, "used_in_fmd": 0}
    for source, event in iter_all_events(args):
        stats["records"] += 1
        if event is None:
            continue
        if source == "fixed":
            stats["parsed_j"] += 1
        else:
            stats["provisional"] += 1
        event_date, lat, lon, mag = event
        if args.date_start is not None and event_date < args.date_start:
            continue
        if args.date_end is not None and event_date > args.date_end:
            continue
        if not (args.lat_min <= lat <= args.lat_max and args.lon_min <= lon <= args.lon_max):
            continue
        # Include an event exactly on the north/east boundary in the final cell.
        ilat = min(int(math.floor((lat - args.lat_min) / args.grid + 1e-10)), nlat - 1)
        ilon = min(int(math.floor((lon - args.lon_min) / args.grid + 1e-10)), nlon - 1)
        stats["in_region"] += 1
        mt = int(round(mag * 10))
        if MAG_MIN_TENTH <= mt <= MAG_MAX_TENTH:
            fmd[ilat, ilon, mt - MAG_MIN_TENTH] += 1
            stats["used_in_fmd"] += 1
        if mag + 1e-9 >= args.fixed_mc:
            n_above[ilat, ilon] += 1
            sum_above[ilat, ilon] += mag

    total = fmd.sum(axis=2)
    mc = np.full(shape, np.nan)
    occupied = total > 0
    # np.argmax selects the lowest magnitude if multiple bins share the mode.
    mc[occupied] = MAG_CENTERS[np.argmax(fmd[occupied], axis=1)]

    mean_mag = np.divide(sum_above, n_above, out=np.full(shape, np.nan), where=n_above > 0)
    denominator = mean_mag - (args.fixed_mc - args.bin_width / 2.0)
    b_value = np.divide(math.log10(math.e), denominator,
                        out=np.full(shape, np.nan), where=denominator > 0)
    valid = choose_valid(mc, args.mc_rule, args.fixed_mc) & (n_above >= args.min_events)
    b_value[~valid] = np.nan

    return {"lat_edges": lat_edges, "lon_edges": lon_edges, "fmd": fmd,
            "total": total, "mc": mc, "n_above": n_above,
            "mean_mag": mean_mag, "b_value": b_value, "valid": valid,
            "stats": stats}


def latest_provisional_date(path: Path | None) -> date | None:
    if path is None:
        return None
    latest: date | None = None
    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                value = datetime.fromisoformat(row["datetime_jst"]).date()
            except (KeyError, TypeError, ValueError):
                continue
            latest = value if latest is None or value > latest else latest
    return latest


def overall_fmd_summary(result: dict[str, object], fixed_mc: float,
                        bin_width: float) -> dict[str, object]:
    counts = result["fmd"].sum(axis=(0, 1), dtype=np.uint64)
    cumulative = np.cumsum(counts[::-1], dtype=np.uint64)[::-1]
    occupied = counts > 0
    mc = float(MAG_CENTERS[np.argmax(counts)]) if occupied.any() else math.nan
    threshold_index = int(round(fixed_mc * 10)) - MAG_MIN_TENTH
    threshold_index = max(0, min(threshold_index, len(counts)))
    selected_counts = counts[threshold_index:]
    selected_mags = MAG_CENTERS[threshold_index:]
    n = int(selected_counts.sum())
    mean = (float(np.dot(selected_mags, selected_counts)) / n) if n else math.nan
    denominator = mean - (fixed_mc - bin_width / 2.0)
    b = math.log10(math.e) / denominator if denominator > 0 else math.nan
    return {"counts": counts, "cumulative": cumulative, "mc": mc,
            "n_above": n, "mean_mag": mean, "b_value": b}


def shift_months(value: date, months: int) -> date:
    """Shift a date by whole calendar months, clipping to the target month."""
    month_index = value.year * 12 + value.month - 1 + months
    year, month_zero = divmod(month_index, 12)
    month = month_zero + 1
    day = min(value.day, calendar.monthrange(year, month)[1])
    return date(year, month, day)


def period_ranges(latest: date, names: list[str]) -> list[tuple[str, date | None, date | None]]:
    starts = {
        "1month": shift_months(latest, -1) + timedelta(days=1),
        "6months": shift_months(latest, -6) + timedelta(days=1),
        "1year": shift_months(latest, -12) + timedelta(days=1),
        "5years": shift_months(latest, -60) + timedelta(days=1),
        "10years": shift_months(latest, -120) + timedelta(days=1),
        "all": None,
    }
    return [(name, starts[name], latest if name != "all" else None) for name in names]


def analyze_periods(args: argparse.Namespace,
                    periods: list[tuple[str, date | None, date | None]]) -> dict[str, dict[str, object]]:
    """Read the catalog once and accumulate every requested time window."""
    lat_edges = np.arange(args.lat_min, args.lat_max + args.grid / 2, args.grid)
    lon_edges = np.arange(args.lon_min, args.lon_max + args.grid / 2, args.grid)
    nlat, nlon = len(lat_edges) - 1, len(lon_edges) - 1
    shape = (nlat, nlon)
    states: dict[str, dict[str, object]] = {}
    for label, start, end in periods:
        states[label] = {
            "date_start": start, "date_end": end,
            "fmd": np.zeros((nlat, nlon, len(MAG_CENTERS)), dtype=np.uint32),
            "n_above": np.zeros(shape, dtype=np.uint32),
            "sum_above": np.zeros(shape, dtype=np.float64),
            "stats": {"records": 0, "parsed_j": 0, "provisional": 0,
                      "in_region": 0, "used_in_fmd": 0},
        }

    source_stats = {"records": 0, "parsed_j": 0, "provisional": 0}
    for source, event in iter_all_events(args):
        source_stats["records"] += 1
        if event is None:
            continue
        source_stats["parsed_j" if source == "fixed" else "provisional"] += 1
        event_date, lat, lon, mag = event
        if not (args.lat_min <= lat <= args.lat_max and args.lon_min <= lon <= args.lon_max):
            continue
        ilat = min(int(math.floor((lat - args.lat_min) / args.grid + 1e-10)), nlat - 1)
        ilon = min(int(math.floor((lon - args.lon_min) / args.grid + 1e-10)), nlon - 1)
        mt = int(round(mag * 10))
        for label, start, end in periods:
            if start is not None and event_date < start:
                continue
            if end is not None and event_date > end:
                continue
            state = states[label]
            state["stats"]["in_region"] += 1
            if MAG_MIN_TENTH <= mt <= MAG_MAX_TENTH:
                state["fmd"][ilat, ilon, mt - MAG_MIN_TENTH] += 1
                state["stats"]["used_in_fmd"] += 1
            if mag + 1e-9 >= args.fixed_mc:
                state["n_above"][ilat, ilon] += 1
                state["sum_above"][ilat, ilon] += mag

    results: dict[str, dict[str, object]] = {}
    for label, state in states.items():
        state["stats"].update(source_stats)
        fmd, n_above, sum_above = state["fmd"], state["n_above"], state["sum_above"]
        total = fmd.sum(axis=2)
        mc = np.full(shape, np.nan)
        occupied = total > 0
        mc[occupied] = MAG_CENTERS[np.argmax(fmd[occupied], axis=1)]
        mean_mag = np.divide(sum_above, n_above, out=np.full(shape, np.nan), where=n_above > 0)
        denominator = mean_mag - (args.fixed_mc - args.bin_width / 2.0)
        b_value = np.divide(math.log10(math.e), denominator,
                            out=np.full(shape, np.nan), where=denominator > 0)
        valid = choose_valid(mc, args.mc_rule, args.fixed_mc) & (n_above >= args.min_events)
        b_value[~valid] = np.nan
        results[label] = {"lat_edges": lat_edges, "lon_edges": lon_edges, "fmd": fmd,
                          "total": total, "mc": mc, "n_above": n_above,
                          "mean_mag": mean_mag, "b_value": b_value, "valid": valid,
                          "stats": state["stats"], "date_start": state["date_start"],
                          "date_end": state["date_end"]}
    return results


def write_overall_fmd_csv(path: Path, summary: dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["magnitude", "count", "cumulative_count"])
        for magnitude, count, cumulative in zip(
                MAG_CENTERS, summary["counts"], summary["cumulative"]):
            if count or cumulative:
                writer.writerow([f"{magnitude:.1f}", int(count), int(cumulative)])


def plot_overall_fmd(path: Path, summary: dict[str, object], args: argparse.Namespace) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required for FMD output") from exc
    counts, cumulative = summary["counts"], summary["cumulative"]
    positive_count = counts > 0
    positive_cumulative = cumulative > 0
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.bar(MAG_CENTERS[positive_count], counts[positive_count], width=args.bin_width * 0.85,
           color="0.70", edgecolor="0.25", linewidth=0.4, label="Incremental count")
    ax.semilogy(MAG_CENTERS[positive_cumulative], cumulative[positive_cumulative],
                "o", ms=3, color="black", label="Cumulative count")
    b = summary["b_value"]
    threshold_idx = int(round(args.fixed_mc * 10)) - MAG_MIN_TENTH
    if np.isfinite(b) and 0 <= threshold_idx < len(cumulative) and cumulative[threshold_idx] > 0:
        last_observed = int(np.flatnonzero(counts)[-1])
        x = MAG_CENTERS[threshold_idx:last_observed + 1]
        prediction = cumulative[threshold_idx] * 10 ** (-b * (x - args.fixed_mc))
        ax.semilogy(x, prediction, color="red", linewidth=1.5,
                    label=f"MLE fit: b={b:.3f} (Mc fixed={args.fixed_mc:.1f})")
    ax.axvline(args.fixed_mc, color="red", linestyle="--", linewidth=1)
    if np.isfinite(summary["mc"]):
        ax.axvline(summary["mc"], color="royalblue", linestyle=":", linewidth=1.2,
                   label=f"Maximum Curvature Mc={summary['mc']:.1f}")
    ax.set(xlabel="Magnitude", ylabel="Number of earthquakes (log scale)",
           title=f"Overall frequency-magnitude distribution ({args.period_label})")
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    ax.legend(fontsize=8)
    if positive_count.any():
        observed = MAG_CENTERS[positive_count]
        ax.set_xlim(observed.min() - 0.2, observed.max() + 0.3)
    ax.set_ylim(bottom=0.8)
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)


def write_grid_csv(path: Path, result: dict[str, object], args: argparse.Namespace) -> None:
    lat_edges = result["lat_edges"]
    lon_edges = result["lon_edges"]
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["lat_min", "lat_max", "lon_min", "lon_max", "n_total",
                         "mc_maximum_curvature", "n_m_ge_2_5", "mean_m_ge_2_5",
                         "b_value", "is_displayed"])
        for i in range(len(lat_edges) - 1):
            for j in range(len(lon_edges) - 1):
                mc = result["mc"][i, j]
                mean = result["mean_mag"][i, j]
                b = result["b_value"][i, j]
                writer.writerow([f"{lat_edges[i]:.1f}", f"{lat_edges[i+1]:.1f}",
                                 f"{lon_edges[j]:.1f}", f"{lon_edges[j+1]:.1f}",
                                 int(result["total"][i, j]),
                                 "" if np.isnan(mc) else f"{mc:.1f}",
                                 int(result["n_above"][i, j]),
                                 "" if np.isnan(mean) else f"{mean:.5f}",
                                 "" if np.isnan(b) else f"{b:.5f}",
                                 int(result["valid"][i, j])])


def write_fmd_csv(path: Path, result: dict[str, object]) -> None:
    lat_edges, lon_edges, fmd = result["lat_edges"], result["lon_edges"], result["fmd"]
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["lat_min", "lon_min", "magnitude", "count", "cumulative_count"])
        for i, j in np.argwhere(fmd.sum(axis=2) > 0):
            counts = fmd[i, j]
            cumulative = np.cumsum(counts[::-1], dtype=np.uint64)[::-1]
            for k in np.flatnonzero(counts):
                writer.writerow([f"{lat_edges[i]:.1f}", f"{lon_edges[j]:.1f}",
                                 f"{MAG_CENTERS[k]:.1f}", int(counts[k]), int(cumulative[k])])


def plot_map(path: Path, result: dict[str, object], args: argparse.Namespace) -> None:
    try:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from cartopy.mpl.gridliner import LatitudeFormatter, LongitudeFormatter
        from matplotlib.colors import Normalize
    except ImportError as exc:
        raise SystemExit(
            "matplotlib and cartopy are required for PNG output; "
            "run: conda env update -f environment.yml --prune"
        ) from exc
    values = np.ma.masked_invalid(result["b_value"])
    cmap = plt.colormaps["coolwarm"].copy()  # low=blue, high=red
    cmap.set_bad("0.70")
    finite = values.compressed()
    if args.b_min is None:
        vmin = float(np.percentile(finite, 2)) if finite.size else 0.0
    else:
        vmin = args.b_min
    if args.b_max is None:
        vmax = float(np.percentile(finite, 98)) if finite.size else 2.0
    else:
        vmax = args.b_max
    if vmax <= vmin:
        vmax = vmin + 1.0
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(10, 9), constrained_layout=True,
                           subplot_kw={"projection": projection})
    ax.set_facecolor("white")
    ax.set_extent([args.lon_min, args.lon_max, args.lat_min, args.lat_max], crs=projection)
    mesh = ax.pcolormesh(result["lon_edges"], result["lat_edges"], values,
                         cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax), shading="flat",
                         transform=projection, zorder=1)
    # Draw reference geography above the colored cells so coastlines stay visible.
    ax.coastlines(resolution=args.map_resolution, linewidth=0.65,
                  color="black", zorder=3)
    ax.add_feature(cfeature.BORDERS.with_scale(args.map_resolution),
                   linewidth=0.45, edgecolor="black", zorder=3)
    ax.gridlines(crs=projection, draw_labels=False, linewidth=0.35,
                 color="black", alpha=0.35, linestyle=":", zorder=2)
    lon_ticks = np.arange(math.ceil(args.lon_min / 5) * 5, args.lon_max + 0.1, 5)
    lat_ticks = np.arange(math.ceil(args.lat_min / 5) * 5, args.lat_max + 0.1, 5)
    ax.set_xticks(lon_ticks, crs=projection)
    ax.set_yticks(lat_ticks, crs=projection)
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol="°"))
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol="°"))
    ax.tick_params(axis="both", labelsize=9, pad=4)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"JMA b-value map ({args.period_label}; M >= {args.fixed_mc:.1f}; "
                 f"N >= {args.min_events})")
    fig.colorbar(mesh, ax=ax, label="b value")
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir", type=Path, help="Directory containing h1998 ... h2023")
    p.add_argument("--output-dir", type=Path, default=Path("outputs"))
    p.add_argument("--first-year", type=int, default=1998)
    p.add_argument("--last-year", type=int, default=2023)
    p.add_argument("--provisional-file", type=Path, default=None,
                   help="CSV made by update_jma_provisional.py; defaults to INPUT_DIR/provisional_daily.csv if present")
    p.add_argument("--lat-min", type=float, default=20.0)
    p.add_argument("--lat-max", type=float, default=50.0)
    p.add_argument("--lon-min", type=float, default=120.0)
    p.add_argument("--lon-max", type=float, default=150.0)
    p.add_argument("--grid", type=float, nargs="+", default=None, metavar="DEGREES",
                   help="Override period-specific defaults for every period")
    p.add_argument("--short-grid", type=float, default=1.0,
                   help="Grid for 1month, 6months and 1year (default: 1.0)")
    p.add_argument("--long-grid", type=float, default=0.5,
                   help="Grid for 5years, 10years and all (default: 0.5)")
    p.add_argument("--calendar-years", type=int, nargs=2, metavar=("START", "END"),
                   help="Produce one map per calendar year, inclusive")
    p.add_argument("--year-grid", type=float, default=0.1,
                   help="Grid width for --calendar-years (default: 0.1)")
    p.add_argument("--bin-width", type=float, default=0.1,
                   help="Magnitude bin width; this implementation requires 0.1")
    p.add_argument("--fixed-mc", type=float, default=2.5,
                   help="Fixed completeness threshold used in the b-value MLE")
    p.add_argument("--mc-rule", choices=("ignore", "below", "at-or-below", "at-or-above"),
                   default="at-or-below",
                   help="Diagnostic-Mc gate; default requires estimated Mc <= fixed Mc")
    p.add_argument("--min-events", type=int, default=100,
                   help="Minimum number of M >= fixed-Mc events required for display")
    p.add_argument("--periods", nargs="+",
                   choices=("1month", "6months", "1year", "5years", "10years", "all"),
                   default=["1month", "6months", "1year", "5years", "10years", "all"],
                   help="Time windows to produce (default: all six)")
    p.add_argument("--save-fmd", action="store_true", help="Also save non-zero FMD bins as CSV")
    p.add_argument("--b-min", type=float, default=0.5,
                   help="Fixed lower color-scale limit (default: 0.5)")
    p.add_argument("--b-max", type=float, default=1.4,
                   help="Fixed upper color-scale limit (default: 1.4)")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--map-resolution", choices=("110m", "50m", "10m"), default="50m",
                   help="Natural Earth coastline resolution used by Cartopy")
    return p


def write_result_bundle(args: argparse.Namespace, result: dict[str, object], suffix: str) -> None:
    grid_csv = args.output_dir / f"bvalue_grid_{suffix}.csv"
    map_png = args.output_dir / f"bvalue_map_{suffix}.png"
    write_grid_csv(grid_csv, result, args)
    plot_map(map_png, result, args)
    if args.save_fmd:
        write_fmd_csv(args.output_dir / f"fmd_by_cell_{suffix}.csv", result)
    summary = overall_fmd_summary(result, args.fixed_mc, args.bin_width)
    overall_csv = args.output_dir / f"overall_fmd_{suffix}.csv"
    overall_png = args.output_dir / f"overall_fmd_{suffix}.png"
    write_overall_fmd_csv(overall_csv, summary)
    plot_overall_fmd(overall_png, summary, args)

    metadata = vars(args).copy()
    metadata["input_dir"] = str(metadata["input_dir"])
    metadata["output_dir"] = str(metadata["output_dir"])
    metadata["provisional_file"] = (str(metadata["provisional_file"])
                                      if metadata["provisional_file"] else None)
    metadata["date_start"] = (metadata["date_start"].isoformat()
                                if metadata["date_start"] else None)
    metadata["date_end"] = (metadata["date_end"].isoformat()
                              if metadata["date_end"] else None)
    metadata["files"] = [str(p) for p in metadata["files"]]
    metadata["run_stats"] = result["stats"]
    metadata["overall_fmd"] = {"maximum_curvature_mc": summary["mc"],
                               "fixed_mc": args.fixed_mc,
                               "n_m_ge_fixed_mc": summary["n_above"],
                               "mean_m_ge_fixed_mc": summary["mean_mag"],
                               "b_value": summary["b_value"]}
    (args.output_dir / f"run_metadata_{suffix}.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {map_png}")
    print(f"Overall FMD: Mc={summary['mc']:.1f}, b={summary['b_value']:.4f}, "
          f"N(M>={args.fixed_mc:.1f})={summary['n_above']}")


def main() -> int:
    args = build_parser().parse_args()
    if not math.isclose(args.bin_width, 0.1):
        raise SystemExit("--bin-width currently must be 0.1")
    if args.b_min is not None and args.b_max is not None and args.b_max <= args.b_min:
        raise SystemExit("--b-max must be greater than --b-min")
    if args.provisional_file is None:
        candidate = args.input_dir / "provisional_daily.csv"
        if candidate.is_file():
            args.provisional_file = candidate
    elif not args.provisional_file.is_file():
        raise SystemExit(f"Provisional file not found: {args.provisional_file}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.calendar_years is not None:
        start_year, end_year = args.calendar_years
        if end_year < start_year:
            raise SystemExit("--calendar-years END must be >= START")
        if args.year_grid <= 0:
            raise SystemExit("--year-grid must be positive")
        for year in range(start_year, end_year + 1):
            run_args = copy.copy(args)
            run_args.grid = args.year_grid
            run_args.date_start = date(year, 1, 1)
            run_args.date_end = date(year, 12, 31)
            run_args.period_label = str(year)
            run_args.files = sorted(
                path for path in args.input_dir.glob(f"h{year}*") if path.is_file())
            if year <= args.last_year and not run_args.files:
                raise SystemExit(f"Fixed catalog for {year} was not found in {args.input_dir}")
            if year > args.last_year and run_args.provisional_file is None:
                raise SystemExit(f"{year} requires provisional_daily.csv")
            if year <= args.last_year:
                # Avoid rescanning the 2024+ provisional file for every historical year.
                run_args.provisional_file = None
            print(f"Calendar year {year}, grid {args.year_grid:g} deg")
            result = analyze(run_args)
            suffix = f"{year}_{grid_label(args.year_grid)}"
            write_result_bundle(run_args, result, suffix)
        return 0

    args.files = find_year_files(args.input_dir, args.first_year, args.last_year)
    latest = latest_provisional_date(args.provisional_file) or date(args.last_year, 12, 31)
    periods = period_ranges(latest, list(dict.fromkeys(args.periods)))
    if args.grid is None:
        short_names = {"1month", "6months", "1year"}
        short_periods = [period for period in periods if period[0] in short_names]
        long_periods = [period for period in periods if period[0] not in short_names]
        jobs = []
        if short_periods:
            jobs.append((args.short_grid, short_periods))
        if long_periods:
            jobs.append((args.long_grid, long_periods))
    else:
        jobs = [(grid, periods) for grid in args.grid]

    for grid, job_periods in jobs:
        if grid <= 0:
            raise SystemExit("Every --grid value must be positive")
        lat_cells = (args.lat_max - args.lat_min) / grid
        lon_cells = (args.lon_max - args.lon_min) / grid
        if not (math.isclose(lat_cells, round(lat_cells), abs_tol=1e-8)
                and math.isclose(lon_cells, round(lon_cells), abs_tol=1e-8)):
            raise SystemExit(f"Grid {grid:g} does not divide the requested map extent exactly")

        run_args = copy.copy(args)
        run_args.grid = grid
        results = analyze_periods(run_args, job_periods)
        for period_label, result in results.items():
            run_args.period_label = period_label
            run_args.date_start = result["date_start"]
            run_args.date_end = result["date_end"]
            dated_period = (period_label if period_label == "all" else
                            f"{period_label}_{result['date_start']:%Y%m%d}_{result['date_end']:%Y%m%d}")
            suffix = f"{dated_period}_{grid_label(grid)}"
            grid_csv = args.output_dir / f"bvalue_grid_{suffix}.csv"
            map_png = args.output_dir / f"bvalue_map_{suffix}.png"
            write_grid_csv(grid_csv, result, run_args)
            plot_map(map_png, result, run_args)
            if args.save_fmd:
                write_fmd_csv(args.output_dir / f"fmd_by_cell_{suffix}.csv", result)
            summary = overall_fmd_summary(result, run_args.fixed_mc, run_args.bin_width)
            overall_csv = args.output_dir / f"overall_fmd_{suffix}.csv"
            overall_png = args.output_dir / f"overall_fmd_{suffix}.png"
            write_overall_fmd_csv(overall_csv, summary)
            plot_overall_fmd(overall_png, summary, run_args)

            metadata = vars(run_args).copy()
            metadata["input_dir"] = str(metadata["input_dir"])
            metadata["output_dir"] = str(metadata["output_dir"])
            metadata["provisional_file"] = (str(metadata["provisional_file"])
                                              if metadata["provisional_file"] else None)
            metadata["date_start"] = (metadata["date_start"].isoformat()
                                        if metadata["date_start"] else None)
            metadata["date_end"] = (metadata["date_end"].isoformat()
                                      if metadata["date_end"] else None)
            metadata["files"] = [str(p) for p in metadata["files"]]
            metadata["run_stats"] = result["stats"]
            metadata["overall_fmd"] = {"maximum_curvature_mc": summary["mc"],
                                       "fixed_mc": run_args.fixed_mc,
                                       "n_m_ge_fixed_mc": summary["n_above"],
                                       "mean_m_ge_fixed_mc": summary["mean_mag"],
                                       "b_value": summary["b_value"]}
            (args.output_dir / f"run_metadata_{suffix}.json").write_text(
                json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            print(f"{period_label}, grid {grid:g}: {json.dumps(result['stats'], ensure_ascii=False)}")
            print(f"Wrote {map_png}")
            print(f"Overall FMD: Mc={summary['mc']:.1f}, b={summary['b_value']:.4f}, "
                  f"N(M>={run_args.fixed_mc:.1f})={summary['n_above']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
