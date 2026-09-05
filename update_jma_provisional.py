#!/usr/bin/env python3
"""Incrementally download JMA daily hypocenter lists into a provisional CSV."""

from __future__ import annotations

import argparse
import csv
import html
import re
import sys
import time
from datetime import date, datetime, timedelta
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


BASE_URL = "https://www.data.jma.go.jp/eqev/data/daily_map/"
FIELDS = ["datetime_jst", "latitude", "longitude", "depth_km", "magnitude",
          "region", "source_url", "status"]
DATE_LINK_RE = re.compile(r'href=["\'](?P<day>\d{8})\.html["\']', re.I)
EVENT_RE = re.compile(
    r"^\s*(?P<year>\d{4})\s+(?P<month>\d{1,2})\s+(?P<day>\d{1,2})\s+"
    r"(?P<hour>\d{1,2}):(?P<minute>\d{2})\s+(?P<second>\d{1,2}(?:\.\d+)?)\s+"
    r"(?P<lat_deg>\d{1,2})°\s*(?P<lat_min>\d{1,2}(?:\.\d+)?)['′]N\s+"
    r"(?P<lon_deg>\d{1,3})°\s*(?P<lon_min>\d{1,2}(?:\.\d+)?)['′]E\s+"
    r"(?P<depth>-?\d+|---)\s+(?P<mag>-?\d+(?:\.\d+)?)\s*(?P<region>.*)$"
)


class TextExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() in {"br", "p", "pre"}:
            self.parts.append("\n")

    def text(self) -> str:
        return html.unescape("".join(self.parts))


def fetch_text(url: str, timeout: float) -> str:
    request = Request(url, headers={"User-Agent": "jma-bvalue-research/1.0"})
    with urlopen(request, timeout=timeout) as response:
        raw = response.read()
        charset = response.headers.get_content_charset() or "utf-8"
    try:
        return raw.decode(charset)
    except (LookupError, UnicodeDecodeError):
        return raw.decode("utf-8", errors="replace")


def available_days(index_html: str) -> list[date]:
    return sorted({datetime.strptime(m.group("day"), "%Y%m%d").date()
                   for m in DATE_LINK_RE.finditer(index_html)})


def parse_daily_html(page_html: str, source_url: str) -> list[dict[str, str]]:
    parser = TextExtractor()
    parser.feed(page_html)
    rows: list[dict[str, str]] = []
    for line in parser.text().splitlines():
        match = EVENT_RE.match(line)
        if not match:
            continue
        item = match.groupdict()
        second = float(item["second"])
        whole_second = int(second)
        microsecond = round((second - whole_second) * 1_000_000)
        origin = datetime(int(item["year"]), int(item["month"]), int(item["day"]),
                          int(item["hour"]), int(item["minute"]), whole_second, microsecond)
        latitude = int(item["lat_deg"]) + float(item["lat_min"]) / 60.0
        longitude = int(item["lon_deg"]) + float(item["lon_min"]) / 60.0
        rows.append({
            "datetime_jst": origin.isoformat(timespec="milliseconds"),
            "latitude": f"{latitude:.6f}",
            "longitude": f"{longitude:.6f}",
            "depth_km": "" if item["depth"] == "---" else item["depth"],
            "magnitude": f"{float(item['mag']):.1f}",
            "region": item["region"].strip(),
            "source_url": source_url,
            "status": "provisional",
        })
    return rows


def latest_stored_day(path: Path) -> date | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    latest: date | None = None
    with path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                day = date.fromisoformat(row["datetime_jst"][:10])
            except (KeyError, TypeError, ValueError):
                continue
            latest = day if latest is None or day > latest else latest
    return latest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalog_dir", type=Path,
                        help="Directory containing the fixed JMA catalogs")
    parser.add_argument("--output-file", type=Path, default=None,
                        help="Default: CATALOG_DIR/provisional_daily.csv")
    parser.add_argument("--start-date", type=date.fromisoformat, default=date(2024, 1, 1))
    parser.add_argument("--end-date", type=date.fromisoformat, default=None,
                        help="Optional upper bound; default is the latest day listed by JMA")
    parser.add_argument("--refresh-last-day", action="store_true",
                        help="Re-download the latest stored day as well as newer days")
    parser.add_argument("--delay", type=float, default=0.05,
                        help="Seconds between requests")
    parser.add_argument("--timeout", type=float, default=30.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output = args.output_file or args.catalog_dir / "provisional_daily.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        index_html = fetch_text(BASE_URL, args.timeout)
    except (HTTPError, URLError, TimeoutError) as exc:
        raise SystemExit(f"Could not fetch JMA daily-list index: {exc}") from exc
    published = available_days(index_html)
    if not published:
        raise SystemExit("No dated daily-list links were found on the JMA index page")

    stored = latest_stored_day(output)
    start = args.start_date
    if stored is not None:
        start = max(start, stored if args.refresh_last_day else stored + timedelta(days=1))
    end = min(args.end_date, published[-1]) if args.end_date else published[-1]
    targets = [day for day in published if start <= day <= end]
    if not targets:
        print(f"Already up to date: {stored or 'no applicable published dates'}")
        return 0

    # Refresh mode rewrites without the selected last day, then appends its current contents.
    if args.refresh_last_day and stored is not None and output.exists():
        temp = output.with_suffix(output.suffix + ".tmp")
        with output.open("r", newline="", encoding="utf-8") as source, \
                temp.open("w", newline="", encoding="utf-8") as destination:
            reader = csv.DictReader(source)
            writer = csv.DictWriter(destination, fieldnames=FIELDS)
            writer.writeheader()
            for row in reader:
                if row.get("datetime_jst", "")[:10] != stored.isoformat():
                    writer.writerow({key: row.get(key, "") for key in FIELDS})
        temp.replace(output)

    new_file = not output.exists() or output.stat().st_size == 0
    total = 0
    with output.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
        for number, day in enumerate(targets, 1):
            url = f"{BASE_URL}{day:%Y%m%d}.html"
            try:
                rows = parse_daily_html(fetch_text(url, args.timeout), url)
            except (HTTPError, URLError, TimeoutError) as exc:
                raise SystemExit(f"Update stopped at {day}: {exc}. Re-run to resume.") from exc
            if not rows:
                raise SystemExit(f"Update stopped: no earthquake rows parsed for {day} ({url})")
            writer.writerows(rows)
            stream.flush()
            total += len(rows)
            print(f"[{number}/{len(targets)}] {day}: {len(rows)} events")
            if args.delay > 0 and number < len(targets):
                time.sleep(args.delay)
    print(f"Added {total} events through {targets[-1]} to {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
