#!/usr/bin/env python3
"""Create web/manifest.json from generated b-value map images."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path


WEB = Path(__file__).resolve().parent / "web"
MAP_RE = re.compile(r"^bvalue_map_(.+)_(0p5|1)deg\.png$")
YEAR_RE = re.compile(r"^(19|20)\d{2}$")
ORDER = {"1month": 0, "6months": 1, "1year": 2, "5years": 3,
         "10years": 4, "all": 5}
LABELS = {"1month": "最近1か月", "6months": "最近6か月", "1year": "最近1年",
          "5years": "最近5年", "10years": "最近10年", "all": "全期間"}


def main() -> None:
    entries = []
    for image in sorted((WEB / "data").rglob("bvalue_map_*.png")):
        match = MAP_RE.match(image.name)
        if not match:
            continue
        identity, grid_code = match.groups()
        parts = identity.split("_")
        period = parts[0]
        is_year = bool(YEAR_RE.fullmatch(period))
        fmd = image.with_name(image.name.replace("bvalue_map_", "overall_fmd_", 1))
        entries.append({
            "id": f"{identity}_{grid_code}deg",
            "group": "yearly" if is_year else "recent",
            "label": f"{period}年" if is_year else LABELS.get(period, period),
            "period": period,
            "date_start": parts[1] if len(parts) >= 3 else None,
            "date_end": parts[2] if len(parts) >= 3 else None,
            "grid_degrees": 0.5 if grid_code == "0p5" else 1.0,
            "map": image.relative_to(WEB).as_posix(),
            "fmd": fmd.relative_to(WEB).as_posix() if fmd.exists() else None,
        })
    entries.sort(key=lambda item: (0 if item["group"] == "recent" else 1,
                                   ORDER.get(item["period"], 99) if item["group"] == "recent"
                                   else -int(item["period"])))
    manifest = {"updated_at": datetime.now(timezone.utc).isoformat(), "entries": entries}
    (WEB / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Manifest: {len(entries)} maps")


if __name__ == "__main__":
    main()
