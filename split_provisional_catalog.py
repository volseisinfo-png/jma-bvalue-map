#!/usr/bin/env python3
"""Split provisional_daily.csv into GitHub-friendly yearly CSV files."""
from __future__ import annotations
import argparse,csv
from pathlib import Path

def main():
 p=argparse.ArgumentParser(); p.add_argument("source",type=Path); p.add_argument("output_dir",type=Path); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
 handles={}; writers={}; counts={}
 try:
  with a.source.open(encoding="utf-8",newline="") as src:
   reader=csv.DictReader(src)
   if not reader.fieldnames: raise SystemExit("CSV header was not found")
   for row in reader:
    try: year=int(row["datetime_jst"][:4])
    except (KeyError,TypeError,ValueError): continue
    if year not in writers:
     path=a.output_dir/f"provisional_{year}.csv"; handles[year]=path.open("w",encoding="utf-8",newline=""); writers[year]=csv.DictWriter(handles[year],fieldnames=reader.fieldnames); writers[year].writeheader(); counts[year]=0
    writers[year].writerow(row); counts[year]+=1
 finally:
  for handle in handles.values(): handle.close()
 for year in sorted(counts): print(f"{year}: {counts[year]} events")
if __name__=="__main__": main()
