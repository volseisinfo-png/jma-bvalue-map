#!/usr/bin/env python3
"""Backtest the research foreshock monitor and write web JSON."""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict,deque
from datetime import datetime,timedelta
from pathlib import Path
from generate_foreshock_alerts import (JST,NATIONAL,REGIONAL,Event,distance,inside,limits,read_events)

def is_aftershock(event,previous,use_magnitude_rule):
 age=(event.time-previous.time).total_seconds()/86400
 if age<=0:return False
 # Small-event filtering follows Hirose et al. Eq. (3). Target-mainshock
 # declustering uses Eqs. (1)-(2), but its preceding mainshock must still be
 # larger; a later larger event cannot be an aftershock of a smaller event.
 if use_magnitude_rule:
  if event.mag>=previous.mag-1.:return False
 elif event.mag>=previous.mag:
  return False
 km,days=limits(previous.mag)
 return age<=days and distance(previous,event)<=km

def decluster(events,use_magnitude_rule,min_candidate_magnitude=2.):
 """Sequentially test against original preceding events, as in equations 1-3."""
 kept=[];removed=[]; buckets=defaultdict(deque); size=.25
 max_km=max((limits(e.mag)[0] for e in events),default=0.); max_days=max((limits(e.mag)[1] for e in events),default=0.); radius=max(1,math.ceil(max_km/17.))
 for event in events:
  rejected=False; bi,bj=math.floor(event.lat/size),math.floor(event.lon/size)
  if event.mag>=min_candidate_magnitude:
   for i in range(bi-radius,bi+radius+1):
    for j in range(bj-radius,bj+radius+1):
     bucket=buckets.get((i,j))
     if not bucket:continue
     while bucket and (event.time-bucket[0].time).total_seconds()/86400>max_days:bucket.popleft()
     for previous in reversed(bucket):
      if is_aftershock(event,previous,use_magnitude_rule):rejected=True;break
     if rejected:break
    if rejected:break
   (removed if rejected else kept).append(event)
  buckets[(bi,bj)].append(event)
 return kept,removed

def cell_indices(event,p):
 step=p.d/2; half=p.d/2
 i0=math.ceil((event.lat-half-(p.lat_min+half))/step-1e-10);i1=math.floor((event.lat+half-(p.lat_min+half))/step+1e-10)
 j0=math.ceil((event.lon-half-(p.lon_min+half))/step-1e-10);j1=math.floor((event.lon+half-(p.lon_min+half))/step+1e-10)
 max_i=math.floor((p.lat_max-p.lat_min-p.d)/step+1e-10);max_j=math.floor((p.lon_max-p.lon_min-p.d)/step+1e-10)
 for i in range(max(0,i0),min(max_i,i1)+1):
  for j in range(max(0,j0),min(max_j,j1)+1):yield round(p.lat_min+half+i*step,10),round(p.lon_min+half+j*step,10)

def alarms_for_profile(events,p):
 cells=defaultdict(list)
 for e in events:
  if e.mag<p.mf or not inside(e.lat,e.lon,p):continue
  if p.depth_min is not None and (e.depth is None or not p.depth_min<=e.depth<=p.depth_max):continue
  for center in cell_indices(e,p):cells[center].append(e)
 alarms=[]
 for (cy,cx),cell in cells.items():
  window=deque()
  for event in cell:
   before=len(window)
   while window and window[0].time<event.time-timedelta(days=p.tf):window.popleft()
   before=len(window);window.append(event)
   if before<p.nf and len(window)>=p.nf:
    alarms.append({"profile_id":p.id,"profile_label":p.label,"source":p.source,"parameters":{"D_degrees":p.d,"Mf_min":p.mf,"Tf_days":p.tf,"Nf":p.nf,"Ta_days":p.ta,"target_M_min":p.target_m},"center":{"latitude":cy,"longitude":cx},"bounds":{"south":cy-p.d/2,"north":cy+p.d/2,"west":cx-p.d/2,"east":cx+p.d/2},"triggered_at":event.time,"expires_at":event.time+timedelta(days=p.ta),"events":list(window)[-p.nf:]})
 return alarms

def profile_for_target(e):
 for p in REGIONAL:
  if inside(e.lat,e.lon,p) and (p.depth_min is None or e.depth is not None and p.depth_min<=e.depth<=p.depth_max):return p
 return NATIONAL if inside(e.lat,e.lon,NATIONAL) else None

def event_json(e):return {"datetime_jst":e.time.isoformat(),"latitude":e.lat,"longitude":e.lon,"depth_km":e.depth,"magnitude":e.mag,"region":e.region}
def alarm_json(a):
 out={k:v for k,v in a.items() if k!="events"};out["triggered_at"]=a["triggered_at"].isoformat();out["expires_at"]=a["expires_at"].isoformat();out["foreshocks"]=[event_json(e) for e in a["events"]];return out

def group_alarm_earthquakes(alarms):
 """Count one alarm earthquake once even when it alarms overlapping grids."""
 grouped={}
 for alarm in alarms:
  ae=alarm["events"][-1]
  key=(alarm["profile_id"],ae.time,ae.lat,ae.lon,ae.mag,ae.depth)
  if key not in grouped:
   grouped[key]={**alarm,"center":None,"cell_bounds":[],"cell_centers":[],"raw_cell_count":0,"events":[]}
  item=grouped[key];item["cell_bounds"].append(alarm["bounds"]);item["cell_centers"].append(alarm["center"]);item["raw_cell_count"]+=1
  seen={(e.time,e.lat,e.lon,e.mag,e.depth) for e in item["events"]}
  for e in alarm["events"]:
   event_key=(e.time,e.lat,e.lon,e.mag,e.depth)
   if event_key not in seen:item["events"].append(e);seen.add(event_key)
 for item in grouped.values():
  item["events"].sort(key=lambda e:e.time); bounds=item["cell_bounds"]
  item["bounds"]={"south":min(x["south"] for x in bounds),"north":max(x["north"] for x in bounds),"west":min(x["west"] for x in bounds),"east":max(x["east"] for x in bounds)}
 return sorted(grouped.values(),key=lambda x:x["triggered_at"])

def target_in_alarm(target,alarm):
 return any(b["south"]<=target.lat<=b["north"] and b["west"]<=target.lon<=b["east"] for b in alarm["cell_bounds"])

def main():
 ap=argparse.ArgumentParser();ap.add_argument("catalogs",type=Path,nargs="+");ap.add_argument("--start",default="2024-01-01");ap.add_argument("--output",type=Path,default=Path("web/data/alerts/backtest.json"));a=ap.parse_args();start=datetime.fromisoformat(a.start).replace(tzinfo=JST)
 events=read_events(a.catalogs);foreshock_catalog,small_removed=decluster(events,True);target_catalog,target_aftershocks=decluster(events,False,5.)
 alarms=[]
 for p in REGIONAL:alarms+=alarms_for_profile(foreshock_catalog,p)
 national=alarms_for_profile(foreshock_catalog,NATIONAL);alarms += [x for x in national if not any(inside(x["center"]["latitude"],x["center"]["longitude"],p) for p in REGIONAL)];alarms=[x for x in alarms if x["triggered_at"]>=start]
 alarm_earthquakes=group_alarm_earthquakes(alarms)
 targets=[]
 for e in target_catalog:
  p=profile_for_target(e)
  if p and e.time>=start and e.mag>=p.target_m:targets.append((e,p))
 matched_alarm_ids=set();hits=[];misses=[]
 for target,p in targets:
  matches=[]
  for i,alarm in enumerate(alarm_earthquakes):
   max_foreshock=max(e.mag for e in alarm["events"])
   if alarm["profile_id"]==p.id and alarm["triggered_at"]<target.time<=alarm["expires_at"] and target.mag>max_foreshock and target_in_alarm(target,alarm):
    matched_alarm_ids.add(i);matches.append(alarm_json(alarm))
  item={"target":event_json(target),"profile_id":p.id,"profile_label":p.label,"matching_alarms":matches}
  (hits if matches else misses).append(item)
 false_alarms=[alarm_json(x) for i,x in enumerate(alarm_earthquakes) if i not in matched_alarm_ids]
 for values in (hits,misses,false_alarms):values.sort(key=lambda x:x.get("target",{}).get("datetime_jst",x.get("triggered_at","")),reverse=True)
 target_count=len(hits)+len(misses); matched_alarm_count=len(matched_alarm_ids); alarm_count=len(alarm_earthquakes)
 out={"generated_at":datetime.now(JST).isoformat(),"period_start":start.isoformat(),"period_end":events[-1].time.isoformat() if events else None,"definitions":{"hit":"対象本震がアラーム地震に対応するいずれかのセル内で警報期間中に発生","miss":"対象本震に対応する事前判定なし","false_alarm":"警報時空間内に対象本震が続かなかったアラーム地震","alarm_earthquake":"各系列で条件を満たすに至ったNf個目の地震（重複格子でも同じ地震は1件）","alarm_rate":"予測された対象本震数 / 全対象本震数","truth_rate":"適中したアラーム地震数 / 全アラーム地震数"},"counts":{"hits":len(hits),"misses":len(misses),"targets":target_count,"successful_alarm_earthquakes":matched_alarm_count,"false_alarm_earthquakes":len(false_alarms),"alarm_earthquakes":alarm_count,"alarms_raw_cells":len(alarms),"small_aftershocks_removed":len(small_removed),"target_aftershocks_removed":len(target_aftershocks)},"rates":{"alarm_rate":len(hits)/target_count if target_count else None,"truth_rate":matched_alarm_count/alarm_count if alarm_count else None},"hits":hits,"misses":misses,"false_alarms":false_alarms};a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");print(json.dumps({**out["counts"],**out["rates"]},ensure_ascii=False))
if __name__=="__main__":main()
