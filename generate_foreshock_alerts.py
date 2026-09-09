#!/usr/bin/env python3
"""Generate the research-only foreshock activity monitor from yearly JMA CSVs."""
from __future__ import annotations
import argparse,csv,json,math
from dataclasses import asdict,dataclass
from datetime import datetime,timedelta,timezone
from pathlib import Path
JST=timezone(timedelta(hours=9))
@dataclass(frozen=True)
class Event: time:datetime; lat:float; lon:float; mag:float; depth:float|None; region:str=""
@dataclass(frozen=True)
class Profile:
 id:str; label:str; source:str; lat_min:float; lat_max:float; lon_min:float; lon_max:float; depth_min:float|None; depth_max:float|None; d:float; mf:float; tf:float; nf:int; ta:float; target_m:float
REGIONAL=[
 Profile("off_iwate","岩手県沖","Hirose et al. (2021)",39.,40.,143.,144.,0.,100.,.5,5.,9.,3,4.,6.),
 Profile("off_miyagi","宮城県沖","Hirose et al. (2021)",38.,39.,142.5,144.,0.,100.,.5,5.,9.,3,4.,6.),
 Profile("off_ibaraki","茨城県沖","Hirose et al. (2021)",36.,36.5,141.25,142.25,0.,100.,.5,5.,3.,2,1.,6.),
 Profile("central_honshu","中部日本","Hirose et al. (2021)",35.6,37.1,137.2,139.,0.,30.,.2,2.,1.,5,5.,5.),
 Profile("izu_islands","伊豆諸島","Hirose et al. (2021)",33.5,35.3,138.6,139.8,0.,50.,.2,3.,1.,2,4.,5.)]
NATIONAL=Profile("national","その他の全国域","Maeda (1996)",20.,50.,120.,150.,None,None,.5,5.,10.,3,5.,6.)
def read_events(paths):
 out=[]
 for path in paths:
  with path.open(encoding="utf-8-sig",newline="") as stream:
   for row in csv.DictReader(stream):
    try:
     mag=float(row["magnitude"])
     if mag<2.: continue
     t=datetime.fromisoformat(row["datetime_jst"]); t=t if t.tzinfo else t.replace(tzinfo=JST)
     depth=float(row["depth_km"]) if row.get("depth_km","").strip() else None
     out.append(Event(t,float(row["latitude"]),float(row["longitude"]),mag,depth,row.get("region","")))
    except (KeyError,TypeError,ValueError): pass
 return sorted(out,key=lambda e:e.time)
def limits(m): return 10**(.5*m-1.8),max(0.,10**((.17+.85*(m-4.))/1.3)-.3)
def distance(a,b):
 p1,p2=math.radians(a.lat),math.radians(b.lat); dp=p2-p1; dl=math.radians(b.lon-a.lon); h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
 return 12742.0176*math.asin(min(1.,math.sqrt(h)))
def current_candidates(events,as_of):
 """Apply equations 1-3 only to events needed by currently active alarms."""
 recent=[e for e in events if as_of-timedelta(days=15)<=e.time<=as_of]; retained=[]; removed=[]
 for event in recent:
  small=False
  for previous in reversed(events):
   if previous.time>=event.time: continue
   if event.mag>=previous.mag-1.: continue
   km,days=limits(previous.mag); age=(event.time-previous.time).total_seconds()/86400
   if age>days: continue
   if distance(previous,event)<=km: small=True; break
  (removed if small else retained).append(event)
 return retained,removed
def centers(lo,hi,step): return [round(lo+i*step,10) for i in range(int(math.floor((hi-lo)/step+1e-9))+1)]
def inside(lat,lon,p): return p.lat_min<=lat<=p.lat_max and p.lon_min<=lon<=p.lon_max
def detect(events,as_of,p):
 half=p.d/2; selected=[e for e in events if e.mag>=p.mf and (p.depth_min is None or e.depth is not None and p.depth_min<=e.depth<=p.depth_max)]; alarms=[]
 for cy in centers(p.lat_min+half,p.lat_max-half,half):
  for cx in centers(p.lon_min+half,p.lon_max-half,half):
   cell=[e for e in selected if abs(e.lat-cy)<=half and abs(e.lon-cx)<=half]
   for i,event in enumerate(cell):
    seq=[e for e in cell[:i+1] if e.time>=event.time-timedelta(days=p.tf)]; prev=[e for e in cell[:i] if e.time>=event.time-timedelta(days=p.tf)]
    if len(seq)>=p.nf and len(prev)<p.nf and event.time+timedelta(days=p.ta)>=as_of:
     alarms.append({"profile_id":p.id,"profile_label":p.label,"source":p.source,"parameters":{"D_degrees":p.d,"Mf_min":p.mf,"Tf_days":p.tf,"Nf":p.nf,"Ta_days":p.ta,"target_M_min":p.target_m,"depth_km":[p.depth_min,p.depth_max]},"center":{"latitude":cy,"longitude":cx},"bounds":{"south":cy-half,"north":cy+half,"west":cx-half,"east":cx+half},"triggered_at":event.time.isoformat(),"expires_at":(event.time+timedelta(days=p.ta)).isoformat(),"events":[{"datetime_jst":e.time.isoformat(),"latitude":e.lat,"longitude":e.lon,"depth_km":e.depth,"magnitude":e.mag,"region":e.region} for e in seq[-p.nf:]]})
 return alarms
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("catalogs",type=Path,nargs="+"); ap.add_argument("--output",type=Path,default=Path("web/data/alerts/current.json")); ap.add_argument("--as-of"); a=ap.parse_args(); now=datetime.fromisoformat(a.as_of) if a.as_of else datetime.now(JST); now=now if now.tzinfo else now.replace(tzinfo=JST)
 events=read_events(a.catalogs); candidates,removed=current_candidates(events,now); alarms=[]
 for p in REGIONAL: alarms+=detect(candidates,now,p)
 national=detect(candidates,now,NATIONAL); alarms += [x for x in national if not any(inside(x["center"]["latitude"],x["center"]["longitude"],p) for p in REGIONAL)]; alarms.sort(key=lambda x:x["expires_at"],reverse=True)
 out={"generated_at":datetime.now(JST).isoformat(),"as_of":now.isoformat(),"catalog_through":events[-1].time.isoformat() if events else None,"method":"Maeda (1996) nationwide + Hirose et al. (2021) regional parameters","profiles":[asdict(p) for p in REGIONAL+[NATIONAL]],"aftershock_filter":{"applied":True,"distance":"log10(L_km) <= 0.5*Mpre - 1.8","time":"log10(ta_days + 0.3) <= (0.17 + 0.85*(Mpre - 4.0))/1.3","magnitude":"Ma < Mpre - 1.0","removed_recent_events":len(removed)},"active_count":len(alarms),"active_areas":alarms}; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"); print(f"Small aftershocks removed: {len(removed)}; active cells: {len(alarms)}")
if __name__=="__main__": main()
