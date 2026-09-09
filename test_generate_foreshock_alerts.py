import unittest
from datetime import datetime,timedelta
from generate_foreshock_alerts import Event,JST,limits,current_candidates
from generate_foreshock_backtest import group_alarm_earthquakes,is_aftershock
class AlertTest(unittest.TestCase):
 def test_paper_examples(self):
  self.assertEqual(round(limits(8)[0]),158); self.assertEqual(round(limits(8)[1]),557); self.assertEqual(round(limits(7)[1]),123); self.assertEqual(round(limits(6)[1]),27)
 def test_small_aftershock_filter(self):
  t=datetime(2026,9,1,tzinfo=JST); main=Event(t,35.,140.,6.,10.); small=Event(t+timedelta(days=1),35.01,140.01,4.9,10.); equal=Event(t+timedelta(days=1),35.01,140.01,5.,10.)
  kept,removed=current_candidates([main,small,equal],t+timedelta(days=2)); self.assertIn(small,removed); self.assertIn(equal,kept)
 def test_later_larger_event_is_not_target_aftershock(self):
  t=datetime(2025,11,9,7,33,tzinfo=JST); previous=Event(t,39.4,143.5,5.1,10.); target=Event(t+timedelta(hours=10),39.4,143.5,6.9,16.)
  self.assertFalse(is_aftershock(target,previous,False))
 def test_later_smaller_event_can_be_target_aftershock(self):
  t=datetime(2025,11,9,7,33,tzinfo=JST); previous=Event(t,39.4,143.5,6.9,10.); target=Event(t+timedelta(hours=10),39.4,143.5,5.1,16.)
  self.assertTrue(is_aftershock(target,previous,False))
 def test_same_alarm_earthquake_in_overlapping_cells_is_counted_once(self):
  t=datetime(2025,11,9,7,15,tzinfo=JST); event=Event(t,39.4,143.5,5.9,10.)
  base={"profile_id":"off_iwate","profile_label":"岩手県沖","source":"test","parameters":{},"triggered_at":t,"expires_at":t+timedelta(days=4),"events":[event]}
  alarms=[{**base,"center":{"latitude":39.25,"longitude":143.25},"bounds":{"south":39.,"north":39.5,"west":143.,"east":143.5}},{**base,"center":{"latitude":39.25,"longitude":143.75},"bounds":{"south":39.,"north":39.5,"west":143.5,"east":144.}}]
  grouped=group_alarm_earthquakes(alarms)
  self.assertEqual(len(grouped),1);self.assertEqual(grouped[0]["raw_cell_count"],2)
if __name__=="__main__": unittest.main()
