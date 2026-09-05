import unittest
from datetime import date

from jma_bvalue_map import (decode_jma_magnitude, grid_label,
                            parse_hypocenter_date, parse_hypocenter_record,
                            period_ranges)
from update_jma_provisional import available_days, parse_daily_html


class JmaParserTest(unittest.TestCase):
    def test_magnitudes(self):
        cases = {b"25": 2.5, b"00": 0.0, b"-1": -0.1, b"-9": -0.9,
                 b"A0": -1.0, b"A9": -1.9, b"B0": -2.0, b"C0": -3.0,
                 b"  ": None}
        for raw, expected in cases.items():
            self.assertEqual(decode_jma_magnitude(raw), expected)

    def test_fixed_columns(self):
        row = bytearray(b" " * 96)
        row[0:1] = b"J"
        row[1:9] = b"20260903"
        row[21:24] = b" 35"
        row[24:28] = b"3000"
        row[32:36] = b" 140"
        row[36:40] = b"1500"
        row[52:54] = b"25"
        self.assertEqual(parse_hypocenter_record(bytes(row)), (35.5, 140.25, 2.5))
        self.assertEqual(str(parse_hypocenter_date(bytes(row))), "2026-09-03")

    def test_non_j_and_short_records_are_ignored(self):
        self.assertIsNone(parse_hypocenter_record(b"J short"))
        self.assertIsNone(parse_hypocenter_record(b"U" + b" " * 95))

    def test_grid_filename_labels(self):
        self.assertEqual(grid_label(0.2), "0p2deg")
        self.assertEqual(grid_label(0.5), "0p5deg")
        self.assertEqual(grid_label(1.0), "1deg")

    def test_daily_list_parser(self):
        sample = """<html><body><pre>
2026  9  3 00:00 50.6  34°10.5'N 139°18.6'E   10     1.5  三宅島近海
2026  9  3 00:12 41.4  32°14.3'N 130°25.3'E   15    -0.1  熊本県天草・芦北地方
</pre></body></html>"""
        rows = parse_daily_html(sample, "https://example.test/20260903.html")
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(float(rows[0]["latitude"]), 34.175)
        self.assertAlmostEqual(float(rows[0]["longitude"]), 139.31)
        self.assertEqual(rows[1]["magnitude"], "-0.1")
        self.assertEqual(rows[1]["region"], "熊本県天草・芦北地方")

    def test_daily_index_dates(self):
        days = available_days('<a href="20240101.html">x</a><a href="20260903.html">y</a>')
        self.assertEqual([str(day) for day in days], ["2024-01-01", "2026-09-03"])

    def test_requested_period_ranges(self):
        ranges = period_ranges(date(2026, 9, 3), ["1month", "6months", "1year",
                                                       "5years", "10years", "all"])
        self.assertEqual(ranges[0], ("1month", date(2026, 8, 4), date(2026, 9, 3)))
        self.assertEqual(ranges[2], ("1year", date(2025, 9, 4), date(2026, 9, 3)))
        self.assertEqual(ranges[-1], ("all", None, None))


if __name__ == "__main__":
    unittest.main()
