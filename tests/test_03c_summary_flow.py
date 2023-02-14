"""
summary
----------------------------------

Tests for `hspf_utils` module.
"""

import shlex
import subprocess
from io import StringIO
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal


class TestDescribe(TestCase):
    def setUp(self):
        self.extract = pd.read_csv("tests/data_03c_summary_flow.csv")

    def test_summary(self):
        args = "hspf_utils summary --uci tests/03lsjr_c/1995_Scenario/3c_SJRB_1995.uci tests/03lsjr_c/1995_Scenario/3c_land.hbn"
        args = shlex.split(args)
        complete = subprocess.run(args, capture_output=True, text=True, check=True)
        out = pd.read_csv(StringIO(complete.stdout))
        assert_frame_equal(out, self.extract)
