"""
    Example run FMBV and export data as .csv.
"""

from usa.analysis.base import GE_US
from usa.analysis.fmbv import fmbv_v0 as fmbv

from usa.interface.fmbv import fmbv_on_single, fmbv_on_list, us_bundle_from_paths

import pandas as pd

def run():
    0

"""
    Load and process data.
"""
PD_PATH = "examples/fmbv/data/wl1_12_dp.nii.gz"
SEG_PATH = "examples/fmbv/data/wl1_12_seg.nii.gz"
KRETZ_PATH = "examples/fmbv/data/wl1_12.vol"

us = us_bundle_from_paths(
    kretz_path=KRETZ_PATH,
    pd_path=PD_PATH,
    seg_path=SEG_PATH
)

# f = fmbv(us, verbose=True)
# f.global_method()
# f.depth_correction_method()



