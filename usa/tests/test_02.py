"""
    Make sure that FMBV is refactored properly. 
"""

from usa.analysis.base import GE_US, resample_GE
from usa.analysis.fmbv import fmbv_v0 as fmbv

# hardcoded values from running on original code.
GLOBAL = 15.600886685145513
DC = 20.024242546299877
VOL = 35189.45134522546

PD_PATH = "./usa/data/wl1_12_dp.nii.gz"
SEG_PATH = "./usa/data/wl1_12_seg.nii.gz"
KRETZ_PATH = "./usa/data/wl1_12.vol"

us = GE_US(minimise_memory=False, verbose=False)
us.load_kretz(KRETZ_PATH)
us.load_pd(PD_PATH)
us.load_seg(SEG_PATH)
us.set_distance()

f = fmbv(us, verbose=True)
f.global_method()
f.depth_correction_method()

assert f.global_fmbv_2 == GLOBAL
assert f.depth_corrected_fmbv == DC
assert f.seg_volume == VOL
print(f.global_fmbv_2)
print(f.depth_corrected_fmbv)
print(f.seg_volume)