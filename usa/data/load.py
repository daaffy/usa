
from usa.analysis.base import GE_US, resample_GE
from usa.analysis.fmbv import fmbv_v0 as fmbv

# make these relative...
PD_PATH = "/Users/jackh/Documents/PIRG/usa/examples/fmbv/data/wl2_2_dp.nii.gz"
SEG_PATH = "/Users/jackh/Documents/PIRG/usa/examples/fmbv/data/wl2_2_seg.nii.gz"
KRETZ_PATH = "/Users/jackh/Documents/PIRG/usa/examples/fmbv/data/wl2_2.vol"

def test(zoom = 1.):
    us = GE_US(minimise_memory=False)
    us.load_kretz(KRETZ_PATH)
    us.load_pd(PD_PATH)
    us.load_seg(SEG_PATH)

    if not zoom == .1:
        return resample_GE(us, 0.5)
    else:
        return us