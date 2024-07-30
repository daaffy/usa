
from usa.analysis.base import GE_US, resample_GE
from usa.analysis.fmbv import fmbv_v0 as fmbv

def test(zoom = 1.):

    PD_PATH = "wl1_12_dp.nii.gz"
    SEG_PATH = "wl1_12_seg.nii.gz"
    KRETZ_PATH = "wl1_12.vol"

    us = GE_US(minimise_memory=False)
    us.load_kretz(KRETZ_PATH)
    us.load_pd(PD_PATH)
    us.load_seg(SEG_PATH)

    if not zoom == .1:
        return resample_GE(us, 0.5)
    else:
        return us