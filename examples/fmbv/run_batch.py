import numpy as np

from usa.analysis.base import GE_US
from usa.analysis.fmbv import fmbv_v0 as fmbv

PD_PATH = "examples/fmbv/data/wl1_12_dp.nii.gz"
SEG_PATH = "examples/fmbv/data/wl1_12_seg.nii.gz"
KRETZ_PATH = "examples/fmbv/data/wl1_12.vol"

us0 = GE_US()
us0.load_pd(PD_PATH)
us0.load_seg(SEG_PATH)

f = fmbv(us0, verbose=True)


