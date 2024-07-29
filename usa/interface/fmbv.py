"""
    Convenient interfaces for processing FMBV.
"""

from usa.analysis.base import GE_US
from usa.analysis.fmbv import fmbv_v0 as fmbv

import numpy as np

def us_bundle_from_paths(
        kretz_path = None,
        pd_path = None,
        seg_path = None,
        **kwargs
):
    us = GE_US(**kwargs)
    us.load_kretz(kretz_path)
    us.load_pd(pd_path)
    us.load_seg(seg_path)
    us.set_distance()

    return us

def fmbv_on_single(
        us : GE_US
):
    0

def fmbv_on_list(
        us_list : list
):  
    assert isinstance(us_list, list)
    assert np.all([isinstance(x, GE_US) for x in us_list])
    0