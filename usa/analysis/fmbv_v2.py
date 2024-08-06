from usa.analysis.base import US_Method, GE_US
import numpy as np
import copy

class fmbv_v2(US_Method):

    def __init__(self, us : GE_US):
        """
            TODO Segmentation to mask pd_array prior to calculating standardisation.
        """
        self.us = us
        self.pd_array = self.us.pd_array
        self.global_fmbv = None

        self.__cumsum()

    def fmbv(self, qs):
        self.global_fmbv = np.mean(self.standardise_image(qs)[self.us.seg_array == self.us.vol_id])

    def standardise_image(self, qs):
        pd_array_std = copy.copy(self.pd_array)

        assert len(qs) == 2
        q_min, q_max = qs[0], qs[1]
        std_min, std_max = self.__get_std(q_min), self.__get_std(q_max)

        self.std_min, self.std_max = std_min, std_max

        # threshold bottom
        pd_array_std = self.pd_array - std_min
        mask_min = pd_array_std < 0 
        pd_array_std[mask_min] = 0

        # normalise
        pd_array_std = pd_array_std / (std_max - std_min)
        # clip above std_max
        mask_max = pd_array_std > 1
        pd_array_std[mask_max] = 1

        return pd_array_std

    def __interp_del(self, y, ys : list):
        assert len(ys) == 2
        return (y - ys[0])/(ys[1] - ys[0])

    def __get_std(self, q):
        assert q >= 0 and q < 1

        i_ = 0
        if (self.cdf[0] >= q):
            return 0
        
        for i in range(len(self.cdf)-1):
            if (self.cdf[i] <= q) and (self.cdf[i+1] > q):
                i_ = i
                break
                
        # need to do some checks here...
        i_true = i_ + self.__interp_del(q, [self.cdf[i], self.cdf[i+1]])

        return i_true

    def __cumsum(self):
        count, _ = np.histogram(self.pd_array, bins=255) 
        self.pdf = count / sum(count) 
        self.cdf = np.cumsum(self.pdf) 