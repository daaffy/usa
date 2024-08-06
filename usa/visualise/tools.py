import matplotlib.pyplot as plt
from usa.analysis.base import GE_US
from usa.analysis.fmbv import fmbv_v0
import numpy as np

class visualise_fmbv_v0():

    def __init__(self, us : GE_US, f : fmbv_v0):
        self.us = us
        self.f = f
        self.pd_array = us.pd_array

        self.j_ = int(.5*self.pd_array.shape[1])

        # dist = us.distance

        self.std1, self.std2 = f.global_standardisation_value_1, f.global_standardisation_value_2

        self.visualise_pd_standardisation()

    def visualise_pd_standardisation(self):
        fig, ax = plt.subplots(3)
        fig.set_figheight(8)
        fig.set_figwidth(5)

        # 1. ---
        ax[0].imshow(self.pd_array[:,self.j_,:],cmap='bone')

        # 2. ---
        ax[1].imshow(self.pd_array[:,self.j_,:],cmap='bone')

        # standardisation 1
        masked1 = np.ma.masked_where(self.pd_array < self.std1,np.ones(self.pd_array.shape))
        ax[1].imshow(masked1[:,int(.5*self.pd_array.shape[1]),:],alpha=0.6,cmap='Blues',vmin=0, vmax=1)

        # standardisation 2
        masked2 = np.ma.masked_where(self.pd_array < self.std2, np.ones(self.pd_array.shape))
        ax[1].imshow(masked2[:,int(.5*self.pd_array.shape[1]),:],alpha=0.6,cmap = 'Reds',vmin=0, vmax=1)

        # 3. ---
        ax[2].imshow(self.f.pd_array_global_standardisation[:,self.j_,:],cmap='bone')


"""
    visualise_fmbv_v2
"""