"""
    Fractional Moving Blood Volume (FMBV) processing.
"""

import copy

from usa.analysis.base import *

# Cumulative Distribution Algorithms
from usa.analysis.cumulative import v0

# --------------------------------------------------------------------------------------------------
# FMBVs
class fmbv_v0(US_Method):

    def __init__(self, 
                 us_bundle : US,
                 **kwargs):
        """
            Original FMBV algorithm using the two-tangent method for standisation.

            Rubin et al. Normalizing fractional moving blood volume estimates with power Doppler US: defining a stable intravascular point with the cumulative power distribution function.
        """
        # Default Parameters
        # self.mode = 0                   # Select processing mode. 
        self.verbose = False            # print what's going on...
        # self.pre_scale = 1
        # self.scale = 1.00               # Scale final standardised image to this value.
        self.pd_low_threshold = 2.0     # Voxels values below this value are discarded in FMBV processing
        # self.combine = True             # Ignore segmentation partitions and treat all non-zero indices as a whole segmentation
        self.vol_id = 1                 # 
        self.depth_correction_layer_thickness = 2           #
        self.max_pixel_value = 1.0      # For constructing pd_array_global_standardisation
        # self.dp_img, self.seg_img = None, None
        # self.vis_information = None
        self.clip = True                # clip standardised pd_array
        # self.std_regress = False
        # self.original_seg_voxel_nums = 0 # this needs better treatment
        self.minimise_memory = False    # process differently depending on memory constraints

        self.__load_kwargs(**kwargs)

        # Initialise
        US_Method.__init__(self, verbose=self.verbose)
        self._vrbs("--- New FMBV(version=0) instance. ---")

        self.us_bundle = us_bundle
        self.__check_us_bundle()

        self.seg_volume = self.us_bundle.seg_volume
    
    def pre_process(self):
        # self.scaled_spacing = self.us_bundle.pd_img.GetSpacing()[0] # assume isotropic scaling
        # self.seg_volume = self.original_seg_voxel_nums * (self.scaled_spacing/self.zoom)**3 # mm^3
        0

    def global_method(self):
        self._vrbs("* I am in global method...")

        # Flatten and mask the data before processing.
        pd_flat = self.pd_array.flatten().astype(float)
        seg_flat = self.seg_array.flatten().astype(int)

        mask = (seg_flat == self.vol_id)

        pd_trim = pd_flat[mask]
        self.mpi = np.mean(pd_trim)

        # Calculate the standardisation values.
        self.global_standardisation_value_1, self.global_standardisation_value_2, self.global_vis_1, self.global_vis_2 = self.std_method(pd_trim[pd_trim > self.pd_low_threshold])

        # Calculate the FMBV values by averaging the normalised data.
        self._vrbs("Calculating FMBVs...")
        self.global_fmbv_1 = v0.normalise(copy.deepcopy(pd_trim), self.global_standardisation_value_1)
        self.global_fmbv_2 = v0.normalise(copy.deepcopy(pd_trim), self.global_standardisation_value_2)

        self._vrbs("Global FMBV: " + str(self.global_fmbv_2))
        
        # Visualise the standardised power Doppler image.
        if not self.minimise_memory: # e.g., if mm=True (we want to minimize memory) we ignore the storage of large arrays.
            self.pd_std_global_array = np.multiply(self.pd_array, self.seg_array/self.global_standardisation_value_2) # Note that we multiply by the segmentation mask here.
            if self.clip:
                self.pd_std_global_array[self.pd_std_global_array >= 1] = 1.0

            self.pd_array_global_standardisation = self.pd_std_global_array*self.max_pixel_value

    def depth_correction_method(self):
        self._vrbs("* I am in depth-correction method...")

        assert self.us_bundle.kretz_supplied

        distance = self.us_bundle.distance
        img_distance = (self.seg_array == self.vol_id) * distance
        if not self.minimise_memory:
            self.layers = img_distance
            self.pd_array_depth_corrected_standardisation = np.zeros(np.shape(self.pd_array))

        pd_flat = self.pd_array.flatten().astype(int)
        seg_flat = img_distance.flatten().astype(int)
        ma = seg_flat[seg_flat != 0]

        self.start_depth = np.min(ma)
        self.end_depth = np.max(ma)
        self.scaled_spacing = self.us_bundle.resolution

        c = self.depth_correction_layer_thickness

        self._vrbs("-> Entering onion...")
        fmbv_vals, vxl_nums = [], []
        self.depth_correction_vis_1, self.depth_correction_vis_2 = [], []
        self.depth_correction_standard_values_1, self.depth_correction_standard_values_2 = [], []

        success_num = 0
        layer_range = np.arange(ma.min(), ma.max(), c)
        self.depths = layer_range + .5*c
        self.numlayers = len(layer_range)
        for i in range(self.numlayers):
            d = layer_range[i]
            self._vrbs('\n Layer '+ str(i+1) + ' out of '+ str(self.numlayers))

            mask = (seg_flat >= d) & (seg_flat <= d + c)
            array_mask = (img_distance >= d) & (img_distance <= d + c)

            pd_fmbv = pd_flat[mask]
            vxl_num = sum(mask)

            std_value_1, std_value_2, vis1, vis2 = self.std_method(pd_fmbv[pd_fmbv > self.pd_low_threshold])

            self.depth_correction_vis_1.append(vis1)
            self.depth_correction_vis_2.append(vis2)
            self.depth_correction_standard_values_1.append(std_value_1)
            self.depth_correction_standard_values_2.append(std_value_2)

            try: 
                fmbv_value = np.nan
                if std_value_2 is not None:
                    fmbv_value = v0.normalise(copy.deepcopy(pd_fmbv), std_value_2)
                    if not self.minimise_memory:
                        self.pd_array_depth_corrected_standardisation = self.pd_array_depth_corrected_standardisation + np.multiply(self.pd_array, array_mask/std_value_2)
                    success_num = success_num + 1
            except Exception as e:
                self._vrbs(str(i))
                self._vrbs(e)

                fmbv_value = np.nan

            fmbv_vals.append(fmbv_value)
            vxl_nums.append(vxl_num)
        
        self._vrbs("<- Exiting onion...")
        self._vrbs("Successful layers: " + str(int(success_num / self.numlayers * 100)) + "%")

        self.vxl_nums = vxl_nums
        self.depth_corrected_fmbv = np.nanmean(fmbv_vals)

        self._vrbs("Depth-Corrected FMBV: " + str(self.depth_corrected_fmbv))

        if not self.minimise_memory:
            if self.clip:
                self.pd_array_depth_corrected_standardisation[self.pd_array_depth_corrected_standardisation >= 1] = 1.
            self.pd_array_depth_corrected_standardisation = self.pd_array_depth_corrected_standardisation*self.max_pixel_value

    def std_method(self, pd_data):
        return v0.two_tangent_standardisation(
            pd_data,
            self.pd_low_threshold,
            visualise = False, 
            verbose = self.verbose
            )

    def __check_us_bundle(self):
        if self.us_bundle.pd_supplied:
            self.pd_array = self.us_bundle.pd_array
        else:
            raise Exception("No power Doppler supplied in ultrasound bundle.")
        
        if self.us_bundle.seg_supplied:
            self.seg_array = self.us_bundle.seg_array
        else:
            raise Exception("No segmentation supplied in ultrasound bundle.")
        
        if (not np.array_equiv(self.pd_array.shape, self.seg_array.shape)):
            raise Exception("Power doppler array size", self.pd_array.shape, " and segmentation array ", self.seg_array.shape, " do not have the same dimensions.")
        
    def __load_kwargs(self, **kwargs):
        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]
        if "mode" in kwargs:
            self.mode = kwargs["mode"]
        if "scale" in kwargs:
            self.scale = kwargs["scale"]
        if "vol_id" in kwargs:
            # We assume here that the user intends on selecting a specific index to perform FMBV on
            self.combine = False
            self.vol_id = kwargs["vol_id"]
        if "max_pixel_value" in kwargs:
            self.max_pixel_value = kwargs["max_pixel_value"]
        if "pd_low_threshold" in kwargs:
            self.pd_low_threshold = kwargs["pd_low_threshold"]
        if "pre_scale" in kwargs:
            self.pre_scale = kwargs["pre_scale"]
        if "clip" in kwargs:
            self.clip = kwargs["clip"]
        if "std_regress" in kwargs:
            self.std_regress = kwargs["std_regress"]
        if "minimise_memory" in kwargs:
            self.mm = kwargs["minimise_memory"]
        if "zoom" in kwargs:
            self.zoom = kwargs["zoom"]
