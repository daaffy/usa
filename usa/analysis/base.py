import SimpleITK as sitk
import numpy as np
from scipy import ndimage

from usa.analysis.ge import GETagAnalyser
from usa.analysis.ge.distance import get_angle
from usa.utils.loader import load_nifti_from_path

"""
    base.py
Base objects for ultrasound image analysis.

=========================================================================================================================
--- US Bundles
US:
Basic ultrasound image bundle containing, e.g., power Doppler array, B-mode array, segmentation array...

GE_US:
GE specific ultrasound image bundle used for traditional FMBV analysis. Kretz (.vol) file metadata is included.

NOTE:
Eventually we will also write a bundle for the UNSW-IDE DICOM files...

=========================================================================================================================
--- Methods
US_Method:
Basic ultrasound image bundle processing object. Includes verbose flag functionality.

=========================================================================================================================
--- Other
resample_GE:
Resamples a GE_US bundle to a different image resolution.

vrbs:
Provides verbose functionality.

=========================================================================================================================
"""

class US():

    def __init__(self, vol_id = 1):
        """
            Bundle of ultrasound data that will be further processed, e.g., calculating FMBV and other image statistics.
        """
        self.bm_array, self.bm_supplied = None, False
        self.pd_array, self.pd_supplied = None, False
        self.seg_array, self.seg_supplied = None, False
        self.distance = None

        assert isinstance(vol_id, int)
        self.vol_id = vol_id

class GE_US(US):

    def __init__(self, 
                 minimise_memory = False,
                 verbose = False
                 ):
        """
            Ultrasound bundle specific to the GE machine with .vol output.

            Traditionally, we convert the raw .vol into B-Mode and Power Doppler .nii.gz files. Distance information is obtained from the .vol metadata.
            We don't usually care about B-Mode information for FMBV calculation, however, functionality can be added if required later on.

            Functionality:
            - Load power Doppler, load kretz (.vol), load segmentation.
            - Calculate volume (mm^3) of the segmentation.
            - Calculate transducer distance information from the kretz metadata. 
        """
        US.__init__(self)
        self.kretz_supplied = False
        self.seg_volume = None      # mm^3
        self.pd_img = None
        self.resolution = None

        self.minimise_memory = minimise_memory
        self.verbose = verbose

        vrbs("Memory minimise: "+str(self.minimise_memory), self.verbose)

        # NOTE Could add some argument handling to automate the loading process...
    
    # Load Functions
    def load_kretz(self,
                   input
                   ):
        
        # Only metadata required from .vol file
        if isinstance(input, str):
            self.kretz_metadata = GETagAnalyser.GETagAnalyser(input)

            self.bmode_spacing = np.frombuffer(self.kretz_metadata.getDataFromTags(
                int('0xc100', 16), int('0x201', 16)), np.float64)[0] * 1000.0
            self.volume_offset = np.frombuffer(self.kretz_metadata.getDataFromTags(
                int('0xc200', 16), int('0x202', 16)), np.float64)[0]
            # other metadata...?

            self.sweep_radius = self.bmode_spacing * self.volume_offset
        else:
            raise Exception("Invalid Kretz input.")
        
        self.kretz_supplied = True

    def load_pd(self,
                input,
                ):
        self.pd_array, self.pd_img = self.__load_img(input)
        self.pd_supplied = True

        if not self.pd_img == None:
            self.resolution = self.pd_img.GetSpacing()[0] # Get resolution for, e.g., volume measurements. Assume isotropic resolution.

    def load_seg(self,
                 input,
                 ):
        self.seg_array, _ = self.__load_img(input)
        self.seg_supplied = True

        self.calculate_volume()

    def load_default_seg(self,
                        mode = "sweep", 
                        floor = 0
                                ):
        self.seg_array = self._get_default_segmentation(mode = mode, floor = floor)
        self.seg_supplied = True

        self.calculate_volume()

    # Other
    def calculate_volume(self):
        assert self.seg_supplied

        seg_voxels = np.sum(self.seg_array == self.vol_id)

        if self.resolution == None:
            self.seg_volume = 0
        else:
            self.seg_volume = seg_voxels * self.resolution**3

    def _get_default_segmentation(self, 
                                  mode = "sweep", 
                                  floor = 0
                                  ):
        """
            Return a default segmentation

            Parameters:
            mode            :
                = "sweep"   : segmentation defined as the sweep of the power Doppler array (non-zero values *usually* define the curvilinear sweep of the ultrasound).
                = "blank"   : segmentation defined as the whole domain.
            floor           : set floor that defines how sweep is defined.
        """
        if mode == "sweep":
            if (self.pd_supplied):
                return np.multiply(np.ones(self.pd_array.shape), (self.pd_array > floor))
            else:
                raise Exception("Power Doppler not supplied to base default segmentation from.")
        elif mode == "blank":
                return np.ones(self.pd_array.shape)

    def __load_img(self, 
            input
            ):
        """
            Load images specific to the current GE processing workflow.
        """

        if isinstance(input, np.ndarray):
            array = input
            img = None
        elif isinstance(input, str):
            array, img = load_nifti_from_path(input)
        else:
            raise TypeError("Input must be a string (path) or a numpy array.")

        return array, img

    def set_distance(self):
        assert not self.pd_img == None and self.kretz_supplied
        self.distance = self.get_distance()

    def get_distance(self):
        assert not self.pd_img == None

        vrbs("Calculating transducer distance matrix...", self.verbose)

        dimensions = self.pd_img.GetSize() # do we have to use pd_img here? can we get this info. from pd_array?
        corner1 = self.pd_img.TransformIndexToPhysicalPoint((0,0,0))
        corner2 = self.pd_img.TransformIndexToPhysicalPoint((dimensions[0]-1,dimensions[1]-1,dimensions[2]-1))

        dimensions_pd = np.shape(self.pd_array)

        x = np.linspace(corner1[0], corner2[0], dimensions_pd[2])
        y = np.linspace(corner1[1], corner2[1], dimensions_pd[1])
        z = np.linspace(corner1[2], corner2[2], dimensions_pd[0])

        xx, yy, zz = np.meshgrid(x, y, z)

        _, _, rDist = get_angle(xx, yy, zz, self.sweep_radius, minimise_memory=self.minimise_memory, verbose=self.verbose)
        distance = np.transpose(rDist, (2, 0, 1))

        return distance

class US_Method:

    def __init__(self, 
                 verbose = False, 
                 **kwargs
                 ):
        """
            Basic US_Method object.
        """
        self.verbose = verbose

    def _vrbs(
            self, 
            msg
            ):
        """
            Handles verbose flag: if self.verbose is True, print message.
        """
        vrbs(msg, self.verbose)

def resample_GE(
    us : GE_US,
    zoom
): 
    """
        Return a resampled US bundle.
    """

    def resample(array):
        if not zoom == 1:
            array = ndimage.zoom(array, (zoom, zoom, zoom))
        return array

    us_resample = GE_US(minimise_memory=us.minimise_memory)
    us_resample.load_pd(resample(us.pd_array))
    us_resample.load_seg(np.round(resample(us.seg_array)).astype(int))
    
    us_resample.resolution = us.resolution/zoom
    us_resample.calculate_volume()
    # us_resample.seg_volume = us.seg_volume
    
    us_resample.pd_img = us.pd_img
    us_resample.kretz_supplied = us.kretz_supplied
    us_resample.kretz_metadata = us.kretz_metadata
    us_resample.bmode_spacing = us.bmode_spacing
    us_resample.sweep_radius = us.sweep_radius

    us_resample.set_distance()

    return us_resample

def vrbs(
        msg, 
        verbose: bool
        ):
    """
        I've found this function to be useful for debugging the development of new methods!
    """
    if verbose:
        try:
            print('[verbose] ' + msg)
        except(TypeError):
            try:
                print('[verbose]... ')
                print(msg)
            except:
                print('[verbose] Trouble printing verbose message.')
