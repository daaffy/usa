from SimpleITK import ReadImage, GetArrayFromImage

def load_nifti_from_path(
        path: str
        ):
    """
        Load a .nii.gz file as a numpy array and a sitk image.

        Parameters:
        path (str)      : Path to .nii.gz file.

        Returns:
        np.array        : Pixel data
        sitk.Image      : Image object
    """

    img = ReadImage(path)
    array = GetArrayFromImage(img)

    return array, img
