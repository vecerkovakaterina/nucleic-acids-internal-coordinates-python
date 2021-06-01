from coordinates.coordinatesCurves import CoordinatesCurves
import numpy as np


class CoordinatesCgDNA(CoordinatesCurves):
    """A class to represent internal coordinates according to the cgDNA definition."""

    def get_theta(self, rotation_matrices):
        """Takes the array of rotation matrices.
        Returns an array of rotation angles."""
        theta_array = [np.arccos((np.trace(x) - 1) / 2) for x in rotation_matrices]
        theta_array = [np.degrees(2 * np.tan(t / 2)) for t in theta_array]
        return theta_array
