from abc import ABC, abstractmethod
import numpy as np


def round_array_to_2_decimals(array):
    return np.round(array, 2)


class Coordinates(ABC):
    """A class to represent internal coordinates."""

    def __init__(self, frames_1, frames_2, origins_1, origins_2, strand_len):
        """Constructs all the necessary attributes for the coordinates object."""
        self.frames_1 = frames_2
        self.frames_2 = frames_1
        self.origins_1 = origins_2
        self.origins_2 = origins_1
        self.intra_len = strand_len
        self.inter_len = self.intra_len - 1
        self.intra_middle_frames = []
        self.intra_middle_frames_origins = []
        self.inter_middle_frames = []
        self.inter_middle_frames_origins = []
        self.shear = None
        self.stretch = None
        self.stagger = None
        self.buckle = None
        self.propeller = None
        self.opening = None
        self.shift = None
        self.slide = None
        self.rise = None
        self.roll = None
        self.tilt = None
        self.twist = None

    @abstractmethod
    def run(self):
        pass

    def round_coords(self):
        self.shear = round_array_to_2_decimals(self.shear)
        self.stretch = round_array_to_2_decimals(self.stretch)
        self.stagger = round_array_to_2_decimals(self.stagger)
        self.buckle = round_array_to_2_decimals(self.buckle)
        self.propeller = round_array_to_2_decimals(self.propeller)
        self.opening = round_array_to_2_decimals(self.opening)
        self.shift = round_array_to_2_decimals(self.shift)
        self.slide = round_array_to_2_decimals(self.slide)
        self.rise = round_array_to_2_decimals(self.rise)
        self.roll = round_array_to_2_decimals(self.roll)
        self.tilt = round_array_to_2_decimals(self.tilt)
        self.twist = round_array_to_2_decimals(self.twist)

    @staticmethod
    def append_to_file(filename, array):
        try:
            file = open(filename, 'a+')
            line_to_append = '\t'.join(map(str, array))
            file.write(line_to_append + '\n')
            file.close()
        except IOError:
            print("Output file path does not exist!")

    def write_to_bp_output_files(self, output_filepath):
        """
        Write intra-bp coordinates to output files.
        :param output_filepath: output file dir and root name
        """
        for i in range(0, self.intra_len):
            output_filepath_name = output_filepath + "_bp_" + str(i + 1) + ".out"
            bp_array = [self.shear[i], self.stretch[i], self.stagger[i], self.buckle[i], self.propeller[i],
                        self.opening[i]]
            Coordinates.append_to_file(output_filepath_name, bp_array)

    def write_to_bp_step_output_files(self, output_filepath):
        """
        Write inter-bp coordinates to output files.
        :param output_filepath: output file dir and root name
        """
        for i in range(0, self.inter_len):
            output_filepath_name = output_filepath + "_step_" + str(i + 1) + ".out"
            bp_step_array = [self.shift[i], self.slide[i], self.rise[i], self.roll[i], self.tilt[i], self.twist[i]]
            Coordinates.append_to_file(output_filepath_name, bp_step_array)

    def write_to_test_accuracy(self, output_filepath, pdb_number):
        """
        Write coordinates to output file for each PDB file separately.
        Comparable to output from 3DNA bp_step.par.
        For testing purposes.
        :param output_filepath: output file fir and root name
        :param pdb_number: Number of PDB file being processed
        """
        output_filepath_name = output_filepath + "_test_" + str(pdb_number) + ".out"
        to_append = [self.shear[0], self.stretch[0], self.stagger[0], self.buckle[0], self.propeller[0],
                    self.opening[0], 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        Coordinates.append_to_file(output_filepath_name, to_append)
        for i in range(0, self.inter_len):
            bp_step = [self.shear[i+1], self.stretch[i+1], self.stagger[i+1], self.buckle[i+1], self.propeller[i+1],
                        self.opening[i+1], self.shift[i], self.slide[i], self.rise[i], self.roll[i], self.tilt[i], self.twist[i]]
            Coordinates.append_to_file(output_filepath_name, bp_step)
