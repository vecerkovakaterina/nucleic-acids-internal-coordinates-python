import pandas as pd
import numpy as np
from functools import reduce


class FrameFitter:
    """A class to represent the reference frame fitting."""

    def __init__(self):
        """Construct all the necessary attributes for the frame fitting procedure."""
        self.frames_1 = []
        self.frames_2 = []
        self.origins_1 = []
        self.origins_2 = []
        self.bases_1 = []
        self.bases_2 = []
        self.experimental_1 = []
        self.experimental_2 = []
        self.number_of_bases = 0
        self.strand_len = 0

    # dict of relevant atoms in nucleotide bases
    ring_atoms = {'G': ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'],
                  'A': ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'],
                  'C': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
                  'U': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
                  'T': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']}

    # dict of standard frames coordinates accessible by nucleotide letter and atom number
    standard_frames = {
        'G': {'N9': [-1.289, 4.551, 0.000],
              'C8': [0.023, 4.962, 0.000],
              'N7': [0.870, 3.969, 0.000],
              'C5': [0.071, 2.833, 0.000],
              'C6': [0.424, 1.460, 0.000],
              'N1': [-0.700, 0.641, 0.000],
              'C2': [-1.999, 1.087, 0.000],
              'N3': [-2.342, 2.364, 0.001],
              'C4': [-1.265, 3.177, 0.000]},
        'C': {'N1': [-1.285, 4.542, 0.000],
              'C2': [-1.472, 3.158, 0.000],
              'N3': [-0.391, 2.344, 0.000],
              'C4': [0.837, 2.868, 0.000],
              'C5': [1.056, 4.275, 0.000],
              'C6': [-0.023, 5.068, 0.000]},
        'A': {'N9': [-1.291, 4.498, 0.000],
              'C8': [0.024, 4.897, 0.000],
              'N7': [0.877, 3.902, 0.000],
              'C5': [0.071, 2.771, 0.000],
              'C6': [0.369, 1.398, 0.000],
              'N1': [-0.668, 0.532, 0.000],
              'C2': [-1.912, 1.023, 0.000],
              'N3': [-2.320, 2.290, 0.000],
              'C4': [-1.267, 3.124, 0.000]},
        'T': {'N1': [-1.284, 4.500, 0.000],
              'C2': [-1.462, 3.135, 0.000],
              'N3': [-0.298, 2.407, 0.000],
              'C4': [0.994, 2.897, 0.000],
              'C5': [1.106, 4.338, 0.000],
              'C6': [-0.024, 5.057, 0.000]},
        'U': {'N1': [-1.284, 4.500, 0.000],
              'C2': [-1.462, 3.131, 0.000],
              'N3': [-0.302, 2.397, 0.000],
              'C4': [0.989, 2.884, 0.000],
              'C5': [1.089, 4.311, 0.000],
              'C6': [-0.024, 5.053, 0.000]}
    }

    def read_pdb(self, pdb_filename):
        """Takes pdb file name and reads the experimental frames."""
        pdb = pd.read_csv(pdb_filename, delim_whitespace=True, header=None)

        # find bases
        pdb = pdb.loc[pd.isnull(pdb[3]) != True]
        pdb = pdb.loc[pdb[3].str.startswith('D')]
        pdb['base'] = pdb[3].str[1]

        # length of one strand
        self.number_of_bases = int(pdb.iloc[-1, 4])
        self.strand_len = int(self.number_of_bases / 2)

        # nucleotide sequence in both strands (complementary 5' --> 3' 5' --> 3)
        bases_both_strands = [pdb.loc[pdb[4] == i, 'base'].reset_index(drop=True)[0] for i in
                              range(1, self.number_of_bases + 1)]
        self.bases_1 = bases_both_strands[:self.strand_len]
        self.bases_2 = bases_both_strands[self.strand_len:]

        experimental_frames = []
        for i in range(1, self.number_of_bases + 1):
            ex_frame = []
            for atom in self.ring_atoms[pdb.loc[pdb[4] == i, 'base'].values[0]]:
                ex_frame.append(pdb.loc[(pdb[4] == i) & (pdb[2] == atom), [2, 5, 6, 7]].values[0])
            experimental_frames.append(pd.DataFrame.from_records(ex_frame).set_index([0]))

        # separate frames into strands
        self.experimental_1 = experimental_frames[:self.strand_len]
        self.experimental_2 = experimental_frames[self.strand_len:]

        # reverse strand 2
        self.bases_2 = self.bases_2[::-1]
        self.experimental_2 = self.experimental_2[::-1]

    @staticmethod
    def get_covariance_matrix(standard, experimental, no_atoms):
        """Takes a standard frame matrix, experimental frame matrix and number of relevant atoms in a base.
        Returns the covariance matrix."""
        ones = np.ones(shape=(no_atoms, 1))
        cm = (1 / (no_atoms - 1)) * (reduce(np.matmul, [standard.T, experimental])
                                     - (1 / no_atoms) * reduce(np.matmul, [standard.T, ones, ones.T, experimental]))
        return cm

    @staticmethod
    def get_symmetric_matrix(cm):
        """Takes the covariance matrix.
        Returns the symmetric matrix."""
        sm = np.array([
            [cm[0, 0] + cm[1, 1] + cm[2, 2], cm[1, 2] - cm[2, 1], cm[2, 0] - cm[0, 2], cm[0, 1] - cm[1, 0]],
            [cm[1, 2] - cm[2, 1], cm[0, 0] - cm[1, 1] - cm[2, 2], cm[0, 1] + cm[1, 0], cm[2, 0] + cm[0, 2]],
            [cm[2, 0] - cm[0, 2], cm[0, 1] + cm[1, 0], -cm[0, 0] + cm[1, 1] - cm[2, 2], cm[1, 2] + cm[2, 1]],
            [cm[0, 1] - cm[1, 0], cm[2, 0] + cm[0, 2], cm[1, 2] + cm[2, 1], -cm[0, 0] - cm[1, 1] + cm[2, 2]]
        ])
        return sm

    @staticmethod
    def get_rotation_matrix(sm):
        """Takes the symmetric matrix.
        Returns the rotation matrix."""
        eigval = np.max(np.linalg.eigvals(sm))
        index = np.where(np.linalg.eigvals(sm) == eigval)
        eigvals, eigvecs = np.linalg.eig(sm)
        eigenvec = -eigvecs.T[index][0]

        rm = np.array([
            [eigenvec[0] ** 2 + eigenvec[1] ** 2 - eigenvec[2] ** 2 - eigenvec[3] ** 2,
             2 * (eigenvec[1] * eigenvec[2] - eigenvec[0] * eigenvec[3]),
             2 * (eigenvec[1] * eigenvec[3] + eigenvec[0] * eigenvec[2])],
            [2 * (eigenvec[2] * eigenvec[1] + eigenvec[0] * eigenvec[3]),
             eigenvec[0] ** 2 - eigenvec[1] ** 2 + eigenvec[2] ** 2 - eigenvec[3] ** 2,
             2 * (eigenvec[2] * eigenvec[3] - eigenvec[0] * eigenvec[1])],
            [2 * (eigenvec[3] * eigenvec[1] - eigenvec[0] * eigenvec[2]),
             2 * (eigenvec[3] * eigenvec[2] + eigenvec[0] * eigenvec[1]),
             eigenvec[0] ** 2 - eigenvec[1] ** 2 - eigenvec[2] ** 2 + eigenvec[3] ** 2]])
        return rm

    @staticmethod
    def get_origin(rm, standard, experiment):
        """Takes the rotation matrix, standard frame matrix and experimental frame matrix.
        Returns the reference point of the fitted frame."""
        e_ave = experiment.mean(axis=0)
        s_ave = standard.mean(axis=0)
        o = e_ave - s_ave.dot(rm.T)
        return o

    @staticmethod
    def get_fitted_coordinates(rm, standard, o):
        """Takes the rotation matrix, the standard frame matrix and the reference point.
        Returns the fitted coordinates."""
        fitted_coordinates = np.array([arr + o for arr in standard.dot(rm.T)])
        return fitted_coordinates

    @staticmethod
    def get_matrices_and_no_atoms(exp_frames, std_frames, bases, index):
        """Takes the experimental frame array, the standard frames array, the base letter and index.
        Returns the experimental frame matrix, the corresponding standard frame and number of atoms in frame."""
        exp_matrix = exp_frames[index].to_numpy()
        std_matrix = np.array(list(std_frames[bases[index]].values()))
        n = len(exp_frames[index].index)
        return exp_matrix, std_matrix, n

    @staticmethod
    def get_frame_and_origin(standard_matrix, experimental_matrix, n):
        """Takes the standard frame matrix, the experimental frame matrix and the number of atoms.
        Returns the rotation matrix and its reference point."""
        covariance_matrix = FrameFitter.get_covariance_matrix(standard_matrix, experimental_matrix, n)
        symmetric_matrix = FrameFitter.get_symmetric_matrix(covariance_matrix)
        rotation_matrix = FrameFitter.get_rotation_matrix(symmetric_matrix)
        origin = FrameFitter.get_origin(rotation_matrix, standard_matrix, experimental_matrix)
        return rotation_matrix, origin

    def get_frames_and_origins(self):
        """Computes the fitted frames."""
        for i in range(0, self.strand_len):
            exp_matrix_1, std_matrix_1, n_1 = FrameFitter.get_matrices_and_no_atoms(self.experimental_1,
                                                                                    self.standard_frames,
                                                                                    self.bases_1, i)
            exp_matrix_2, std_matrix_2, n_2 = FrameFitter.get_matrices_and_no_atoms(self.experimental_2,
                                                                                    self.standard_frames,
                                                                                    self.bases_2, i)

            frame_1, origin_1 = FrameFitter.get_frame_and_origin(std_matrix_1, exp_matrix_1, n_1)
            self.frames_1.append(frame_1)
            self.origins_1.append(origin_1)

            frame_2, origin_2 = FrameFitter.get_frame_and_origin(std_matrix_2, exp_matrix_2, n_2)
            self.frames_2.append(frame_2)
            self.origins_2.append(origin_2)

    def rotate_strand_2(self):
        """Rotates frame in strand 2 180 degrees around x-axis."""
        x_rot = np.array([
            [1, 0, 0],
            [0, -1, 0],
            [0, 0, -1]
        ])
        self.frames_2 = [np.matmul(frame, x_rot) for frame in self.frames_2]

    def run(self, pdb_filename):
        """Computes the fitted frames in both strands."""
        self.read_pdb(pdb_filename)
        self.get_frames_and_origins()
        self.rotate_strand_2()
