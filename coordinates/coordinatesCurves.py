from coordinates.coordinates import Coordinates
import numpy as np


class CoordinatesCurves(Coordinates):
    """A class to represent internal coordinates according to the Curves+ definition."""

    @staticmethod
    def get_rotation_matrices(frames_1, frames_2):
        """Takes two arrays of frames.
        Returns an array of rotation matrices."""
        rm = [np.matmul(b1, b2.T) for b1, b2 in zip(frames_1, frames_2)]
        return rm

    def get_theta(self, rotation_matrices):
        """Takes an array of rotation matrices.
        Returns an array of rotation angles."""
        theta_array = [
            np.degrees(np.arccos((np.trace(x) - 1) / 2)) for x in rotation_matrices
        ]
        return theta_array

    @staticmethod
    def get_eigenvectors_nonsymmetric(rotation_matrices):
        """Takes an array of rotation matrices.
        Returns an array of eigenvectors."""
        w = []
        for x in rotation_matrices:
            w.append([x[1, 2] - x[2, 1], x[2, 0] - x[0, 2], x[0, 1] - x[1, 0]])
        return w

    @staticmethod
    def get_rotation_matrices_rodrigues(rotation_angle, u):
        """Takes an array of rotation angles and unit rotation vectors.
        Returns an array of rotation matrices."""
        q = []

        id_matrix = np.identity(3)
        for angle, vector in zip(rotation_angle, u):
            u_a_x = np.array(
                [
                    [0, -vector[2], vector[1]],
                    [vector[2], 0, -vector[0]],
                    [-vector[1], vector[0], 0],
                ]
            )

            q.append(
                id_matrix * (np.cos(angle))
                + (1 - np.cos(angle)) * np.outer(vector, vector)
                + (np.sin(angle)) * u_a_x
            )
        return q

    @staticmethod
    def gram_schmidt_columns(matrix):
        """Takes a matrix, returns the matrix orthonormalized."""
        n = matrix.shape[1]
        for j in range(n):
            for k in range(j):
                matrix[:, j] -= np.dot(matrix[:, k], matrix[:, j]) * matrix[:, k]
            matrix[:, j] = matrix[:, j] / np.linalg.norm(matrix[:, j])
        return matrix

    @staticmethod
    def get_translational_coords(lambda_vec, middle_frame):
        """Takes an array of displacement vectors between frame origins and an array of middle frames.
        Returns three arrays of translational coordinates."""
        coor1, coor2, coor3 = [], [], []

        for vec, matrix in zip(lambda_vec, middle_frame):
            x, y, z = np.dot(vec, matrix)
            coor1.append(x)
            coor2.append(y)
            coor3.append(z)

        return coor1, coor2, coor3

    @staticmethod
    def get_rotational_coords(theta, unit_vector, middle_frame):
        """Takes arrays of rotation angle, rotation vectors and middle frames.
        Returns three arrays of rotational coordinates."""
        coor1, coor2, coor3 = [], [], []

        for angle_t, u_vector, mid_frame in zip(theta, unit_vector, middle_frame):
            x, y, z = np.dot([angle_t * u for u in u_vector], mid_frame)
            coor1.append(x)
            coor2.append(y)
            coor3.append(z)

        return coor1, coor2, coor3

    @staticmethod
    def get_middle_frames(frames_1, frames_2, origins_1, origins_2):
        """Takes arrays of frames in two strands and their reference points.
        Returns an array of middle frames and their reference points."""
        middle_frames = [
            (frame1 + frame2) / 2 for frame1, frame2 in zip(frames_1, frames_2)
        ]
        origins_middle_frames = [
            (origin1 + origin2) / 2 for origin1, origin2 in zip(origins_1, origins_2)
        ]
        middle_frames = [mf / np.linalg.norm(mf) for mf in middle_frames]
        middle_frames = [CoordinatesCurves.gram_schmidt_columns(mf) for mf in middle_frames]

        return middle_frames, origins_middle_frames

    def get_intra_coords(self):
        """Computes the base frame coordinates."""
        intra_rotation_matrices = CoordinatesCurves.get_rotation_matrices(
            self.frames_1, self.frames_2
        )
        theta_a = self.get_theta(intra_rotation_matrices)
        lambda_a = [(x2 - x1) for x1, x2 in zip(self.origins_1, self.origins_2)]
        w = CoordinatesCurves.get_eigenvectors_nonsymmetric(intra_rotation_matrices)
        unit_rotation_vector_a = [[v / np.linalg.norm(vec) for v in vec] for vec in w]
        self.shear, self.stagger, self.stretch = CoordinatesCurves.get_translational_coords(
            lambda_a, self.intra_middle_frames
        )
        self.propeller, self.buckle, self.opening = CoordinatesCurves.get_rotational_coords(
            theta_a, unit_rotation_vector_a, self.intra_middle_frames
        )

    def get_inter_coords(self):
        """Computes the base pair frames coordinates."""
        inter_rotation_matrices = CoordinatesCurves.get_rotation_matrices(
            self.intra_middle_frames, self.intra_middle_frames[1:]
        )
        theta_e = self.get_theta(inter_rotation_matrices)
        w_e = CoordinatesCurves.get_eigenvectors_nonsymmetric(inter_rotation_matrices)
        unit_rotation_vector_e = [[v / np.linalg.norm(vec) for v in vec] for vec in w_e]
        lambda_e = [
            (x2 - x1)
            for x1, x2 in zip(
                self.intra_middle_frames_origins, self.intra_middle_frames_origins[1:]
            )
        ]
        self.shift, self.slide, self.rise = CoordinatesCurves.get_translational_coords(
            lambda_e, self.inter_middle_frames
        )
        self.roll, self.tilt, self.twist = CoordinatesCurves.get_rotational_coords(
            theta_e, unit_rotation_vector_e, self.inter_middle_frames
        )

    def run(self):
        """Compute the internal coordinates."""
        # intra
        (
            self.intra_middle_frames,
            self.intra_middle_frames_origins,
        ) = CoordinatesCurves.get_middle_frames(
            self.frames_1, self.frames_2, self.origins_1, self.origins_2
        )

        self.get_intra_coords()

        # inter
        (
            self.inter_middle_frames,
            self.inter_middle_frames_origins,
        ) = CoordinatesCurves.get_middle_frames(
            self.intra_middle_frames,
            self.intra_middle_frames[1:],
            self.intra_middle_frames_origins,
            self.intra_middle_frames_origins[1:],
        )

        self.get_inter_coords()
