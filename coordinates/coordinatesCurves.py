from coordinates.coordinates import Coordinates
import numpy as np


class CoordinatesCurves(Coordinates):
    def get_rotation_matrices(self, frames_1, frames_2):
        rm = [np.matmul(b1, b2.T) for b1, b2 in zip(frames_1, frames_2)]
        return rm

    def get_theta(self, rotation_matrices):
        theta_array = [
            np.degrees(np.arccos((np.trace(x) - 1) / 2)) for x in rotation_matrices
        ]
        return theta_array

    def get_eigenvectors_nonsymmetric(self, rotation_matrices):
        w = []
        for x in rotation_matrices:
            w.append([x[1, 2] - x[2, 1], x[2, 0] - x[0, 2], x[0, 1] - x[1, 0]])
        return w

    def get_rotation_matrices_rodrigues(self, rotation_angle, u):
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

    def gram_schmidt_columns(self, matrix):
        n = matrix.shape[1]
        for j in range(n):
            for k in range(j):
                matrix[:, j] -= np.dot(matrix[:, k], matrix[:, j]) * matrix[:, k]
            matrix[:, j] = matrix[:, j] / np.linalg.norm(matrix[:, j])
        return matrix

    def get_translational_coords(self, lambda_vec, middle_frame):
        coor1, coor2, coor3 = ([] for i in range(3))

        for vec, matrix in zip(lambda_vec, middle_frame):
            x, y, z = np.dot(vec, matrix)
            coor1.append(x)
            coor2.append(y)
            coor3.append(z)

        return coor1, coor2, coor3

    def get_rotational_coords(self, theta, unit_vector, middle_frame):
        coor1, coor2, coor3 = ([] for i in range(3))

        for angle_t, u_vector, mid_frame in zip(theta, unit_vector, middle_frame):
            x, y, z = np.dot([angle_t * u for u in u_vector], mid_frame)
            coor1.append(x)
            coor2.append(y)
            coor3.append(z)

        return coor1, coor2, coor3

    def get_middle_frames(self, frames_1, frames_2, origins_1, origins_2):
        middle_frames = [
            (frame1 + frame2) / 2 for frame1, frame2 in zip(frames_1, frames_2)
        ]
        origins_middle_frames = [
            (origin1 + origin2) / 2 for origin1, origin2 in zip(origins_1, origins_2)
        ]
        middle_frames = [mf / np.linalg.norm(mf) for mf in middle_frames]
        middle_frames = [self.gram_schmidt_columns(mf) for mf in middle_frames]

        return middle_frames, origins_middle_frames

    def get_intra_coords(self):
        intra_rotation_matrices = self.get_rotation_matrices(
            self.frames_1, self.frames_2
        )
        theta_a = self.get_theta(intra_rotation_matrices)
        lambda_a = [(x2 - x1) for x1, x2 in zip(self.origins_1, self.origins_2)]
        w = self.get_eigenvectors_nonsymmetric(intra_rotation_matrices)
        unit_rotation_vector_a = [[v / np.linalg.norm(vec) for v in vec] for vec in w]
        shear, stagger, stretch = self.get_translational_coords(
            lambda_a, self.intra_middle_frames
        )
        propeller, buckle, opening = self.get_rotational_coords(
            theta_a, unit_rotation_vector_a, self.intra_middle_frames
        )

        return shear, stagger, stretch, propeller, buckle, opening

    def get_inter_coords(self):
        inter_rotation_matrices = self.get_rotation_matrices(
            self.intra_middle_frames, self.intra_middle_frames[1:]
        )
        theta_e = self.get_theta(inter_rotation_matrices)
        w_e = self.get_eigenvectors_nonsymmetric(inter_rotation_matrices)
        unit_rotation_vector_e = [[v / np.linalg.norm(vec) for v in vec] for vec in w_e]
        lambda_e = [
            (x2 - x1)
            for x1, x2 in zip(
                self.intra_middle_frames_origins, self.intra_middle_frames_origins[1:]
            )
        ]
        shift, slide, rise = self.get_translational_coords(
            lambda_e, self.inter_middle_frames
        )
        roll, tilt, twist = self.get_rotational_coords(
            theta_e, unit_rotation_vector_e, self.inter_middle_frames
        )

        return shift, slide, rise, roll, tilt, twist

    def run(self):
        # intra
        (
            self.intra_middle_frames,
            self.intra_middle_frames_origins,
        ) = self.get_middle_frames(
            self.frames_1, self.frames_2, self.origins_1, self.origins_2
        )

        (
            self.shear,
            self.stretch,
            self.stagger,
            self.buckle,
            self.propeller,
            self.opening,
        ) = self.get_intra_coords()

        # inter
        (
            self.inter_middle_frames,
            self.inter_middle_frames_origins,
        ) = self.get_middle_frames(
            self.intra_middle_frames,
            self.intra_middle_frames[1:],
            self.intra_middle_frames_origins,
            self.intra_middle_frames_origins[1:],
        )

        (
            self.shift,
            self.slide,
            self.rise,
            self.tilt,
            self.roll,
            self.twist,
        ) = self.get_inter_coords()
