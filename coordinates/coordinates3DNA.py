from coordinates.coordinates import Coordinates
import math
import numpy as np
import pandas as pd


class Coordinates3DNA(Coordinates):

    def get_hinge_axis_array(self):
        axis_array = np.zeros(shape=(self.intra_len, 3))
        for i in range(self.intra_len):
            axis = np.cross(self.frames_1[i].T[2], self.frames_2[i].T[2])
            axis_normalized = axis / np.linalg.norm(axis)
            axis_array[i] = axis_normalized
        return axis_array

    def get_angle_array(self, col):
        angle_array = np.zeros(shape=self.intra_len)
        for i in range(self.intra_len):
            angle = np.dot(self.frames_1[i].T[col], self.frames_2[i].T[col]) / (
                    np.linalg.norm(self.frames_1[i].T[col]) * np.linalg.norm(self.frames_2[i].T[col]))
            angle = np.degrees(np.arccos(angle))
            angle_array[i] = angle
        return angle_array

    @staticmethod
    def create_rotation_matrix(angle, hinge_axis):
        r = pd.DataFrame([
            [math.cos(math.radians(angle)) + (1 - math.cos(math.radians(angle))) * (hinge_axis[0] ** 2),
             (1 - math.cos(math.radians(angle))) * hinge_axis[0] * hinge_axis[1] - hinge_axis[2] * math.sin(
                 math.radians(angle)),
             (1 - math.cos(math.radians(angle))) * hinge_axis[0] * hinge_axis[2] + hinge_axis[1] * math.sin(
                 math.radians(angle))],
            [(1 - math.cos(math.radians(angle))) * hinge_axis[0] * hinge_axis[1] + hinge_axis[2] * math.sin(
                math.radians(angle)),
             math.cos(math.radians(angle)) + (1 - math.cos(math.radians(angle))) * (hinge_axis[1] ** 2),
             (1 - math.cos(math.radians(angle))) * hinge_axis[1] * hinge_axis[2] - hinge_axis[0] * math.sin(
                 math.radians(angle))],
            [(1 - math.cos(math.radians(angle))) * hinge_axis[0] * hinge_axis[2] - hinge_axis[1] * math.sin(
                math.radians(angle)),
             (1 - math.cos(math.radians(angle))) * hinge_axis[1] * hinge_axis[2] + hinge_axis[0] * math.sin(
                 math.radians(angle)),
             math.cos(math.radians(angle)) + (1 - math.cos(math.radians(angle))) * (hinge_axis[2] ** 2)]
        ])
        return r

    def get_rotation_matrix_array(self, multiply, axis_array, angle_array):
        rotation_matrices = np.zeros(shape=(len(axis_array), 3, 3))
        for i in range(len(axis_array)):
            r_matrix = self.create_rotation_matrix(angle_array[i] * multiply, axis_array[i])
            rotation_matrices[i] = r_matrix
        return rotation_matrices

    def rotate_frames(self, strand, rotation_matrices):
        rotated_matrices = np.zeros(shape=(self.intra_len, 3, 3))
        for i in range(len(strand)):
            rotated_matrices[i] = np.matmul(rotation_matrices[i], strand[i])
        return rotated_matrices

    def create_middle_frames_array(self, rotated_frames_1, rotated_frames_2):
        self.intra_middle_frames = np.zeros(shape=(self.intra_len, 3, 3))
        for i in range(self.intra_len):
            self.intra_middle_frames[i] = np.mean([rotated_frames_1[i], rotated_frames_2[i]], axis=0)

    def create_middle_frames_origins_array(self):
        self.intra_middle_frames_origins = np.zeros(shape=(self.intra_len, 3))
        for i in range(self.intra_len):
            self.intra_middle_frames_origins[i] = np.mean([self.origins_1[i], self.origins_2[i]], axis=0)

    def get_translational_coordinates(self):
        translational_coordinates = np.zeros(shape=(self.intra_len, 3))
        for i in range(self.intra_len):
            translational_coordinates[i] = np.matmul((self.origins_2[i] - self.origins_1[i]),
                                                     self.intra_middle_frames[i])
        return translational_coordinates

    def get_opening_angle_array(self, rotated_1, rotated_2):
        angle_array = np.zeros(shape=self.intra_len)
        for i in range(self.intra_len):
            angle = np.dot(rotated_1[i].T[1], rotated_2[i].T[1]) / (
                    np.linalg.norm(rotated_1[i].T[1]) * np.linalg.norm(rotated_2[i].T[1]))
            opening_sign = np.dot(np.cross(rotated_1[i].T[1],
                                           rotated_2[i].T[1]),
                                  self.intra_middle_frames[i].T[2])
            angle = np.degrees(np.arccos(angle))
            angle_array[i] = angle * np.sign(opening_sign)
        return angle_array

    @staticmethod
    def get_phase_angle_array(hinge_axis_array, middle_frames):
        angle_array = np.zeros(shape=len(hinge_axis_array))
        for i in range(len(hinge_axis_array)):
            angle = np.dot(hinge_axis_array[i], middle_frames[i].T[1]) / (
                    np.linalg.norm(hinge_axis_array[i]) * np.linalg.norm(middle_frames[i].T[1]))
            angle_sign = np.dot(np.cross(hinge_axis_array[i],
                                         middle_frames[i].T[1]),
                                middle_frames[i].T[2])
            angle = np.degrees(np.arccos(angle))
            angle_array[i] = angle * np.sign(angle_sign)
        return angle_array

    @staticmethod
    def get_rotational_coordinates_array(angle, phase_angle):
        coordinates_array = np.zeros(shape=(len(angle), 2))
        for i in range(len(angle)):
            coordinates = [angle[i] * math.cos(np.radians(phase_angle[i])),
                           angle[i] * math.sin(np.radians(phase_angle[i]))]
            coordinates_array[i] = coordinates
        return coordinates_array

    def compute_intra_coordinates(self):
        hinge_axes = self.get_hinge_axis_array()
        buckle_propeller_angle_array = self.get_angle_array(2)

        rotation_matrices_2 = self.get_rotation_matrix_array((-0.5), hinge_axes, buckle_propeller_angle_array)
        rotation_matrices_1 = self.get_rotation_matrix_array(0.5, hinge_axes, buckle_propeller_angle_array)

        rotated_base_frames_strand_2 = self.rotate_frames(self.frames_2, rotation_matrices_2)
        rotated_base_frames_strand_1 = self.rotate_frames(self.frames_1, rotation_matrices_1)

        self.create_middle_frames_array(rotated_base_frames_strand_1, rotated_base_frames_strand_2)
        self.create_middle_frames_origins_array()

        self.shear, self.stretch, self.stagger = np.transpose(self.get_translational_coordinates())
        self.opening = self.get_opening_angle_array(rotated_base_frames_strand_1, rotated_base_frames_strand_2)

        phase_array_intra = self.get_phase_angle_array(hinge_axes, self.intra_middle_frames)

        self.propeller, self.buckle = np.transpose(self.get_rotational_coordinates_array(buckle_propeller_angle_array,
                                                                                         phase_array_intra))

    def get_hinge_axes_array_inter(self):
        hinge_axes_array = np.zeros(shape=(self.inter_len, 3))
        for i in range(self.inter_len):
            hinge_axis = np.cross(self.intra_middle_frames[i].T[2], self.intra_middle_frames[i + 1].T[2])
            hinge_axis_normalized = hinge_axis / np.linalg.norm(hinge_axis)
            hinge_axes_array[i] = hinge_axis_normalized
        return hinge_axes_array

    def get_angle_array_inter(self, col):
        angle_array = np.zeros(shape=self.inter_len)
        for i in range(self.inter_len):
            angle = np.dot(self.intra_middle_frames[i].T[col], self.intra_middle_frames[i + 1].T[col]) / (
                    np.linalg.norm(self.intra_middle_frames[i].T[col]) *
                    np.linalg.norm(self.intra_middle_frames[i + 1].T[col]))
            angle = np.degrees(np.arccos(angle))
            angle_array[i] = angle
        return angle_array

    def rotate_basepair_frames(self, angle, axis, index):
        r1_matrix = self.create_rotation_matrix(angle * 0.5, axis)
        rotated_frame_1 = np.matmul(r1_matrix, self.intra_middle_frames[index])
        r2_matrix = self.create_rotation_matrix((angle * (-0.5)), axis)
        rotated_frame_2 = np.matmul(r2_matrix, self.intra_middle_frames[index + 1])

        return rotated_frame_1, rotated_frame_2

    def get_middle_basepair_frames(self, angle_array, hinge_axes_array):
        for i in range(self.inter_len):
            rotated_frame_1, rotated_frame_2 = self.rotate_basepair_frames(angle_array[i], hinge_axes_array[i], i)

            self.inter_middle_frames.append(np.mean([rotated_frame_1, rotated_frame_2], axis=0))

    def get_middle_basepair_frames_origins(self):
        for i in range(self.inter_len):
            self.inter_middle_frames_origins.append(np.mean([self.intra_middle_frames_origins[i],
                                                             self.intra_middle_frames_origins[i + 1]], axis=0))

    def get_basepair_translational_coordinates(self):
        coordinates = np.zeros(shape=(self.inter_len, 3))
        for i in range(self.inter_len):
            coordinates[i] = np.matmul(
                (self.intra_middle_frames_origins[i + 1] - self.intra_middle_frames_origins[i]).T,
                self.inter_middle_frames[i])
        return coordinates

    def get_twist_angle_array(self, angle_array, hinge_axes_array):
        twist_angles_array = np.zeros(shape=self.inter_len)
        for i in range(self.inter_len):
            rotated_frame_1, rotated_frame_2 = self.rotate_basepair_frames(angle_array[i], hinge_axes_array[i], i)
            twist_angle = np.dot(rotated_frame_1[1], rotated_frame_2[1]) / (
                    np.linalg.norm(rotated_frame_1[1]) * np.linalg.norm(rotated_frame_2[1]))
            twist_sign = np.dot(np.cross(rotated_frame_1[1], rotated_frame_2[1]), self.inter_middle_frames[i].T[2])
            twist_angle = np.degrees(np.arccos(twist_angle))
            twist_angles_array[i] = twist_angle * np.sign(twist_sign)
        return twist_angles_array

    def compute_inter_coordinates(self):
        hinge_axes = self.get_hinge_axes_array_inter()
        roll_tilt_angle_array = self.get_angle_array_inter(2)
        self.get_middle_basepair_frames(roll_tilt_angle_array, hinge_axes)
        self.get_middle_basepair_frames_origins()
        self.shift, self.slide, self.rise = np.transpose(self.get_basepair_translational_coordinates())

        self.twist = self.get_twist_angle_array(roll_tilt_angle_array, hinge_axes)
        phase_array_inter = self.get_phase_angle_array(hinge_axes, self.inter_middle_frames)
        self.roll, self.tilt = np.transpose(self.get_rotational_coordinates_array(roll_tilt_angle_array,
                                                                                  phase_array_inter))

    def run(self):
        self.compute_intra_coordinates()
        self.compute_inter_coordinates()
