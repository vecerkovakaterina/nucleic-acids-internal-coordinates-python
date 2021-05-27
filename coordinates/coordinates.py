from abc import ABC, abstractmethod


class Coordinates(ABC):
    def __init__(self, frames_1, frames_2, origins_1, origins_2, strand_len):
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

    def write_to_output_files(self):
        pass
