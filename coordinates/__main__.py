import numpy as np

from coordinates.coordinates3DNA import Coordinates3DNA
from coordinates.frameFitter import FrameFitter

# todo nacitani argumentu
# todo komentare
# todo vytvareni souboru
if __name__ == "__main__":
    ff = FrameFitter()
    ff.run("/Users/katerina/Desktop/Bakalářská práce/fitovani/teplota.300.pdb.1")
    coords_type = Coordinates3DNA(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
    coords_type.run()
    [print(#np.round(a, 2),
           #np.round(b, 2),
           #np.round(c, 2),
           #np.round(d, 2),
           np.round(e, 2),
           #np.round(f, 2)
           ) for
     a, b, c, d, e, f in
     zip(coords_type.shift, coords_type.slide, coords_type.rise, coords_type.roll, coords_type.tilt,
         coords_type.twist)]
