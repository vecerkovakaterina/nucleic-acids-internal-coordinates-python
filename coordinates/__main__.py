import numpy as np
import click

from coordinates.coordinates import Coordinates
from coordinates.coordinates3DNA import Coordinates3DNA
from coordinates.coordinatesCgDNA import CoordinatesCgDNA
from coordinates.coordinatesCurves import CoordinatesCurves
from coordinates.frameFitter import FrameFitter


# todo nacitani argumentu
# todo komentare
# todo vytvareni souboru
def attach_number_to_input_filename(ctx, param, value):
    if value is not None:
        value += str(1)
    return value


@click.command()
@click.option('-t', '--coords-type', required=True,
              type=click.Choice(['3dna', 'curves', 'cgdna'], case_sensitive=False),
              help='Type of coordinates definition.')
@click.option('-i', '--input-file',
              callback=attach_number_to_input_filename, required=True, type=click.Path(),
              help='Path to input pdb file.')
@click.option('-n', '--number-of-snapshots', required=True, type=int, help='Number of snapshots to be processed.')
def compute_coordinates(coords_type, input_file, number_of_snapshots):
    input_file = input_file[:-1]
    ff = FrameFitter()
    for i in range(1, number_of_snapshots + 1):
        ff.run(input_file + str(i))
        if coords_type == "3dna":
            coords_type = Coordinates3DNA(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
            print("3dna")
        elif coords_type == "curves":
            coords_type = CoordinatesCurves(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
            print("curves")
        elif coords_type == "cgdna":
            coords_type = CoordinatesCgDNA(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
            print("cgdna")
        else:
            raise ValueError("Invalid coordinate type.")

        coords_type.run()
        coords_type.write_to_output_files()  # todo


compute_coordinates()

