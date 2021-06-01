import click

from coordinates.coordinates3DNA import Coordinates3DNA
from coordinates.coordinatesCgDNA import CoordinatesCgDNA
from coordinates.coordinatesCurves import CoordinatesCurves
from coordinates.frameFitter import FrameFitter


def attach_number_to_input_filename(value):
    """Callback function takes value of input file path
    and appends 1 before checking if path exists."""
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
@click.option('-o', '--output-file',
              required=True, type=str,
              help='Path to output file.')
@click.option('-n', '--number-of-snapshots', required=True, type=int, help='Number of snapshots to be processed.')
def compute_coordinates(coords_type, input_file, output_file, number_of_snapshots):
    """Takes the command lien arguments coordinate type, input file path,
    output files path and number of snapshots to be processed.
    Computes the internal coordinates according to definition type selected."""
    input_file = input_file[:-1]
    for i in range(1, number_of_snapshots + 1):
        ff = FrameFitter()
        ff.run(input_file + str(i))
        if coords_type == "3dna":
            coordinates = Coordinates3DNA(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
        elif coords_type == "curves":
            coordinates = CoordinatesCurves(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
        elif coords_type == "cgdna":
            coordinates = CoordinatesCgDNA(ff.frames_1, ff.frames_2, ff.origins_1, ff.origins_2, ff.strand_len)
        else:
            raise ValueError("Invalid coordinate type.")

        coordinates.run()
        coordinates.round_coords()
        coordinates.write_to_bp_output_files(output_file)
        coordinates.write_to_bp_step_output_files(output_file)


compute_coordinates()
