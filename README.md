# Efficient computation of rigid base coordinates in Python

Python module for computation of internal coordinates in DNA and RNA. Takes data from MD simulations in the form of a PDB file and creates output files containing base  coordiantes *shear, stretch, stagger, buckle, propeller, opening* and base step coordinates *shift, slide, rise, roll, tilt, twist*. The coordinates can be computed in either of the three definitions - 3DNA, Curves+ or cgDNA.

This project was created as a a part of my bachelor thesis.


# Requirements

The module can be run with Python 3.9.0 interpreter. Requirements are listed in requirements.txt file and can be installed simply using

    pip -r requirements.txt

# Usage
The usage of the module can be displayed using `--help`

    Usage: python -m coordinates [OPTIONS]
    
      Takes the command lien arguments coordinate type, input file path, output
      files path and number of snapshots to be processed. Computes the internal
      coordinates according to definition type selected.
    
    Options:
      -t, --coords-type [3dna|curves|cgdna]
                                      Type of coordinates definition.  [required]
      -i, --input-file PATH           Path to input pdb file.  [required]
      -o, --output-file TEXT          Path to output file.  [required]
      -n, --number-of-snapshots INTEGER
                                      Number of snapshots to be processed.
                                      [required]
      --help                          Show this message and exit.

It is meant to be run as

    python -m coordinates -t definitio_type -n number_of_pdb_files -i /input/path/filename -o /output/path/filename

where input file name is a common name for a PDB files numbered 1 to n. And the output file name is a common name for the output files. For example

    python -m coordinates -t 3dna -n 10 -i /my_dir/simulation -o /my_other_dir/coordinates

would compute base and base step coordinates according to the 3DNA definition for PDB files named simulation.pdb.1 through simulation.pdb.10.

# Output files
The output files are created for base and base step coordinates separately. For base coordinates, each base pair in the simulated oligomer will have its own output file, where tab separated columns represent *shear, stretch, stagger, buckle, propeller and opening* and rows represent its coordinates through the snapshots. The same goes for base step coordinates, each pair of base pairs will have its own file with columns *shift, slide, rise, roll, tilt, twist* and rows with coordinates for every snapshot. In the example above for an oligomer of length 5 with 10 snapshots the output files created are coordiantes_base_prm_1.out through coordiantes_base_prm_5.out each with 10 rows and coordinates_base_step_1.out through coordinates_base_step_4.out also with 10 rows.
