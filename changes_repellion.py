#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# changes_repellion.py
#
# 28 de febr. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def preprocess_rmout(
    input_rmout: str,
    output_modified_rmout: str, ):
    """
    Preprocess RM.out files. The 16th field is inconsistent; some rows contain
    an asterisk, while other rows have an empty 16th field. Make all rows
    consistent (asterisk -> True, empty -> False).

    Input
    -----

    + input_rmout: String pointing to a file.

    + output_modified_rmout: String where a new file will be written to. If
      there exists a file with the same filename, it will be overwritten.
    """
    with open(input_rmout) as inp, open(output_modified_rmout, "w") as out:
        # Skip first three lines (strangely formatted header).
        for i in range(3):
            out.write(inp.readline())
        # Process the rest of the file.
        for num, line in enumerate(inp):
            overlapping_bool = (line.split()[-1] == "*")
            out_line = str(line).strip("\n*") + "\t" +\
                str(overlapping_bool) + "\n"
            # Some rows might not have all 16 fields!
            if len(out_line.split()) != 16:
                print("ERROR: The amount of fields/columns in the line",
                      int(num)+4, "is not equal to 16. Make sure that "+
                      "this line is not missing any fields.")
                print("The problematic line is...")
                print(line)
                # Stop writing file if an error is encountered.
                return None
            out.write(out_line)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.


