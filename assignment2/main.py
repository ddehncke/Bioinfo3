import argparse
import sys
import os
from Bio.PDB import *



def main():
    parser = argparse.ArgumentParser(description='pdb files')
    parser.add_argument('-i', dest="folder", metavar="FILE", help='input folder', required=True)

    args = parser.parse_args()

    all_paths = os.listdir(args.folder)

    # parse the pdb file

    alpha_helix_count = 0
    three_ten_helix_count = 0
    sheet_count = 0
    chain_length_counter = 0
    i = 0

    for filename in all_paths:
        print(i)
        i += 1
        parser = PDBParser()
        path_to_current_file = os.path.join(args.folder, filename)
        structure = parser.get_structure('temp', path_to_current_file)

        print(filename)

        # additionally read in the secondary structure
        read_secondary_structure = open(path_to_current_file, 'r')
        for line in read_secondary_structure:
            ary_line = line.split()

            if ary_line[0] == 'SEQRES':
                chain_length_counter += len(ary_line) - 2

            if len(ary_line) < 6:
                continue

            # only count for first conformation
            if ary_line[4] == 'A' or ary_line[5] == 'A':
                if ary_line[0] == 'HELIX':

                    if line[38:40].strip() == '1':
                        try:
                            alpha_helix_count += int(ary_line[8]) - int(ary_line[5])
                        except ValueError:
                            end = [int(s) for s in ary_line[8].split() if s.isdigit()]
                            start = [int(s) for s in ary_line[5].split() if s.isdigit()]
                            alpha_helix_count += end - start
                    if line[38:40].strip() == '5':
                        three_ten_helix_count += int(ary_line[8]) - int(ary_line[5])
                if ary_line[0] == 'SHEET':
                    try:
                        end = int(ary_line[9])
                        start = int(ary_line[6])
                    except ValueError:
                        end = int(ary_line[8])
                        start = int(ary_line[5])
                    except IndexError:
                        end = int(ary_line[8])
                        start = int(ary_line[5])
                    sheet_count += end - start

        # print(chain_length_counter, alpha_helix_count, three_ten_helix_count, sheet_count)
    print('alphaHelixContent:',alpha_helix_count/chain_length_counter)
    print('three_ten_helix_count:', three_ten_helix_count / chain_length_counter)
    print('sheet_count:', sheet_count / chain_length_counter)


main()