import argparse
import sys
import os
from collections import OrderedDict
from operator import itemgetter
from Bio.PDB import *

dic_helix = {'CYS': 0, 'ASP': 0, 'SER': 0, 'GLN': 0, 'LYS': 0,
     'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0,
     'GLY': 0, 'HIS': 0, 'LEU': 0, 'ARG': 0, 'TRP': 0,
     'ALA': 0, 'VAL': 0, 'GLU': 0, 'TYR': 0, 'MET': 0, 'MSE': 0}

dic_sheet = {'CYS': 0, 'ASP': 0, 'SER': 0, 'GLN': 0, 'LYS': 0,
     'ILE': 0, 'PRO': 0, 'THR': 0, 'PHE': 0, 'ASN': 0,
     'GLY': 0, 'HIS': 0, 'LEU': 0, 'ARG': 0, 'TRP': 0,
     'ALA': 0, 'VAL': 0, 'GLU': 0, 'TYR': 0, 'MET': 0, 'MSE': 0}

def main():
    parser = argparse.ArgumentParser(description='pdb files')
    parser.add_argument('-i', dest="folder", metavar="FILE", help='input folder', required=True)

    args = parser.parse_args()

    all_paths = os.listdir(args.folder)

    # parse the pdb file

    alpha_helix_count = 0
    three_ten_helix_count = 0
    sheet_count = 0
    nr_all_chains = 0
    i = 0

    for filename in all_paths:
        # print(i)
        i += 1
        parser = PDBParser()
        path_to_current_file = os.path.join(args.folder, filename)
        structure = parser.get_structure('temp', path_to_current_file)

        # print(filename)
        nr_Conformations = 0
        current_Conformation = ''
        aa_ary = []

        # additionally read in the secondary structure
        read_secondary_structure = open(path_to_current_file, 'r')
        for line in read_secondary_structure:
            ary_line = line.split()

            # count residues in file
            if line[0:6] == 'SEQRES':
                if len(ary_line[4]) < 2:
                    continue
                else:
                    aa_ary.extend(ary_line[4:])

            if line[0:5] == 'HELIX':
                # need to keep track of conformations to correctly compute the percentage. Will not
                # be computed for sheets as this would be double times
                if line[19] != current_Conformation:
                    current_Conformation = line[19]
                    nr_Conformations += 1

                # print(line[72:76])
                # alpha helix
                if line[39:40].strip() == '1':
                    alpha_helix_count += int(line[72:76].strip())

                    # the official residues in the pdb and the aa in the aminoacid dont match?
                    # print('start aa: ',line[15:18])
                    # print('end aa',line[27:30])
                    # print(start, end)
                    # print(aa_ary)

                    start = int(line[22:25])
                    end = int(line[34:37])

                    for aa in range(start, end):
                        if end < len(aa_ary):
                            if dic_helix.get(aa_ary[aa]) is not None:
                                dic_helix[aa_ary[aa]] = dic_helix[aa_ary[aa]] + 1


                # threeTen Helix
                if line[39:40].strip() == '5':
                    three_ten_helix_count += int(line[72:76].strip())
            else:
                if line[0:5] == 'SHEET':
                    sheet_count += int(line[33:37].strip()) - int(line[22:26].strip())

                    start = int(line[22:26])
                    end = int(line[34:37])
                    # print(start, end)
                    # print(len(ary_line))


                    for aa in range(start, end):
                        #  skipping ranges of sheets that are obviously not possible: end longer than sequence
                        if end < len(aa_ary):
                            if dic_sheet.get(aa_ary[aa]) is not None:
                                dic_sheet[aa_ary[aa]] = dic_sheet[aa_ary[aa]] + 1

        # all different secondary structures in file have been counted, therefore
        # the cain hast to be counted multiple times
        nr_all_chains += len(aa_ary)

    # print_a_and_b(alpha_helix_count, three_ten_helix_count, sheet_count, nr_all_chains, dic_helix, dic_sheet)

    # for model in structure:
    #     for chain in model:
    #         residues = chain.get_residues()
    #         # for i in range(0,len(residues)):
    #         #     print(residues[i])
    #         for residue in chain:
    #             id = residue.get_full_id()
    #             print(id[3][1])
    #             # for atom in residue:
    #                 # print(atom.get_parent())

    #
    # for model in structure:
    #     for chain in model:
    #         residues = chain.get_residues()
    #         # for i in range(0,len(residues)):
    #         #     print(residues[i])
    #         for residue in chain:
    #             id = residue.get_full_id()
    #             print(id[3][1])
    #             # for atom in residue:
    #                 # print(atom.get_parent())



def print_a_and_b(alpha_helix_count, three_ten_helix_count, sheet_count, nr_all_chains, dic_helix, dic_sheet):
    print('The percentages of the secondary structures measures by the sequences.')
    print('alphaHelixContent:',alpha_helix_count/nr_all_chains)
    print('three_ten_helix_count:', three_ten_helix_count / nr_all_chains)
    print('sheet_count:', sheet_count / nr_all_chains)

    print('The number of Aminoacids in the present helix')
    dic_helix_sorted = OrderedDict(sorted(dic_helix.items(), key=itemgetter(1)))
    for i in dic_helix_sorted:
        print(i, dic_helix_sorted[i])


    print('The number of Aminoacids in the present sheets')
    dic_sheet_sorted = OrderedDict(sorted(dic_sheet.items(), key=itemgetter(1)))
    for i in dic_sheet_sorted:
        print(i, dic_sheet_sorted[i])



main()