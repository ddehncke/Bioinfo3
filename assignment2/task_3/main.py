#!/usr/bin/env python

import argparse
import sys
import os
from collections import OrderedDict
from operator import itemgetter
from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py

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
    dist_O_N = []
    dist_O_N_Helix = []

    for filename in all_paths:
        print(filename)
        path_to_current_file = os.path.join(args.folder, filename)
        parser = PDBParser()
        structure = parser.get_structure('temp', path_to_current_file)

        # add the distances from the current file
        dist_O_N.extend(getStructureDistancesO_N(path_to_current_file))

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

                    # computes the distance for the given Helix
                    if end < len(aa_ary):
                        id_chain = line[19:20]
                        # dist_O_N_Helix.extend(getHelixDistances(id_chain, start, end, structure))


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

    # histogram(dist_O_N_Helix, 'distance O-N in helices')
    histogram(dist_O_N, 'distance O-N in whole structure')

def getHelixDistances(id_chain, start, end, structure):
    list_dist = []

    # compute the distance between the fourth N and the first O using and counter
    wait_four_residues = 1
    stack_o = []
    for model in structure:
        for chain in model:
            for residue in chain:

                # checks if the current residue is bigger than the start of the helix and smaller
                # than the end of the helix. Then the distance to the fourth atom is computed.

                if start <= int(residue.get_full_id()[3][1]) < end and id_chain == chain.id:

                    # need to skip residues that are incomplete. Also need to skip conformations that got the same
                    # chain identifier.
                    if 'O' not in residue or 'N' not in residue:
                        wait_four_residues = 1
                        stack_o = []
                        continue

                    stack_o.append(residue)

                    if int(residue.get_full_id()[3][1]) < int(stack_o[0].get_full_id()[3][1]):
                        wait_four_residues = 1
                        stack_o = []
                        continue

                    if wait_four_residues >= 4:
                        temp = np.linalg.norm(stack_o.pop(0)['O'].get_vector() - residue['N'].get_vector())
                        list_dist.append(temp)

                    wait_four_residues += 1

    return list_dist


def getStructureDistancesO_N(path):
    temp_dist = []
    parser = PDBParser()
    structure = parser.get_structure('temp', path)
    # compute the distance between the fourth N and the first O using and counter
    wait_four_residues = 1
    stack_o = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # need to skip residues that are incorrect
                if 'O' not in residue or 'N' not in residue:
                    wait_four_residues = 1
                    stack_o = []
                    continue

                stack_o.append(residue)

                if int(residue.get_full_id()[3][1]) < int(stack_o[0].get_full_id()[3][1]):
                    wait_four_residues = 1
                    stack_o = []
                    continue

                if wait_four_residues >= 4:
                    # temp_dist.append(np.linalg.norm(stack_O.pop(0) - residue['N'].get_vector()))
                    temp_dist.append(np.linalg.norm(stack_o.pop(0)['O'].get_vector() - residue['N'].get_vector()))


                wait_four_residues += 1

    return temp_dist


def histogram(distances, name):
    # remove values that are to big. Dont know why, tried to fix
    a = np.array(distances)
    b = a[np.logical_and(a>=0, a<9.5)]
    print(len(b)/len(a))
    plt.hist(b)
    plt.title(name)
    plt.show()

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