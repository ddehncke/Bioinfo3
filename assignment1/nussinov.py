import argparse
import sys
import os

# initialize the DP matrix

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

def main():
    parser = argparse.ArgumentParser(description='imput fasta file')
    parser.add_argument('-i', dest="filename", metavar="FILE", help='input file', required=True,type=lambda x: is_valid_file(parser, x))

    args = parser.parse_args()

    input_file = args.filename

    sequence = ''

    with open('input_files/supplement01/test.fasta', 'r') as file_handle:
        for line in file_handle:
            line = line.strip()

            if line[0] == '>':
                continue

            sequence = line.strip()

    print(sequence)

main()