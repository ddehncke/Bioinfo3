import argparse
import sys
import os
import Bio

def main():
    parser = argparse.ArgumentParser(description='pdb files')
    parser.add_argument('-i', dest="folder", metavar="FILE", help='input file', required=True)

    args = parser.parse_args()

    all_paths = os.listdir(args.folder)

    for filename in all_paths:
        print(filename)

main()