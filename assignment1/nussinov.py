import argparse
import sys
import os

# initialize the DP matrix

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

# helper function so i don't mess up the indices
rangeInclusive = lambda start, end: range(start, end + 1)
rangeExclusive = lambda start, end: range(start + 1, end)

def main():

    parser = argparse.ArgumentParser(description='imput fasta file')
    parser.add_argument('-i', dest="filename", metavar="FILE", help='input file', required=True,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("--min-loop-length", type=int, dest='min_loop_length', default=1)
    parser.add_argument("-score_AU", type=int, dest='score_AU', default=10)
    parser.add_argument("-score_CG", type=int, dest='score_CG', default=1)
    parser.add_argument("-score_GU", type=int, dest='score_GU', default=0)

    args = parser.parse_args()
    # mein temp kommentar

    score_AU = args.score_AU
    score_CG = args.score_CG
    score_GU = args.score_GU
    min_loop_length = args.min_loop_length
    input_file = args.filename

    sequence = ''

    with open('input_files/supplement01/test.fasta', 'r') as file_handle:
        for line in file_handle:
            line = line.strip()
            if line[0] == '>':
                continue
            sequence = line.strip()

    # transform seuence to dict because of indices
    dic_sequence = {}
    for i in rangeInclusive(1,len(sequence)):
        dic_sequence[i] = sequence[i-1]

    dp_dic = {}
    traceback_dic = {}

    # create 2 dimensional dictionatries for dp and traceback
    for i in rangeInclusive(1,len(sequence)):
        dp_dic[i] = {}
        traceback_dic[i] = {}

    # fill with 0
    for i in rangeInclusive(2,len(sequence)):
        dp_dic[i][i-1] = 0
    for i in rangeInclusive(1,len(sequence)):
        dp_dic[i][i] = 0

    dic_score = {('A','U'):score_AU, ('U','A'):score_AU,
                 ('C','G'):score_CG, ('G','C'):score_CG,
                 ('G','U'):score_AU, ('U','G'):score_AU,
                 ('U', 'U'): 0, ('G', 'G'): 0,
                 ('A', 'A'): 0, ('C', 'C'): 0,
                 ('U', 'C'): 0, ('C', 'U'): 0,
                 ('A', 'C'): 0, ('C', 'A'): 0,
                 ('A', 'G'): 0, ('G', 'A'): 0}

    fillMatrices(dp_dic, dic_score, traceback_dic, dic_sequence)

def fillMatrices(dp_dic, dic_score, traceback_dic, dic_sequence):

    def getPairScore(i, j):
        return dic_score[(dic_sequence[i], dic_sequence[j])]

    # initialization recursion
    for I in rangeInclusive(2, len(dic_sequence)):
        for j in rangeInclusive(I, len(dic_sequence)):
            i = j - I + 1

            case1 = dp_dic[i + 1][j]
            case2 = dp_dic[i][j - 1]

            case3 = dp_dic[i + 1][j - 1] + getPairScore(i, j)

            case4 = -1
            temp_trace = ()
            # TODO hier noch minimum loop length implementieren
            for k in rangeExclusive(i,j):
                if dp_dic[i][k] + dp_dic[k+1][j] > case4:
                    temp_trace = ((i,k),(k+1,j))

            if case1 > case2 and case1 > case3 and case1 > case4:
                dp_dic[i][j] = case1
                traceback_dic[i][j] = ((i + 1,j))
            elif case2 > case3 and case2 > case4:
                dp_dic[i][j] = case2
                traceback_dic[i][j] = ((i, j - 1))
            elif case4 > case3: # switched here with case 4 to not confuse stuff
                dp_dic[i][j] = case4
                traceback_dic[i][j] = temp_trace
            else:
                dp_dic[i][j] = case3
                traceback_dic[i][j] = ((i + 1,j - 1))



    def traceback(i,j):
        if i < j:
            if dp_dic[i][j] == dp_dic[i+1][j]:
                traceback(i + 1, j)
            elif dp_dic[i][j] == dp_dic[i][j-1]:
                traceback(i, j - 1)
            elif dp_dic[i][j] == dp_dic[i + 1][j - 1] + getPairScore(i, j):
                print(dic_sequence[i],dic_sequence[j])
                traceback(i + 1, j - 1)
            else:
                for k in rangeInclusive(i+1,j-1):
                    if dp_dic[i][j] == dp_dic[i][k] + dp_dic[k+1][j]:
                        traceback(i, k)
                        traceback(k + 1,j)
                        break

    traceback(1, 30)



def printMatrix(dp_dic):
    for k, v in dp_dic.items():
        print(v)




main()