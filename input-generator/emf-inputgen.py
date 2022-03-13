"""
Author: Jessica M. Gonzalez-Delgado
        North Carolina State University

This script creates an input file for the eMF software from individual
relaxation data files.

The relaxation data files have to have the following format:
    Line 1: Type of relaxation data (values = noe, r1, r2)
    Line 2: Magnetic field (MHz)
    Lines 3+: residue_number value error cluster_number

For the distributed model, the distribution file name must be manually
added after creating the input file.

The relaxation data files must be named as:
    *r1.dat, *r2.dat, *noe.dat

Run as: python emf-inputgen.py
"""

import os
import sys


# get names of relaxation data files from current directory
def get_data_files():
    # get all data files in directory
    buffer = os.listdir(os.getcwd())
    files = []

    # keep only data files ending in r1.dat, r2.dat, or noe.dat
    for i in range(0, len(buffer)):
        if buffer[i][-4:] == '.dat':
            files.append(buffer[i])

    return files


# extract data from data files
def get_data(datafiles):
    files = get_data_files()
    fields = []
    noe = []
    r1 = []
    r2 = []

    for file in files:
        exp_flag = [False]
        with open(file) as fo:
            for line in fo:
                line = line.strip().split()
                if line[0][0] == "#":
                    continue
                else:
                    # check that the data file has the correct format
                    if len(line) > 1 and len(line) < 4 or len(line) > 4:
                        print(file, "\nError: Data file has incorrect format.\n"
                              "Exiting now...\n")
                        sys.exit()

                    # get relaxation data type
                    elif len(line) == 1 and exp_flag[0] is False:
                        if line[0] == "noe":
                            exp_flag[0] = True
                            exp_flag.append('noe')
                        elif line[0] == "r1":
                            exp_flag[0] = True
                            exp_flag.append('r1')
                        elif line[0] == "r2":
                            exp_flag[0] = True
                            exp_flag.append('r2')
                        else:
                            print(file, "\nError: Experiment type missing or "
                                  "incorrect file format.\nExiting now...")
                            sys.exit()

                    # get field info
                    elif len(line) == 1 and exp_flag[0] is True:
                        try:
                            fields.append(float(line[0]))
                        except ValueError:
                            print(file, "\nError: Experiment type missing or "
                                  "incorrect file format.\nExiting now...")
                            sys.exit()

                    # get relaxation data and edit data arrays for convenience
                    else:
                        line.insert(1, line.pop())
                        line.insert(2, fields[-1])
                        if exp_flag[1] == "noe":
                            noe.append(line)
                        elif exp_flag[1] == "r1":
                            r1.append(line)
                        elif exp_flag[1] == "r2":
                            r2.append(line)

    return noe, r1, r2


# check that all data files are the same length
def check_data(noe, r1, r2):
    if len(noe) != len(r1) or len(noe) != len(r2):
        print("Error: Amount of residues in data files do not match.\n"
              "Verify that relaxation data is present for the same residues "
              "in all files.\n")
        sys.exit()


# sort data arrays by residues in ascending order
def sort_data(noe, r1, r2):
    for i in range(0, len(noe)):
        noe[i][0] = int(noe[i][0])
        r1[i][0] = int(r1[i][0])
        r2[i][0] = int(r2[i][0])
    noe.sort()
    r1.sort()
    r2.sort()

    # make sure residues match between data sets
    for i in range(0, len(noe)):
        if noe[i][0] != r1[i][0] or noe[i][0] != r2[i][0]:
            print("Error: Residue mismatch between data files.\n"
                  "Verify that relaxation data is present for the same"
                  "residues in all files.")
            sys.exit()

    return noe, r1, r2


# write input file
def write_input_file(noe, r1, r2):
    input_file = "input.conf"
    header = \
'''# column 1   : residue number
# column 2   : cluster number, used for inclusion or exclusion of data
# column 3   : B0, field strength in MHz
# column 4-5 : R1 relaxation rate and error (1/s)
# column 6-7 : R2 relaxation rate and error (1/s)
# column 8-9 : {1H}-X NOE and error
# column 10  : (only for distributed) tc distribution file name \n
'''

    for i in range(0, len(noe)):
        r2[i] = r2[i][-2:]
        noe[i] = noe[i][-2:]
        r1[i] += r2[i] + noe[i]

    fo = open(input_file, 'w')
    fo.write(header)
    for i in range(0, len(r1)):
        buffer = ""
        # aesthetics - not perfect but it does the trick
        if r1[i][0] < 10:
            buffer += str(r1[i][0]) + "    "
            for j in range(1, len(r1[0])-1):
                buffer += str(r1[i][j])
                buffer += "    "
            buffer += '\n'
        elif r1[i][0] < 100:
            buffer += str(r1[i][0]) + "   "
            for j in range(1, len(r1[0])-1):
                buffer += str(r1[i][j])
                buffer += "    "
            buffer += '\n'
        else:
            buffer += str(r1[i][0]) + "  "
            for j in range(1, len(r1[0])-1):
                buffer += str(r1[i][j])
                buffer += "    "
            buffer += '\n'

        fo.write(buffer)
    fo.close()


# MAIN PROGRAM
def main():
    dfiles = get_data_files()
    data = get_data(dfiles)
    check_data(data[0], data[1], data[2])
    data = sort_data(data[0], data[1], data[2])
    write_input_file(data[0], data[1], data[2])


# RUN PROGRAM
main()
