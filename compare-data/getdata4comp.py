"""
Author: Jessica M. Gonzalez-Delgado
        North Carolina State University

This script creates a file containing relaxation rates that are available for
comparison between two files, e.g. the file will contain spin #1234 only if
data for spin #1234 is present in file #1, file #2 and file #3.

This is useful when you have multiple relaxation or eMF digested data files for
the spin system that might not contain data for the same spins and you want to
see the differences in the data files without having to manually weed out the
spins that are not present in both files.

The input data files can be in two formats:
    (1) Format required to run "emf-inputgen.py" from the eMF software
    (2) Format from the eMF digested data files


The output file is named 'data4comp.out' and has the following general format:
    Line 1: Header wit the names of the files used
    Lines 2+: spin-number     data-from-file-1    data-from-file-2

Run as: python getdata4comp.py file-1 file-2
"""

import sys


# check input files provided
def get_dat_files():

    files = sys.argv[1:]
    emf = False

    if len(files) == 0:
        print("Error! No DAT or OUT file specified.\nExiting now...")
        sys.exit()

    if len(files) == 1:
        print("Error! You need to provide two files.\nExiting now...")
        sys.exit()

    if len(files) > 2:
        print("Error! More than two files specified.\nExiting now...")
        sys.exit()

    if len(files) > 1 and len(files) < 3:
        for file in files:
            if file[-4:] != ".dat" and file[-4:] != ".out":
                print(file, "Error! File is not DAT or OUT format.\n"
                      "Exiting now...")
                sys.exit()

            # iterate through files and make sure their formats are okay even
            # though it is requiered that they are in emf-inputgen.py format
            with open(file) as fo:
                for line in fo:
                    line = line.strip().split()

                    # assuming that file is already in the correct format
                    # ignore the first few lines that have length of 1
                    if line[0][0] == "#" or len(line) == 1:
                        continue

                    # check that the rest of lines have the correct format
                    else:
                        # NOE, R1, R2 data files
                        if len(line) > 1 and len(line) < 4:
                            print(file, "\nError: Data file has incorrect"
                                  "format.\nExiting now...\n")

                        # eMF digested oputput data files
                        if len(line) > 4 and len(line) < 16 or len(line) > 16:
                            print(file, "\nError: Data file has incorrect"
                                  "format.\nExiting now...\n")
                            sys.exit()
                        if len(line) == 16:
                            emf = True

    return files, emf


# get spins available in both files and relaxation data for all spins
def get_data(files, emf):

    # initialize sets and lists
    res_a = set()
    res_b = set()
    data_a = []
    data_b = []

    i = 0
    for file in files:
        with open(file) as fo:
            for line in fo:
                line = line.strip().split()
                if len(line) > 1 and line[0][0] != '#':

                    if emf is False:
                        buff = line[0:2]
                    else:
                        buff = str(line[0]) + ',' + str(line[1]) + ',' + \
                               str(line[7]) + ',' + str(line[9])
                        buff = buff.split(',')

                    if i == 0:
                        data_a.append(buff)
                        res_a.add(int(buff[0]))
                    else:
                        data_b.append(buff)
                        res_b.add(int(buff[0]))
        i += 1

    return (res_a & res_b), data_a, data_b


# remove data for residues that are not present in both data sets
def filter_data(data):

    res = data[0]
    data_a = data[1]
    data_b = data[2]
    data_buff = []

    for i in range(0, len(data_a)):
        if int(data_a[i][0]) in res:
            data_buff.append(data_a[i])
        if i == len(data_a)-1:
            data_a.clear()
            data_a = data_buff[:]
    data_buff.clear()

    for i in range(0, len(data_b)):
        if int(data_b[i][0]) in res:
            data_buff.append(data_b[i][1:])
        if i == len(data_b)-1:
            data_b.clear()
            data_b = data_buff[:]
    data_buff.clear()

    return data_a, data_b


# write output data file
def write_data_file(files, emf, data_a, data_b):
    data_file = "data4comp.out"
    fo = open(data_file, 'w')

    # header
    tab = "    "
    buff = "# Data files: " + str(files[0]) + tab + str(files[1]) + '\n'

    # body for eMF digested files
    if emf is True:
        buff += "# Res    S2-1     te-1     Rex-1     S2-2    te-2    Rex-2\n"
        fo.write(buff)
        for i in range(0, len(data_a)):
            buff = str(data_a[i][0]) + tab + tab + str(data_a[i][1]) + tab + \
                   str(data_a[i][2]) + tab + str(data_a[i][3]) + tab + \
                   str(data_b[i][0]) + tab + str(data_b[i][1]) + tab + \
                   str(data_b[i][2]) + "\n"
            fo.write(buff)
        fo.close()

    # body for NOE, R1, R2 files
    else:
        buff += "# Res    val-1    val-2\n"
        fo.write(buff)
        for i in range(0, len(data_a)):
            buff = str(data_a[i][0]) + tab + tab + str(data_a[i][1]) + tab + \
                   str(data_b[i][0]) + "\n"
            fo.write(buff)
        fo.close()


# MAIN PROGRAM
def main(dfiles):
    data = get_data(dfiles[0], dfiles[1])
    data = filter_data(data)
    write_data_file(dfiles[0], dfiles[1], data[0], data[1])


# RUN PROGRAM
main(get_dat_files())

