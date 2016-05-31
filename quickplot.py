#!/usr/bin/python

from . import helplib as hl
import matplotlib

import matplotlib.pyplot as plt

from . import fitlib as fl

import argparse
import os

#using the argparse module to make use of command line options
parser = argparse.ArgumentParser(description="Quickly plot an energy scan")

#we need at least one thing to fit
filegroup = parser.add_mutually_exclusive_group(required=True)
filegroup.add_argument("--folder", help="Specify a folder. All files contained in the folder will be fitted")
filegroup.add_argument("--filelist", help="Specify a file that includes a list of filenames to be fitted")
filegroup.add_argument("--filename", help="Specify a filename to be fitted")

parser.add_argument('--version', action='version', version='r1')

#parse it
args = parser.parse_args()

#we need this to make a filelist list, that contains filenames
filelist = []

#here we go through the 3 cases: filename, folder, filelist

#this is the easiest case - add the filename to the filelist
if args.filename is not None:
    filelist.append([args.filename, os.path.basename(args.filename)])

elif args.folder is not None:

    #adjust file path in case we're running on fucking windows
    args.folder = os.path.normcase(args.folder)

    #lets go through that folder and add every filename to the filelist
    caldirlist = os.listdir(args.folder)

    for file in caldirlist:
        filelist.append([os.path.join(args.folder, file), file])

elif args.filelist is not None:

    #in this we have to read the list of filenames in a file
    f = hl.openfile(args.filelist)
    
    for line in f:
        #we split the array by whitespaces or tabs (split tries both)
        line = line.strip('\r\n').split()
        #first argument should be the filename + path, therefore appending it to a temp array together with the filename
        filelist.append(line[0], os.path.basename(line[0]))

#let's walk our filelist
for file in filelist:
    #this variable is set to false if we encounter a non-readable file
    usefulfile = True

    try:
        data = hl.readfile(file[0])
    except IOError:
        usefulfile = False
        print('Could not read: ' + file[0])

    if usefulfile is True:
        fig1 = plt.figure()
        fl.plotES(data, file[1])

plt.show()
