#! /usr/bin/python3

import sys

def listToString(List):
    string=''
    for word in List:
        string += word + ' '
    return string + '\n'
    
for filename in ['bbn1', 'bbn2', 'bbn3', 'bbn4']:
    dat = open(filename, "r")
    res = dat.readlines()
    dat.close()
    res = sys.argv[1] + ' ' + filename[-1] + ' ' + listToString(res[9].split()[1:])
    plot = open("plot.dat", "a")
    plot.write(res)
    plot.close()
    

