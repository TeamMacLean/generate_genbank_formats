#!/usr/bin/env python

import sys, os, re

fh=open(sys.argv[1])
basecount=sys.argv[2]

firstline=fh.readline().rstrip().replace("DNA", "DNA\tcircular")
print(firstline)
for line in fh:
    line=line.rstrip()
    if line.startswith("ORIGIN"):
        line=line.replace("ORIGIN", "BASE COUNT\t" + basecount + "\nORIGIN")
        print(line)
    else:
        print(line)

fh.close()
