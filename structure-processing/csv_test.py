#/usr/bin/env python2
import csv
import sys

with open('1A2K_TP0_backrub_1.sc', 'r') as inf:
    # #next(inf)
    # dialect = csv.Sniffer().sniff(inf.read(), delimiters=' ')
    # inf.seek(0)
    # next(inf)
    # reader = csv.reader(inf, dialect)
    # for row in reader:
    #     print row
    lines = inf.readlines()
    lines = [i.strip('\n').split() for i in lines]
    for line in lines[1:]:
        print line[1], line[-1]
