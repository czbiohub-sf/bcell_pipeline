# Gerry Meixiong's script to parse the distance threshold from shazam tuning output.

import csv
import sys

tsv = open(sys.argv[1], 'r')
values = list(csv.reader(tsv, dialect="excel-tab"))
threshold = [y for x,y in values if x == 'threshold']
print(float(threshold[0]))