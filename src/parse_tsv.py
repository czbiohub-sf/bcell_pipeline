# Gerry Meixiong's script to parse the distance threshold from shazam tuning output.

import sys

import pandas as pd

tsv = pd.read_csv(sys.argv[1], delimiter="\t")
print(tsv["threshold"][0])
