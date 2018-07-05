import re
import datetime
import argparse
import yaml

# Parse the sample name from a fastq file uploaded to s3.
def parseReadFile(filepath):
	path = filepath.split('/')
	filename = path[-1]
	sample = filename.split('.')
	return sample[0]

parser = argparse.ArgumentParser()
parser.add_argument('read')
parser.add_argument('run_name')
parser.add_argument('author')
parser.add_argument('version')
parser.add_argument('description')
parser.add_argument('out_file')
args = parser.parse_args()

data = dict(
	title = "pRESTO Report: " + args.run_name,
	author = args.author,
	version = args.version,
	description = args.description,
	sample = parseReadFile(args.read),
	run = args.run_name,
	date = datetime.datetime.today().strftime('%Y-%m-%d')
)

with open (args.out_file, 'w') as outfile:
	yaml.dump(data, outfile, default_flow_style=False)