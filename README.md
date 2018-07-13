## Immcantation B-cell Repertoire Sequencing Pipeline adapted for Reflow. 

##Author: Gerry Meixiong, CZ Biohub
##Date: 7.13.18

## Contents ##
1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Usage](#usage)
4. [Questions](#questions)

## Introduction

The Immcantation framework is a start-to-finish pipeline going from raw reads 
to repertoire analysis. It encompasses pre-processing of fastq files, pRESTO, 
changeo-igblast, tigger, shazam, changeo-clone, and alakazam. The steps leading
up to analysis are included in bcell.rf. The resulting output from changeo-clone,
the end stage of bcell.rf, is used as input for alakazam analysis steps. The reflow
pipeline currently uses the 1.10.2 version of kleinstein/immcantation.



## Setup

To run locally, install reflow from https://github.com/grailbio/reflow. Reflow is
implemented in Go, so you want to install that as well. 

Make sure you have AWS credentials saved to your environment by running aws configure or
exporting them like this:
export AWS_ACCESS_KEY_ID=[your key id]
export AWS_SECRET_ACCESS_KEY=[your secret key]

Then run the following: 
AWS_SDK_LOAD_CONFIG=1 reflow setup-ec2
AWS_SDK_LOAD_CONFIG=1 reflow setup-s3-repository czbiohub-reflow-quickstart-cache
AWS_SDK_LOAD_CONFIG=1 reflow setup-dynamodb-assoc czbiohub-reflow-quickstart

Your reflow config should now be configured with the biohub's s3 bucket and ready to run.


To run on an ec2 instance, use aegea to launch the czbiohub-reflow packer image.
aegea launch --iam-role S3fromEC2 --ami-tags Name=czbiohub-reflow -t t2.micro  [yourname]-reflow
aegea start [yourinstance]

Copy your read fastqs and primers fastas to the reflow/ directory in this repo. Then copy the entire repo to the launched instance using scp or aegea scp as the reflow/ directory. ssh into the instance. You should now be able to run the pipeline. Be sure to aegea stop your instance when you are done using reflow.



## Usage

For single runs (two read files with corresponding primers), cd to the reflow/ directory. 

reflow [-cache off] run bcell.rf -read1_file <read1.fastq> -read2_file <read2.fastq> \
-read1_primers <r1_primers.fasta> -read2_primers <r2_primers.fasta> -run_name <name> \
-results_bucket s3://<bucket> 

Results will be saved to the s3 bucket under the run_name directory. The -cache off option 
prevents caching and forces a full run. 


For batch runs, cd to the experiments/bcell_initial_runs/ directory. Modify bcell_test_run.csv
to match your input files. 

reflow runbatch [--reset]

The reset option will retry all the failed runs, if there were any previously. 


It is important to note that the read1_file and read2_file inputs do not necessarily correspond 
to your R1 and R2 fastq files from the sequencer. The read1_file input must be the fastq sequence 
beginning with the C-region or J-segment. The read2_file input much be the matching fastq sequence 
beginning with the leader V-segment. read1_primers and read2_primers enumerate the primer sequences 
for read1_file and read2_file, respectively. 


For analysis steps, save the germ-pass tsv file in changeo-clone/ in s3. 

reflow run alakazam.rf -changeo_file <germ-pass.tab> -run_name <name> -results_bucket s3://<bucket>

Results will be saved to the s3 bucket under the run_name/alakazam/ directory. 


## Questions

For questions, please contact gerry.meixiong@czbiohub.org or aaron.mcgeever@czbiohub.org. 