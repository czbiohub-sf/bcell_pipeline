## Immcantation B-cell Repertoire Sequencing Pipeline adapted for Reflow. 

Author: Eric Waltari, Gerry Meixiong, Aaron McGeever, CZ Biohub

Date: 2.22.19

## Contents ##
1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Usage](#usage)
4. [Questions](#questions)

## Introduction

The Immcantation framework is a start-to-finish pipeline going from raw reads 
to repertoire analysis. It encompasses pre-processing of fastq files, pRESTO, 
Change-O Igblast, TigGER, SHazaM, Change-O Clone, and Alakazam. All of these steps
up to Alakazam summaries are included in bcell.rf. The reflow
pipeline currently uses the 2.3 version of kleinstein/immcantation.



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

It is important to note that the read1_file and read2_file inputs do not necessarily correspond 
to your R1 and R2 fastq files from the sequencer. The read1_file input must be the fastq sequence 
beginning with the C-region or J-segment. The read2_file input much be the matching fastq sequence 
beginning with the leader V-segment. read1_primers and read2_primers enumerate the primer sequences 
for read1_file and read2_file, respectively. 
Results will be saved to the s3 bucket under the run_name directory (i.e. Alakazam results will be 
under the run_name/alakazam/ directory). The -cache off option prevents caching and forces a full run. 



## Questions

For questions, please contact eric.waltari@czbiohub.org or aaron.mcgeever@czbiohub.org. 
