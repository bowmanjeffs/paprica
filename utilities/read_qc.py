# -*- coding: utf-8 -*-
"""
Created on Tue Feb 02 18:07:14 2016

@author: jeff

Call as:
python read_qc.py -in [fastq or fastq.gz] -phred [desired mean phred score] -n [fraction of bases that may be 'n']

Starting with the lowest scoring bases script will convert these to 'n'
characters with a NaN score until the desired mean phred score is reached.  The
resulting sequences are written to either fasta or fasta.gz, depending on the
format of the input file.

The -n flag indicates the fraction of bases that are allowed to be 'n'.  If this
fraction is exceeded the read will not be kept.  Must be in the range of 0-1.

"""

from Bio import SeqIO

import gzip
import sys

import pandas as pd
import numpy as np

command_args = {}

## default values for testing

command_args['phred'] = 35
command_args['n'] = 0.1
command_args['in'] = 'ERR164409.fastq.gz'

## override defaults

for i,arg in enumerate(sys.argv):
    if arg.startswith('-'):
        arg = arg.strip('-')
        command_args[arg] = sys.argv[i + 1]
        
phred = int(command_args['phred'])
n = float(command_args['n'])

## Create containers for useful statistics on the QC.

kept = 0
discarded = 0
all_scores = []

## Define a function to carry out the QC

def qc(record, fasta_out, kept, discarded):
    
    scores = record.letter_annotations['phred_quality']
    seq = pd.Series(list(str(record.seq)))
    
    scores_series = pd.Series(scores)            
    q = scores_series.min()
    
    while scores_series.mean() < phred:
        scores_series[scores_series == q] = np.nan
        q = q + 1
        
    seq[pd.notnull(scores_series) == False] = 'n'                            
    new_seq = ''.join(seq)
    
    if float(len(scores_series[scores_series > 0])) / len(seq) >= 1 - n:
    
        print >> fasta_out, '>' + record.id + '\n' + new_seq
        #print record.id, len(seq), len(scores_series[scores_series > 0]), scores_series.mean()
        kept = kept + 1
        all_scores.append(scores_series.mean())
        
    else:
        discarded = discarded + 1
        
    return kept, discarded, scores
    
## Execute the function if the file ends with fastq.gz.
    
if command_args['in'].endswith('.fastq.gz'):
    
    fastq_gz = command_args['in']
    name = fastq_gz.rstrip('.fastq.gz')
    
    with gzip.open(fastq_gz, 'rb') as fastq_in, gzip.open(name + '.fasta.gz', 'wb') as fasta_out, open(name + '.qc.txt', 'w') as summary:
        for record in SeqIO.parse(fastq_in, 'fastq'):
            kept, discarded, all_scores = qc(record, fasta_out, kept, discarded)
            
    print name, 'kept=' + str(kept), 'discarded=' + str(discarded), 'mean.score=' + str(pd.Series(all_scores).mean())
    
    print >> summary, 'kept=' + str(kept)
    print >> summary, 'discarded=' + str(discarded)
    print >> summary, 'mean.score=' + str(pd.Series(all_scores).mean())
    
## Execute the function if the file ends with fastq.
            
elif command_args['in'].endswith('.fastq'):
    
    fastq = command_args['in']
    name = fastq.rstrip('.fastq')
    
    with open(fastq, 'rb') as fastq_in, open(name + '.fasta', 'w') as fasta_out, open(name + '.qc.txt', 'w') as summary:
        for record in SeqIO.parse(fastq_in, 'fastq'):
            kept, discarded, all_scores = qc(record, fasta_out, kept, discarded)
            
    print name, 'kept=' + str(kept), 'discarded=' + str(discarded), 'mean.score=' + str(pd.Series(all_scores).mean())
    
    print >> summary, 'kept=' + str(kept)
    print >> summary, 'discarded=' + str(discarded)
    print >> summary, 'mean.score=' + str(pd.Series(all_scores).mean())
    
## Return an error if the file does not end with fastq or fastq.gz.
            
else:
    print 'extension error, input file must be fastq or fastq.gz'
        