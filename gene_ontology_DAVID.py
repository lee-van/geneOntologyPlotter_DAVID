#!/bin/python2.7
#
#  Copyright (c) 2016 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

import argparse
import subprocess
import os
import sys
import datetime
import time
import re
import os.path
import tempfile
import collections
import shutil

parser = argparse.ArgumentParser(description="Take output from DAVID and plots a GO heatmap from p-value, benjamini, and foldchange")
parser.add_argument('--input_filenames','-i',action='store', dest='input_filenames', nargs='+', required=True, help='DAVID GO charts, variable length space delimited list')
parser.add_argument('--output_tags','-t',action='store', dest='output_tags', nargs='+', required=True, help='Short tags for each input file, same length as input_filenames list')
parser.add_argument('--output_dir', '-o', action='store', dest='output_dir', nargs=1, required=True, help='Output directory')
parser.add_argument('--output_prefix', '-p', action='store', dest='output_prefix', nargs=1, required=True, help='Prefix for output files')
parser.add_argument('--levels', '-l', action='store', dest='levels', nargs=1, required=True, type=str, help='Datafile containing GO levels')
parser.add_argument('--sig_threshold', '-s', action='store', dest='sig_threshold', nargs='?', default=0.05, help='threshold of significance (optional)')
parser.add_argument('--min_samples_sig', '-m', action='store', dest='min_samples_sig', nargs='?', type=int, default=1, help='minimum number of samples with signficant enrichment for GO term to be plotted (optional)')
parser.add_argument('--boundaries', '-b', action='store', dest='boundaries', nargs='?', help='levels for plotting script (optional)')
#parser.add_argument('--transform', '-tr', action='store', dest='transform', nargs='?', default='NA', help='log transform output (optional). Options: log2, neglog2, log10, neglog10, NA')
args=parser.parse_args()

#check to make sure input and output tags are of equal length
if len(args.input_filenames) != len(args.output_tags):
    sys.exit("number of output tags does not equal number of input filenames")

#define global variables
delimiter = ";"
null = "NA"
expected_line = re.compile(r'GOTERM')

#locations of called scripts
parent_dir=os.path.abspath(os.path.join(__file__, os.pardir))
plot_heatmap=parent_dir+"/"+"plot_heatmap.R" #C-script

####################################################################################################################################################################
### Parse GO annotation and create a hash of levels
print "###Parsing GO annotations..."

levels = collections.OrderedDict()
namespaces = collections.OrderedDict()
names = collections.OrderedDict()
annotated = collections.OrderedDict()
altids = collections.OrderedDict()

infn=args.levels[0]
with open(infn) as f:
    next(f) #trim header
    for line in f:
        GO_ID, namespace, name, alt, parents, ancestors, level = line.rstrip().split("\t")
        levels[GO_ID] = level
        namespaces[GO_ID] = namespace
        names[GO_ID] = name
        annotated[GO_ID] = 1
        #if alternate ids exist, replicate data for main ID into each alt ID
        if alt != "NA":
            alts = alt.split(",")
            for a in alts:
                levels[a] = level
                namespaces[a] = namespace
                names[a] = name
                annotated[a] = 1

####################################################################################################################################################################
### Run through datasets and grab associated p, benjamini, fold. Use i as an indicator variable for which dataset is being loaded
print "###Parsing DAVID output..."

pvalues = {}
folds = {}
benjaminis = {}

#run through datasets and grab associated p, benjamini, fold. Use i as an indicator variable for which dataset is being loaded
i=-1
for infn in args.input_filenames:
    i += 1
    with open(infn) as f:
        next(f) #skip header
        for line in f:
            if expected_line.search(line) is not None:
                Category, Term, Count, percent, PValue, Genes, List_Total, Pop_Hits, Pop_Total, Fold_Enrichment, Bonferroni, Benjamini, FDR = line.rstrip().split("\t")
                GO_ID, name = Term.split("~")
                if not pvalues.has_key(GO_ID):
                    pvalues[GO_ID] = ['NA'] * len(args.input_filenames)
                pvalues[GO_ID][i] = PValue
                if not benjaminis.has_key(GO_ID):
                    benjaminis[GO_ID] = ['NA'] * len(args.input_filenames)
                benjaminis[GO_ID][i]=Benjamini
                if not folds.has_key(GO_ID):
                    folds[GO_ID] = ['NA'] * len(args.input_filenames)
                folds[GO_ID][i]=Fold_Enrichment

####################################################################################################################################################################
### Initialize output directory and print summary file. Also print tmp files for plotting

#set significance threshold
threshold = float(0.05)
if args.sig_threshold:
    threshold = args.sig_threshold
min_samples = 1
if args.min_samples_sig:
    min_samples = args.min_samples_sig

#Check for output directory and make it if neccessary
print "###Writing output..."
output_folder = re.sub('\/$', '', args.output_dir[0])
if  os.path.isdir(args.output_dir[0]): #if no out dir then make one
    print "existing output folder detected, will overwrite internal files..."
else:
    subprocess.check_call(['mkdir', output_folder])

#Write master lists of GO terms (by p-value, benjamini, and fold). Also write to tmp files in which nonsig output is removed
header = 'GO_ID\tname\tnamespace\tlevel\t'+'\t'.join(args.output_tags)+'\n'

#summary files
outp=open(output_folder+'/'+args.output_prefix[0]+'.pValueSummary.txt','w')
outb=open(output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.txt','w')
outf=open(output_folder+'/'+args.output_prefix[0]+'.foldEnrichmentSummary.txt','w')
outp.write(header)
outb.write(header)
outf.write(header)

#temp files, one for p-value and benjamini, plus associated fold changes for sig terms as determined either by p or benjamini (2 x 2 files total)
tempp = output_folder+'/'+args.output_prefix[0]+'.pValueSummary.tmp'
tempb = output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.tmp'
tempfoldp = output_folder+'/'+args.output_prefix[0]+'.foldEnrichmentSummary.p.tmp'
tempfoldb = output_folder+'/'+args.output_prefix[0]+'.foldEnrichmentSummary.b.tmp'
outtempp=open(tempp,'w')
outtempb=open(tempb,'w')
outtempfoldp=open(tempfoldp,'w')
outtempfoldb=open(tempfoldb,'w')
outtempp.write(header)
outtempb.write(header)
outtempfoldp.write(header)
outtempfoldb.write(header)

for GO_ID in pvalues:
    #print out levels if GO ID is also in GO annotation file (not always the case)
    if GO_ID in annotated: 
        #full summary files
        poutput = '\t'.join(pvalues[GO_ID])
        boutput = '\t'.join(benjaminis[GO_ID])
        foutput = '\t'.join(folds[GO_ID])
        outp.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+poutput+'\n')
        outb.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+boutput+'\n')
        outf.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+foutput+'\n')
        #print tempfile output if at least min_samples are signficant. Also convert all nonsig values in remaining entries (both p and foldchange) to NA
        count=0
        for i in range(0,len(args.input_filenames)):
            if pvalues[GO_ID][i] != 'NA':
                if float(pvalues[GO_ID][i])<=threshold:
                    count +=1
        if count>=min_samples: #if at least min_samples elements in p values are <= sig threshold
            fold = list(folds[GO_ID])
            for i in range(0,len(args.input_filenames)):
                if pvalues[GO_ID][i] != 'NA':
                    if float(pvalues[GO_ID][i])>threshold:
                        pvalues[GO_ID][i] = 'NA'
                        fold[i] = 'NA'
            poutput = '\t'.join(pvalues[GO_ID])
            pfoldoutput = '\t'.join(fold)
            outtempp.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+poutput+'\n')
            outtempfoldp.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+pfoldoutput+'\n')
        count=0
        for i in range(0,len(args.input_filenames)):
            if benjaminis[GO_ID][i] != 'NA':
                if float(benjaminis[GO_ID][i])<=threshold:
                    count +=1
        if count>=min_samples: #if at least min_samples elements in benjamini adjusted p values are <= sig threshold
            fold = list(folds[GO_ID])
            for i in range(0,len(args.input_filenames)):
                if benjaminis[GO_ID][i] != 'NA':
                    if float(benjaminis[GO_ID][i])>threshold:
                       benjaminis[GO_ID][i] = 'NA'
                       fold[i] = 'NA'
            boutput = '\t'.join(benjaminis[GO_ID])
            bfoldoutput = '\t'.join(fold)
            outtempb.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+boutput+'\n')
            outtempfoldb.write(GO_ID+'\t'+names[GO_ID]+'\t'+namespaces[GO_ID]+'\t'+levels[GO_ID]+'\t'+bfoldoutput+'\n')
    else:
        print "Warning: id "+GO_ID+" not found in GO annotation, most likely obsolete";

outp.close()
outb.close()
outf.close()
outtempp.close()
outtempb.close()
outtempfoldp.close()
outtempfoldb.close()

####################################################################################################################################################################
### Print output to R for plotting

if args.boundaries:
    if min_samples > 1:
        lower, upper = args.boundaries.split(",")
        plotP = output_folder+'/'+args.output_prefix[0]+'.pValueSummary.min'+str(min_samples)+'samples.levels'+lower+'-'+upper+'.pdf'
        plotFoldP = output_folder+'/'+args.output_prefix[0]+'.pValueAssociatedFoldEnrichmentSummary.min'+str(min_samples)+'samples.levels'+lower+'-'+upper+'.pdf'
        plotB = output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.min'+str(min_samples)+'samples.levels'+lower+'-'+upper+'.pdf'
        plotFoldB = output_folder+'/'+args.output_prefix[0]+'.benjaminiAssociatedFoldEnrichmentSummary.min'+str(min_samples)+'samples.levels'+lower+'-'+upper+'.pdf'
    else:
        lower, upper = args.boundaries.split(",")
        plotP = output_folder+'/'+args.output_prefix[0]+'.pValueSummary.levels'+lower+'-'+upper+'.pdf'
        plotFoldP = output_folder+'/'+args.output_prefix[0]+'.pValueAssociatedFoldEnrichmentSummary.levels'+lower+'-'+upper+'.pdf'
        plotB = output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.levels'+lower+'-'+upper+'.pdf'
        plotFoldB = output_folder+'/'+args.output_prefix[0]+'.benjaminiAssociatedFoldEnrichmentSummary.levels'+lower+'-'+upper+'.pdf'

else:
    lower = 'NA'
    upper = 'NA'
    if min_samples > 1:
        plotP = output_folder+'/'+args.output_prefix[0]+'.pValueSummary.min'+str(min_samples)+'samples.allLevels.pdf'
        plotFoldP = output_folder+'/'+args.output_prefix[0]+'.pValueAssociatedFoldEnrichmentSummary.min'+str(min_samples)+'samples.allLevels.pdf'
        plotB = output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.min'+str(min_samples)+'samples.allLevels.pdf'
        plotFoldB = output_folder+'/'+args.output_prefix[0]+'.benjaminiAssociatedFoldEnrichmentSummary.min'+str(min_samples)+'samples.allLevels.pdf'
    else:
        plotP = output_folder+'/'+args.output_prefix[0]+'.pValueSummary.allLevels.pdf'
        plotFoldP = output_folder+'/'+args.output_prefix[0]+'.pValueAssociatedFoldEnrichmentSummary.allLevels.pdf'
        plotB = output_folder+'/'+args.output_prefix[0]+'.benjaminiSummary.allLevels.pdf'
        plotFoldB = output_folder+'/'+args.output_prefix[0]+'.benjaminiAssociatedFoldEnrichmentSummary.allLevels.pdf'

subprocess.check_call(['Rscript',plot_heatmap,tempp,upper,lower,'neglog10', plotP])
subprocess.check_call(['Rscript',plot_heatmap,tempfoldp,upper,lower,'NA',plotFoldP])
subprocess.check_call(['Rscript',plot_heatmap,tempb,upper,lower,'neglog10',plotB])
subprocess.check_call(['Rscript',plot_heatmap,tempfoldb,upper,lower,'NA',plotFoldB])

os.remove(tempp)
os.remove(tempfoldp)
os.remove(tempb)
os.remove(tempfoldb)



