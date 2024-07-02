"""
Script to automatically create a manifest file for Qiime2
"""
import glob
import os
import pathlib
import argparse


# Parse script inputs
parser = argparse.ArgumentParser(description='Generate manifest file for QIIME2')
parser.add_argument('-d', required=False, type=str, default='./', help='Path to directory containing the .fastq read files')
parser.add_argument('-o', required=False, type=str, default='./QIIME2_manifest.txt', help='Output manifest file')
parser.add_argument('--frtag', required=False, type=str, default='*R1.fastq', help='Tag')
parser.add_argument('--rrtag', required=False, type=str, default='*R2.fastq', help='Tag')
args = parser.parse_args()


fwd_seq_files = glob.glob(args.d + '/' + args.frtag)
rev_seq_files = glob.glob(args.d + '/' + args.rrtag)

abs_path = os.getcwd()

fwd_seq_files.sort()
rev_seq_files.sort()

smp = []
for i in fwd_seq_files:
  iname = pathlib.Path(i).name[:-len(args.frtag)]
  smp.append(iname)

fwd_dict = dict(zip(smp,fwd_seq_files))
rev_dict = dict(zip(smp,rev_seq_files))

header = 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath'

manifest_list = []
manifest_list.append(header)

for id in smp:
  tmp = []
  tmp.append(id)
  fwd_fl = fwd_dict[id]
  #fwd_pth = abs_path + '/' + fwd_fl
  fwd_pth = fwd_fl
  tmp.append(fwd_pth)
  rev_fl = rev_dict[id]
  # rev_pth = abs_path + '/' + rev_fl
  rev_pth = rev_fl
  tmp.append(rev_pth)
  tmpj = "\t".join(tmp)
  manifest_list.append(tmpj)

manifest_j = "\n".join(manifest_list)
manifest_file = open(args.o, "w")
manifest_file.write(manifest_j)
manifest_file.close()
