#! /usr/bin/env python3

import argparse
import subprocess
import sys
import os

script_dir=os.path.dirname(os.path.realpath(__file__))

def cmd(command,stdout=sys.stderr,stderr=sys.stderr):
	print(file=sys.stderr)
	print(command,file=sys.stderr)
	print(file=sys.stderr)
	subprocess.Popen(command, shell=True, universal_newlines=True, stdout=stdout, stderr=stderr).wait()
	print(file=sys.stderr)
	print(file=sys.stderr)


parser = argparse.ArgumentParser(description='One command exk matching.')
parser.add_argument(
		'-k',
		type=int,
		metavar='int',
		dest='k',
		required=True,
		help='k-mer length',
	)
parser.add_argument(
		'in_fasta',
		type=str,
		help='Input FASTA reference.',
	)
parser.add_argument(
		'in_fq',
		type=str,
		help='Reads to be matched.',
	)

args = parser.parse_args()

bwa=os.path.join(script_dir,"bwa")
exk=os.path.join(script_dir,"exk")

cmd('"{bwa}" index "{fa}"'.format(bwa=bwa,fa=args.in_fasta))
cmd('"{exk}" index -k {k} "{fa}"'.format(exk=exk,fa=args.in_fasta,k=args.k))
cmd('"{exk}" match -v -k {k} "{fa}" "{fq}"'.format(exk=exk,fa=args.in_fasta,fq=args.in_fq,k=args.k),stdout=sys.stdout)
