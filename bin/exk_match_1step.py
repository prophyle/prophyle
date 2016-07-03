#! /usr/bin/env python3

import argparse
import subprocess
import sys
import os

script_dir=os.path.dirname(os.path.realpath(__file__))
bwa=os.path.join(script_dir,"bwa")
exk=os.path.join(script_dir,"exk")

def cmd(command,stdout=sys.stderr,stderr=sys.stderr):
	print(file=sys.stderr)
	print("[exk_match_1step.py]",command,file=sys.stderr)
	#print(file=sys.stderr)
	subprocess.Popen(command, shell=True, universal_newlines=True, stdout=stdout, stderr=stderr).wait()
	#print(file=sys.stderr)
	#print(file=sys.stderr)

def create_bwa_index(fa):
	cmd('"{bwa}" index "{fa}"'.format(bwa=bwa,fa=fa))

def create_klcp(fa, k):
	cmd('"{exk}" index -k {k} "{fa}"'.format(exk=exk,fa=fa,k=k))

def match(fa, fq, k, s=False, u=False):
	params=""
	if s:
		params+=" -s"
	if u:
		params+=" -u"
	cmd('"{exk}" match {params} -v -k {k} "{fa}" "{fq}"'.format(exk=exk,fa=fa,fq=fq,k=k, params=params),stdout=sys.stdout)

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
		'-u',
		action='store_true',
		help='use rolling window',
	)
parser.add_argument(
		'-s',
		action='store_true',
		help='skip k-1 k-mers after failing matching k-mer',
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

fa=args.in_fasta
fq=args.in_fq
k=args.k
u=args.u
s=args.s

create_bwa_index(fa)

if u:
	create_klcp(fa, k)
	#cmd('"{exk}" index -k {k} "{fa}"'.format(exk=exk,fa=args.in_fasta,k=args.k))

#cmd('"{exk}" match -v -k {k} "{fa}" "{fq}"'.format(exk=exk,fa=args.in_fasta,fq=args.in_fq,k=args.k),stdout=sys.stdout)
match(fa, fq, k, s=s, u=u)
