#! /usr/bin/env python3
"""Old RefSeq tree building script

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import argparse
import os
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro

GITDIR = os.path.basename(sys.argv[0])[-3:] == ".py"
if GITDIR:
    C_D = os.path.abspath(os.path.dirname(sys.argv[0]))
else:
    C_D = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

LIBRARIES = ['bacteria', 'viruses', 'plasmids']
FTP_NCBI = 'https://ftp.ncbi.nlm.nih.gov'
RPT_PARSER = os.path.join(C_D, "prophyle_parse_rpt.py")
TREE_BUILDER = os.path.join(C_D, "prophyle_ncbi_tree.py")


def download_rpt(library, library_dir):

    if library == "all":
        for l in LIBRARIES:
            download_rpt(l, library_dir)
        return
    else:
        assert library in LIBRARIES

    d = os.path.join(library_dir, library)
    #os.makedirs(d, exist_ok=True)
    #pro.makedirs(d)
    # if it does not exist, exit, there are no fna files to add to the tree!

    if library == 'bacteria':
        cmd = [
            'cd', d, '&&', 'curl', FTP_NCBI + '/genomes/archive/old_refseq/Bacteria/all.rpt.tar.gz', '|', 'tar', 'xz'
        ]
        pro.run_safe(cmd)

    elif library == 'viruses':
        cmd = ['cd', d, '&&', 'curl', FTP_NCBI + '/genomes/Viruses/all.rpt.tar.gz', '|', 'tar', 'xz']
        pro.run_safe(cmd)

    elif library == 'plasmids':
        cmd = [
            'cd', d, '&&', 'curl', FTP_NCBI + '/genomes/archive/old_refseq/Plasmids/plasmids.all.rpt.tar.gz', '|',
            'tar', 'xz', '--strip', '5'
        ]
        pro.run_safe(cmd)

    else:
        raise ValueError('Unknown library "{}"'.format(library))


def fasta_idx(library, library_dir):

    if library == "all":
        for l in LIBRARIES:
            fasta_idx(l, library_dir)
        return
    else:
        assert library in LIBRARIES

    cmd = [
        'find',
        os.path.join(library_dir, library), '-name', '*.fna', '|'
        'parallel', '--no-notice', '--verbose', 'samtools', 'faidx', '{}'
    ]
    pro.run_safe(cmd)


def parse_rpt(library, library_dir):

    if library == "all":
        for l in LIBRARIES:
            parse_rpt(l, library_dir)
        return
    else:
        assert library in LIBRARIES

    cmd = [RPT_PARSER, os.path.join(library_dir, library), '>', library + '_taxamap.tsv']
    pro.run_safe(cmd)


def build_tree(library, library_dir):

    if library == "all":
        for l in LIBRARIES:
            build_tree(l, library_dir)
        return
    else:
        assert library in LIBRARIES

    root = "Bacteria" if library == 'plasmids' else library.title()

    cmd = [
        TREE_BUILDER, library, library_dir, library + '.nw', library + '_taxamap.tsv', '-l', library + '.log', '-u',
        root
    ]
    pro.run_safe(cmd)


def main():
    parser = argparse.ArgumentParser(description='Build trees for the old refseq')

    parser.add_argument(
        'library',
        metavar='<library>',
        choices=LIBRARIES + ['all'],
        help='genomic library {}'.format(LIBRARIES + ['all']),
    )

    parser.add_argument(
        '-l',
        metavar='DIR',
        dest='library_dir',
        type=str,
        help='directory with the library sequences',
        default=None,
        # required=True,
    )

    args = parser.parse_args()

    library_dir = os.path.abspath(args.library_dir)

    download_rpt(library=args.library, library_dir=library_dir)

    fasta_idx(library=args.library, library_dir=library_dir)

    parse_rpt(library=args.library, library_dir=library_dir)

    build_tree(library=args.library, library_dir=library_dir)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # pipe error (e.g., when head is used)
        sys.stderr.close()
        exit(0)
    except KeyboardInterrupt:
        pro.message("Error: Keyboard interrupt")
        pro.close_log()
        exit(1)
