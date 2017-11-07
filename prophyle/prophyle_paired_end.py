#! /usr/bin/env python3

"""Merge paired-end FASTA or FASTQ files (possibly in gzip format)

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT
"""

import argparse
import os
import sys
import errno

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def read_id(read_1, read_2):
    # two possibilities:
    # <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>   <pair idx>:<is filtered>:<control number>:<index sequence>
    # <instrument>:<flowcell lane>:<tile>:<x-pos>:<y-pos>:<multiplex idx>/<pair idx>
    for sep in ['/', ' ']:
        try:
            id1 = sep.join(read_1.split(sep)[:-1])
            id2 = sep.join(read_2.split(sep)[:-1])
            if id1 != id2 or len(id1) < 0.3 * len(read_1):
                raise ValueError('different id')
            else:
                return id1
        except (KeyError, ValueError):
            pass
    # cannot find id, just use the one of the first read
    return read_1


def merge_reads(f_reads_1, f_reads_2, out_file):
    try:
        reads_1 = pro.open_gzip(f_reads_1)
        reads_2 = pro.open_gzip(f_reads_2)

        first_read_1 = reads_1.readline()
        first_read_2 = reads_2.readline()

        first_char_1 = first_read_1[0]
        first_char_2 = first_read_2[0]
        assert first_char_1 == first_char_2, 'paired-end files of different format'
        if first_char_1 == '@':
            fastq, fasta = True, False
        elif first_char_1 == '>':
            fastq, fasta = False, True
        else:
            print('[prophyle_paired_reads] Error: unknown read format (does not start with > or @)',
                file=sys.stderr)
            reads_1.close()
            reads_2.close()
            sys.exit(1)

        rid=read_id(first_read_1.strip(), first_read_2.strip())
        out_file.write(rid+'\n')

        if fastq:
            # file line
            i = 0
            for next_line_1, next_line_2 in zip(reads_1, reads_2):
                i += 1
                l = i % 4
                if l == 0:
                    assert next_line_1.startswith('@') and next_line_2.startswith('@'), \
                        "malformed fastq files (no id at line {})".format(i)
                    rid=read_id(next_line_1.strip(), next_line_2.strip())
                    out_file.write(rid+'\n')
                elif l == 1:
                    out_file.write(next_line_1.strip() + 'NNN' + next_line_2.strip() + '\n')
                elif l == 2:
                    out_file.write('+\n')
                elif l == 3:
                    out_file.write(next_line_1.strip() + '!!!' + next_line_2.strip() + '\n')

        elif fasta:
            prev_read_1 = ''
            prev_read_2 = ''
            for next_line_1, next_line_2 in zip(reads_1, reads_2):
                if next_line_1.startswith('>'):
                    out_file.write(prev_read_1 + 'NNN' + prev_read_2 + '\n')
                    prev_read_1 = ''
                    prev_read_2 = ''
                    rid=read_id(next_line_1.strip(), next_line_2.strip())
                    out_file.write(rid+'\n')
                else:
                    prev_read_1 += next_line_1.strip()
                    prev_read_2 += next_line_2.strip()
            out_file.write(prev_read_1 + 'NNN' + prev_read_2 + '\n')

        out_file.flush()
        out_file.close()
        if reads_1.readline().strip() != '' or reads_2.readline().strip() != '':
            print('[prophyle_paired_reads] Warning: files of different length (merged till the end of the shortest)',
                file=sys.stderr)

    except IOError as e:
        if e.errno == errno.EPIPE:
            print('[prophyle_paired_reads]'+reads_1.readline().strip()+reads_2.readline().strip(),file=sys.stderr)
            if reads_1.readline().strip() == '' and reads_2.readline().strip() == '':
                # pipe closed but both files correctly processed
                pass
            else:
                print('[prophyle_paired_reads] Error: pipe closed before EOF',
                    file=sys.stderr)
                reads_1.close()
                reads_2.close()
                sys.exit(1)
        else:
            reads_1.close()
            reads_2.close()
            raise

    reads_1.close()
    reads_2.close()


def main():
    desc = """\
        Program: prophyle_paired_end.py

        Merge paired-end FASTA or FASTQ files (possibly in gzip format).
    """

    parser = argparse.ArgumentParser(
        description=desc
    )

    parser.add_argument(
        'reads_1',
        type=str,
        help='1st FASTA or FASTQ file',
    )

    parser.add_argument(
        'reads_2',
        type=str,
        help='2nd FASTA or FASTQ file'
    )

    parser.add_argument('-o', '--output-file',
        type=argparse.FileType('w'),
        metavar='out_file',
        dest='out_file',
        help='output file [stdout]',
        default=sys.stdout
    )

    args = parser.parse_args()
    merge_reads(args.reads_1, args.reads_2, args.out_file)


if __name__ == '__main__':
    main()
