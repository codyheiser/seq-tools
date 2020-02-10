# -*- coding: utf-8 -*-
"""
barcode_id.py
@author: Cody Heiser

usage: barcode_id.py [-h] [-p PATTERN] [--n-id N_ID] [--n-seq N_SEQ]
                     file outfile

Trim a file between two strings

positional arguments:
  file                  input .txt file
  outfile               name of output file

optional arguments:
  -h, --help            show this help message and exit
  -p PATTERN, --pattern PATTERN
                        string to match following n bp ID
  --n-id N_ID           number of bp to return before pattern
  --n-seq N_SEQ         number of bp to trim sequence to
"""
import argparse
import sys
import re
import pandas as pd


def locate_id(seq, pattern="AAACAAAC", n_id=6):
    """
    get n_id bp before chosen pattern to return as barcode ID

    Parameters:
        seq (str): sequence to test.
        pattern (regex): sequence to find in seq.
            n_id characters before this pattern will be returned.
        n_id (int): number of bp to take before pattern as ID.

    Returns:
        id (str): n bp before pattern detected in seq.
    """
    i = [m.start(0) for m in re.finditer(pattern, seq)]
    if not i:
        return ""
    else:
        return seq[i[0]-n_id : i[0]]


def compile_id(file, n_seq=40, **kwargs):
    """
    perform barcode compilation from all lines in file

    Parameters:
        file (str): file with sequences to test.
        n_seq (int): number of bp to trim sequence to (from front of string).
            use -1 to keep entire sequence
        **kwargs: to pass to locate_id()

    Returns:
        ids (pd.DataFrame): df containing output. "seq" column containing
            first n_seq bp of sequence, "id" column containing IDs.
    """
    df = pd.read_csv(file, header=None)
    df["id"] = df[0].apply(locate_id, **kwargs)
    df.rename(columns={0:"seq"}, inplace=True)
    df["seq"] = df["seq"].apply(lambda x: x[:n_seq])
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim a file between two strings')
    parser.add_argument('file', help='input .txt file')
    parser.add_argument('outfile', help='name of output file')
    parser.add_argument('-p','--pattern', help='string to match following n bp ID', required=False, default="AAACAAAC")
    parser.add_argument('--n-id', help='number of bp to return before pattern', required=False, default=6)
    parser.add_argument('--n-seq', help='number of bp to trim sequence to', required=False, default=40)
    args = parser.parse_args()

    out = compile_id(args.file, n_seq=args.n_seq, n_id=args.n_id, pattern=args.pattern)
    out.to_csv(args.outfile, sep=",", index=False)
