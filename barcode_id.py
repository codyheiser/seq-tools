# -*- coding: utf-8 -*-
"""
barcode_id.py
@author: Cody Heiser

usage: barcode_id.py [-h] [--id-pattern ID_PATTERN] [--bc-pattern BC_PATTERN]
                     [--n-id N_ID] [--n-bc N_BC]
                     R1_file R2_file outfile

Trim a file between two strings

positional arguments:
  R1_file               input .txt file with read 1 (perform rev. comp.)
  R2_file               input .txt file with read 2
  outfile               name of output file

optional arguments:
  -h, --help            show this help message and exit
  --id-pattern ID_PATTERN
                        string to match following n bp ID
  --bc-pattern BC_PATTERN
                        string to match following n bp ID
  --n-id N_ID           number of bp to return before pattern
  --n-bc N_BC           number of bp to trim barcode sequence to
"""
import argparse
import sys
import re
import pandas as pd


def reverse_comp(sequence, nucleic_acid):
    """
    Generate reverse compliment to given sequence, s,
    in desired nucleic acid format (i.e. DNA or RNA).
    Output is reported 3' to 5'.

    Parameters:
        sequence (str): sequence to test.
        nucleic_acid (str): "DNA" or "RNA" defining desired output
    """
    # define base complements
    bases = {
        "A": {"DNA": "T", "RNA": "U"},
        "T": {"DNA": "A", "RNA": "A"},
        "G": {"DNA": "C", "RNA": "C"},
        "C": {"DNA": "G", "RNA": "G"},
        "N": {"DNA": "N", "RNA": "N"},
        "a": {"DNA": "T", "RNA": "U"},
        "t": {"DNA": "A", "RNA": "A"},
        "g": {"DNA": "C", "RNA": "C"},
        "c": {"DNA": "G", "RNA": "G"},
        "n": {"DNA": "n", "RNA": "n"},
        " ": {"DNA": " ", "RNA": " "},
        "\t": {"DNA": "\t", "RNA": "\t"},
        "\n": {"DNA": "\n", "RNA": "\n"},
        ",": {"DNA": ",", "RNA": ","},
    }
    # return reverse complimentary string (3' - 5')
    return "".join([bases[x][nucleic_acid] for x in sequence])[::-1]


def locate_id(seq, id_pattern="AAACAAAC", n_id=6):
    """
    get n_id bp before chosen pattern to return as barcode ID

    Parameters:
        seq (str): sequence to test.
        id_pattern (regex): sequence to find in seq.
            n_id characters before this pattern will be returned.
        n_id (int): number of bp to take before pattern as ID.

    Returns:
        id (str): n bp before pattern detected in seq.
    """
    i = [m.start(0) for m in re.finditer(id_pattern, seq)]
    if not i:
        return ""
    else:
        return seq[i[0] - n_id : i[0]]


def locate_bc(seq, bc_pattern="GGGGGT", n_bc=30):
    """
    get n_bc bp after chosen pattern to return as barcode

    Parameters:
        seq (str): sequence to test.
        pattern (regex): sequence to find in seq.
            n_bc characters after this pattern will be returned.
        n_bc (int): number of bp to take after pattern as barcode.

    Returns:
        id (str): n_bc bp after pattern detected in seq.
    """
    i = [m.start(0) for m in re.finditer(bc_pattern, seq)]
    if not i:
        return ""
    else:
        return seq[i[0] + 3 : i[0] + 3 + n_bc]


def compile_id(file, rev_comp=False, id_kwargs={}, bc_kwargs={}):
    """
    perform barcode compilation from all lines in file

    Parameters:
        file (str): file with sequences to test.
        rev_comp (bool): take reverse compliment of each sequence first?
        id_kwargs: to pass to locate_id()
        bc_kwargs: to pass to locate_bc()

    Returns:
        ids (pd.DataFrame): df containing output. "seq" column containing
            first n_seq bp of sequence, "id" column containing IDs.
    """
    df = pd.read_csv(file, header=None)
    if rev_comp:
        df[0] = df[0].apply(reverse_comp, nucleic_acid="DNA")
    df["id"] = df[0].apply(locate_id, **id_kwargs)
    df["barcode"] = df[0].apply(locate_bc, **bc_kwargs)
    df.rename(columns={0: "seq"}, inplace=True)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim a file between two strings")
    parser.add_argument(
        "R1_file", help="input .txt file with read 1 (perform rev. comp.)"
    )
    parser.add_argument("R2_file", help="input .txt file with read 2")
    parser.add_argument("outfile", help="name of output file")
    parser.add_argument(
        "--id-pattern",
        help="string to match following n bp ID",
        type=str,
        required=False,
        default="AAACAAAC",
    )
    parser.add_argument(
        "--bc-pattern",
        help="string to match following n bp ID",
        type=str,
        required=False,
        default="GGGGGT",
    )
    parser.add_argument(
        "--n-id",
        help="number of bp to return before pattern",
        type=int,
        required=False,
        default=6,
    )
    parser.add_argument(
        "--n-bc",
        help="number of bp to trim barcode sequence to",
        type=int,
        required=False,
        default=30,
    )
    args = parser.parse_args()

    out1 = compile_id(
        args.R1_file,
        id_kwargs={"n_id": args.n_id, "id_pattern": args.id_pattern},
        bc_kwargs={"n_bc": args.n_bc, "bc_pattern": args.bc_pattern},
        rev_comp=True,
    )
    out2 = compile_id(
        args.R2_file,
        id_kwargs={"n_id": args.n_id, "id_pattern": args.id_pattern},
        bc_kwargs={"n_bc": args.n_bc, "bc_pattern": args.bc_pattern},
    )
    out = pd.concat([out1, out2], join="outer", axis=0)
    out.to_csv(args.outfile, sep=",", index=False)
