# -*- coding: utf-8 -*-
"""
trim_txt.py
@author: Cody Heiser

usage: trim_txt.py [-h] file start end

Trim a file between two strings

positional arguments:
  file        input .txt file
  start       string to start grabbing text after
  end         string to finish grabbing text before

optional arguments:
  -h, --help  show this help message and exit
"""
import argparse
import sys
import re


def trim_file(filepath, startstr, endstr):
    '''read .txt file from path and print everything between startstr and endstr'''
    F=open(filepath)
    text=F.read()
    reg=re.compile(r'{}(.*){}'.format(startstr, endstr),re.DOTALL)

    for match in reg.finditer(text):
        print (match.groups()[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim a file between two strings')
    parser.add_argument('file', help='input .txt file')
    parser.add_argument('start', help='string to start grabbing text after')
    parser.add_argument('end', help='string to finish grabbing text before')
    args = parser.parse_args()

    trim_file(args.file, args.start, args.end)
