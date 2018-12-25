#!/usr/bin/env python

"""
Simple script for rename the fastq.gz files to specific format.

usage: python rename.py DATADIR

NOTE: All fastq file must be gziped
"""

import os
from os.path import join
import re
import sys

if sys.version_info < (3, 0, 0):
    input = raw_input

def print_usage():
    print(__doc__)

name_patterns = [
    "(.*?)_(1|2).fq.gz",
    "(.*?)_R(1|2).fq.gz",
    "(.*?)_(1|2).fastq.gz",
]

target_pattern = "(.*?)_R(1|2).fastq.gz"

target_format  = "{id}_R{pid}.fastq.gz"

def new_name(id, pid):
    target_name = target_format.format(id=id, pid=pid)
    return target_name

def match(fname):
    for p in name_patterns:
        match = re.match(p, fname)
        if match:
            id, pid = match.groups()
            n = new_name(id, pid)
            print("{} -> {}".format(fname, n))
            return fname, n
    else:
        match = re.match(target_pattern, fname)
        if match:
            print("correct partten: {}".format(fname))
            return fname
        else:
            return None


if __name__ == "__main__":
    argv = sys.argv

    if len(argv) != 2:
        print_usage()
        sys.exit(1)

    data_dir = argv[1]

    fname_list = os.listdir(data_dir)

    rename_lst = []
    for fname in fname_list:
        m = match(fname)
        if isinstance(m, tuple):
            rename_lst.append(m)

    if rename_lst:
        reply = input("Are you sure rename them? ([no]/yes) ")
        if reply == 'yes' or reply == 'y':
            for old, new in rename_lst:
                old, new = join(data_dir, old), join(data_dir, new)
                os.rename(old, new)
