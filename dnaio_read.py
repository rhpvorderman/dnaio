#!/usr/bin/env python3

import dnaio
import sys
import xopen

import functools

opener = functools.partial(xopen.xopen, threads=0)

if __name__ == "__main__":
    with dnaio.open(sys.argv[1], opener=opener) as records:
        for record in records:
            pass
