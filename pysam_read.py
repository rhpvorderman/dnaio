import pysam
import sys

with pysam.FastxFile(sys.argv[1]) as f:
    for record in f:
        pass
