#!/usr/bin/env python
import sys
try:
    f, s = sys.argv[1:]
except:
    sys.exit('Expecting two files: the lines of each will be compared.')
fc = set([i.strip() for i in open(f, 'rU').readlines()])
sc = set([i.strip() for i in open(s, 'rU').readlines()])
i = list(fc.intersection(sc))
i.sort()
sys.stdout.write('\n'.join(i))
sys.stdout.write('\n')