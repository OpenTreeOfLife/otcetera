#!/usr/bin/env python
import sys
try:
    from peyotl.supertree import OtcPipelineSimplifiedDump
except ImportError:
    sys.exit('Depends on some "bleeding edge" feature on the peyotl supertree branch')
import os

SCRIPT_NAME = os.path.split(sys.argv[0])[1]
inp, out = sys.argv[1:]
d = OtcPipelineSimplifiedDump(inp)
d.solve(out)