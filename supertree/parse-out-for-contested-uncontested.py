#!/usr/bin/env python
import sys
import re
try:
  inpfn, contestedfn, uncontestedfn = sys.argv[1:]
except:
  sys.exit('''Expecting 3 arguments:
  1. the filename of the stderr from the otc-decompose-uncontested
      run with the INFO level debugging enabled.
  2. a filepath for the list of contested IDs
  3. a filepath for a list of the uncontested internal node OTT Ids
      that were encountered.
''')
try:
  lineit = iter(open(inpfn, 'rU'))
except:
  sys.exit('Could not open "{}".'.format(inpfn))

prefix = r'^[-0-9:, ]+ INFO +\[default\] +'
check_pat = re.compile(prefix + r'exportOrCollapse for ott([0-9a-zA-Z]+)\s.*$')
result_pat = re.compile(prefix + r'(\S+ontested)\s*$')
contested, uncontested = [], []
while True:
  try:
    line = lineit.next()
  except:
    break
  m = check_pat.match(line)
  if m:
    ott_id = m.group(1)
    n = lineit.next()
    m = result_pat.match(n)
    if not m:
      sys.exit('Expecting Contested/Uncontested line. Got:\n{}\n'.format(n))
    r = m.group(1).lower()
    if r == 'contested':
      contested.append(ott_id)
    else:
      if r != 'uncontested':
        sys.exit('Expecting Contested/Uncontested line. Got:\n{}\n'.format(n))
      uncontested.append(ott_id)

def write_id_list(c, cfn):
  with open(cfn, 'w') as out:
    out.write('\n'.join(c))
    out.write('\n')

write_id_list(contested, contestedfn)
write_id_list(uncontested, uncontestedfn)
