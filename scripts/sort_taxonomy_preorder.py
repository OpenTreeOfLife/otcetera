#!/usr/bin/env python
import sys

def main(taxonomy_fp, root_id, out):
    by_par_id = {}
    newly_emitted = set()
    emitted = set()
    with open(taxonomy_fp, 'r') as inp:
        liter = iter(inp)
        header = next(liter)
        out.write(header)
        for line in inp:
            if not line:
                continue
            ls = line.split('\t')
            line_id = ls[0].strip()
            if line_id == root_id:
                out.write(line)
                newly_emitted.add(root_id)
            else:
                assert ls[1] == '|'
                par = ls[2].strip()
                by_par_id.setdefault(par, []).append((line_id, line))
    while newly_emitted:
        last_emitted = set(newly_emitted)
        emitted.update(newly_emitted)
        newly_emitted = set()
        leml = list(last_emitted)
        leml.sort()
        for par_id in leml:
            kids = by_par_id.get(par_id)
            if not kids:
                continue
            kids.sort()
            for kid in kids:
                kid_id, kid_line = kid
                if kid_id in emitted:
                    raise RuntimeError(f"{kid_id} emitted > 1 time. second time from {par_id}")
                newly_emitted.add(kid_id)
                emitted.add(kid_id)
                out.write(kid_line)
    sys.stderr.write(f"{len(emitted)} records emitted\n")



if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.stdout)