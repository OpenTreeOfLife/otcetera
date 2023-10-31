#!/usr/bin/env python3
import sys

def mapping_from_filepath(fp):
    md = {}
    with open(fp, 'r') as inp:
        for line in inp:
            ls = line.strip()
            if not ls:
                continue
            by_pipe = ls.split('|')
            assert len(by_pipe) == 3
            key = int(by_pipe[0].strip())
            val = int(by_pipe[1].strip())
            assert key not in md
            md[key] = val
    return md

def main(old_fp, new_fp):
    old_map = mapping_from_filepath(old_fp)
    new_map = mapping_from_filepath(new_fp)
    for ok, ov in old_map.items():
        try:
            nv = new_map[ok]
        except:
            print(f'split\t{ok} from {ov}')
            continue
        if nv == ov:
            continue
        try:
            nov = new_map[ov]
        except:
            print(f'split\t{ok} from {ov}')
            print(f'merge\t{ok} -> {nv}')
            continue
        if nov != nv:
            print(f'old[{ok}] -> {ov}\tnew[{ok}]->{nv}  new_map[{ov}]->{nov}')
    for nk, nv in new_map.items():
        if nk  in old_map:
            continue
        print(f'merge\t{nk} -> {nv}')
        
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])


