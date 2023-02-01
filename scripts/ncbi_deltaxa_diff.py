#!/usr/bin/env python3
import sys

def ids_from_filepath(fp):
    idset = set()
    with open(fp, 'r') as inp:
        for line in inp:
            ls = line.strip()
            if not ls:
                continue
            assert ls.endswith('|')
            id_str = ls[:-1].strip()
            id_int = int(id_str)
            assert id_int not in idset
            idset.add(id_int)
    return idset

def main(old_fp, new_fp):
    old_ids = ids_from_filepath(old_fp)
    new_ids = ids_from_filepath(new_fp)
    undel = old_ids.difference(new_ids)
    if undel:
        print("Undeleted {}".format(', '.join([str(i) for i in undel])))
    newdel = list(new_ids.difference(old_ids))
    newdel.sort()
    print('\n'.join([str(i) for i in newdel]))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])


