#!/usr/bin/env python
import sys
import re
inp = open(sys.argv[1], 'rU')
leafset = set()
ps = []
for line in inp:
    if not line.startswith('(('):
        if line.strip():
            sys.exit(line)
            assert(False)
        continue
    ls = line[2:].strip()
    assert(ls.endswith(');'))
    ls = ls[:-2]
    assert(ls.startswith('ott'))
    cs = ls.split(')')
    assert(len(cs) == 2)
    assert(cs[1].startswith(','))
    cs[1] = cs[1][1:]
    ing = frozenset(cs[0].split(','))
    outg = frozenset(cs[1].split(','))
    for el in ing:
        assert el.startswith('ott')
    for el in outg:
        assert el.startswith('ott')
    leafset.update(ing)
    leafset.update(outg)
    ps.append((ing, outg))

def only_one_in_set(f, s, setToCheck):
    if f in setToCheck:
        if s not in setToCheck:
            return True
    elif s in setToCheck:
        return True
    return False
def pair_redundant(f, s, statements):
    for el in statements:
        ing, outg = el
        if only_one_in_set(f, s, ing):
            return False
        if len(ing) == 2 and f in ing:
            return False
        if only_one_in_set(f, s, outg):
            return False
    return True

def merge_redundant(p, ls, statements):
    sset = set()
    to_keep, to_delete = p
    sys.stderr.write('deleting {} redundant with {} \n'.format(to_delete, to_keep))
    red_statements = []
    ls.remove(to_delete)
    for s in statements:
        ing, outg = set(s[0]), set(s[1])
        if to_keep in ing:
            ing.remove(to_delete)
        if to_keep in outg:
            outg.remove(to_delete)
        news = (frozenset(ing), frozenset(outg))
        if news not in sset:
            sset.add(news)
            red_statements.append(news)
    return ls, red_statements

lsl = list(leafset)
foid = lsl.pop()
while lsl:
    red_pair = None
    for soid in lsl:
        if soid == foid:
            continue
        if pair_redundant(foid, soid, ps):
            red_pair = (foid, soid)
            break
        if red_pair is not None:
            break
    if red_pair is None:
        foid = lsl.pop()
    else:
        leafset, ps = merge_redundant(red_pair, leafset, ps)
        lsl.remove(soid)
        sys.stderr.write('foid = {} ls size = {} # statements = {}\n'.format(foid, len(leafset), len(ps)))


for statement in ps:
    ings = ','.join(statement[0])
    outs = ','.join(statement[1])
    print '(({}),{});'.format(ings, outs)