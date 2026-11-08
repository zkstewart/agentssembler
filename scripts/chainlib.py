# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Shared helpers: lift unitig coordinates through a chain AGP.'''
import collections
def read_agp(path):
    '''unitig -> list of (chain, chain_start0, comp_start0, comp_end0_exclusive, strand)'''
    P = collections.defaultdict(list); L = {}
    for l in open(path):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        L[s[0]] = max(L.get(s[0], 0), int(s[2]))
        if s[4] in ("U", "N"): continue
        P[s[5]].append((s[0], int(s[1]) - 1, int(s[6]) - 1, int(s[7]), s[8]))
    return P, L
def lift(P, u, x):
    '''unitig 0-based position -> (chain, position), or (None, None) outside the kept range'''
    for c, cs, a, b, st in P.get(u, ()):
        if a <= x < b: return c, (cs + x - a if st == "+" else cs + b - 1 - x)
    return None, None
