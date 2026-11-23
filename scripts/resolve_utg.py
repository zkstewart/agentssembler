#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Write an edit table with every "utg:<unitig>" name resolved to the chain holding that unitig (chainlib.utg_resolver), so that an
edit naming a connector or arm by its unitig survives the renumbering a reset causes (TRAPS T16). Comment lines are dropped; only
--cols are resolved (1-based, default 1,2), the rest pass unchanged, and with --keep only the first --keep columns are written.
Usage: resolve_utg.py chains.agp edits.tsv [--cols 1,2] [--keep N] > resolved.tsv'''
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chainlib import utg_resolver
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("edits"); p.add_argument("--cols", default="1,2"); p.add_argument("--keep", type=int); a = p.parse_args()
res = utg_resolver(a.agp); cols = {int(x) - 1 for x in a.cols.split(",")}
for l in open(a.edits):
    if l[0] == "#" or not l.strip(): continue
    s = l.rstrip("\n").split("\t"); s = [res(x) if i in cols else x for i, x in enumerate(s)]
    print("\t".join(s[:a.keep] if a.keep else s))
