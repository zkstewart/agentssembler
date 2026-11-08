#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
FASTA of the high-copy chains that round2.py removed. They stay in the assembly as unplaced sequence, named repeat_<chain>.
Usage: extract_repeat_chains.py repeat_chains.txt chains.fasta out.fasta
'''
import sys
want = {l.split()[0] for l in open(sys.argv[1]) if l.strip() and l[0] != "#"}
w = False
with open(sys.argv[3], "w") as o:
    for l in open(sys.argv[2]):
        if l[0] == ">":
            c = l[1:].split()[0]; w = c in want
            if w: o.write(f">repeat_{c}\n")
        elif w: o.write(l)
