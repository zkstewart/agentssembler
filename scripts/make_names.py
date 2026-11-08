#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Names for groups, and which to write in the reference's direction. Both are conventions and do not change any layout.
Takes the outputs of eval_groups.py and eval_order.py; the haplotype suffix comes from the baseline if one was given to
eval_groups.py, otherwise it is arbitrary. Check the names are unique before using them.
Usage: make_names.py eval_groups.txt eval_order.txt out.tsv
'''
import re, sys
ud = {}
for l in open(sys.argv[2]):
    m = re.match(r"(g\d+) .*up/down (\d+) (\d+)", l)
    if m: ud[m[1]] = (int(m[2]), int(m[3]))
names = {}
with open(sys.argv[3], "w") as o:
    o.write("#group\tname\tflip\n")
    for l in open(sys.argv[1]):
        s = l.rstrip("\n").split("\t")
        if re.match(r"g\d+$", s[0]):
            names[s[0]] = s[5]; u, d = ud.get(s[0], (1, 0))
            o.write(f"{s[0]}\t{s[5]}\t{'flip' if d > u else 'keep'}\n")
print(len(names), "groups,", len(set(names.values())), "unique names", file=sys.stderr)
