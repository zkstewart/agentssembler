'''
Helpers shared by the pipeline scripts.

Group names are '<chrom>h<hap>' (e.g. chr4h3) and final pseudochromosome names are
'<chrom>_h<hap>' (e.g. chr4_h3), where <chrom> is the reference sequence name and <hap> the
haplotype number (1..ploidy).
'''
import re

def natural_key(name):
    '''Sort key that orders embedded numbers numerically: chr2 < chr10, chr4_h3 < chr4_h10.'''
    return [ int(x) if x.isdigit() else x for x in re.split(r"(\d+)", name) ]

def split_group(group):
    '''chr4h3 -> ("chr4", "3")'''
    chrom, hap = group.rsplit("h", 1)
    return chrom, hap

def re_site_regex(sites):
    '''"GATC,GANTC" -> a lookahead regex counting every site, with N matching any base.'''
    alts = [ s.upper().replace("N", "[ACGT]") for s in sites.split(",") if s ]
    return re.compile(r"(?=(" + "|".join(alts) + r"))")
