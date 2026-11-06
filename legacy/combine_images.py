#!/usr/bin/env python3
'''
Combine contact-map PNGs into one image.

    combine_images.py grid <out.png> <a.png> <b.png> ...   images in a square-ish grid, in the
                                                         order given, each resized to 700 px
    combine_images.py row  <out.png> <a.png> <b.png> ...   images side by side (e.g. before/after)
'''
import sys, math
from PIL import Image

def grid(out, files, size=700):
    ims = [ Image.open(f).convert("RGB").resize((size, size)) for f in files ]
    cols = math.ceil(math.sqrt(len(ims)))
    rows = math.ceil(len(ims) / cols)
    g = Image.new("RGB", (cols * size, rows * size), "white")
    for k, im in enumerate(ims):
        g.paste(im, ((k % cols) * size, (k // cols) * size))
    g.save(out)

def row(out, files):
    ims = [ Image.open(f).convert("RGB") for f in files ]
    g = Image.new("RGB", (sum(im.width for im in ims), max(im.height for im in ims)), "white")
    x = 0
    for im in ims:
        g.paste(im, (x, 0)); x += im.width
    g.save(out)

if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] in ("grid", "row"):
        (grid if sys.argv[1] == "grid" else row)(sys.argv[2], sys.argv[3:])
    else:
        sys.exit(__doc__)
