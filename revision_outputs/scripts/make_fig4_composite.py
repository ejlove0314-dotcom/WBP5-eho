#!/usr/bin/env python3
"""
Fig. 4 - assembled from the GSE198315 Run B panel outputs (2026-10-25)
=====================================================================
Panel  Source (Run B output, 2026-08-05)      Embedded as
  A    Fig4A_runB.pdf                          vector
  B    Fig4B_runB.png  (600 dpi)               raster, native pixels
  C    Fig4C_runB.png  (600 dpi)               raster, native pixels
  D    Fig4D_runB.png  (600 dpi)               raster, native pixels

These are the panels used in the manuscript since the original submission
(identical content to manuscript v21, Fig. 4A-D). Nothing inside a panel is
redrawn; this script only places the panels and adds the letters (A)-(D).
The panel values were independently recomputed from obs.csv on 2026-10-24
(revision_outputs/gse198315/J1_Fig4A_proxy_compartment_composition.csv).

Usage: python make_fig4_composite.py <panel_dir> <out_dir>
"""
import os
import sys

from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
pdfmetrics.registerFont(TTFont("LSB", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"))

MM = 72 / 25.4
PAGE_W = 180 * MM
LEFT = 6 * MM            # room for the panel letter
GAP = 4 * MM             # vertical gap between panels
LAB = 5 * MM             # height reserved above each panel for its letter

# panel widths (mm) chosen to keep the text in each panel at a similar size
pA = PdfReader(os.path.join(SRC, "Fig4A_runB.pdf")).pages[0]
aw_pt, ah_pt = float(pA.mediabox.width), float(pA.mediabox.height)
A_W = PAGE_W - LEFT
A_S = A_W / aw_pt
A_H = ah_pt * A_S

rasters = []
for lab, fn, w_mm in [("B", "Fig4B_runB.png", 120), ("C", "Fig4C_runB.png", 150),
                      ("D", "Fig4D_runB.png", 174)]:
    im = Image.open(os.path.join(SRC, fn)).convert("RGB")
    w = w_mm * MM
    rasters.append((lab, im, w, w * im.height / im.width))

PAGE_H = LAB + A_H + sum(GAP + LAB + h for _, _, _, h in rasters) + 2 * MM

base = os.path.join(OUT, "_fig4_base.pdf")
cv = canvas.Canvas(base, pagesize=(PAGE_W, PAGE_H), initialFontName="LSB")
cv.setFont("LSB", 12)
y = PAGE_H - LAB
cv.drawString(0, y + 1.2 * MM, "(A)")
y -= A_H
yA = y                                           # bottom of panel A
for lab, im, w, h in rasters:
    y -= GAP + LAB
    cv.drawString(0, y + 1.2 * MM, f"({lab})")
    y -= h
    cv.drawImage(ImageReader(im), LEFT, y, width=w, height=h)
    print(f"panel {lab}: {im.size[0]}x{im.size[1]} px -> {w / MM:.0f} mm wide, "
          f"{im.size[0] / (w / 72):.0f} ppi")
cv.save()

page = PdfReader(base).pages[0]
page.merge_transformed_page(pA, Transformation().scale(A_S).translate(LEFT, yA))
wr = PdfWriter()
wr.add_page(page)
with open(os.path.join(OUT, "Figure_4.pdf"), "wb") as fh:
    wr.write(fh)
os.remove(base)
print(f"Figure_4.pdf: {PAGE_W / MM:.0f} x {PAGE_H / MM:.0f} mm")
