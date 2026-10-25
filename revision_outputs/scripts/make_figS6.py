#!/usr/bin/env python3
"""
Fig. S6 - structural modeling and de novo candidates (submission PDF), 2026-10-25

(A) AlphaFold2 model of WBP5 with the DoGSiteScorer P0 pocket, as exported
    from the ProteinsPlus viewer (Fig7A.png, 3828 x 1785 px, 600 dpi).
    Embedded at its native resolution; no resampling.
(B) 2D structures drawn with RDKit from the SMILES deposited in
    all_candidates_corrected.smi, labelled with the AutoDock Vina mean
    binding energy (5 seeds) from docking_P0/docking_replicates_P0_summary.csv.
    Drawn as vector graphics.

Usage: python make_figS6.py <Fig7A.png> <repo_root> <out_dir>
"""
import csv
import os
import re
import sys

import cairosvg
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont('LibSansBold', '/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf'))

FIG_A, REPO, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

# ---- inputs read from the repository (nothing typed by hand) -------------
smiles = {}
with open(os.path.join(REPO, "all_candidates_corrected.smi")) as fh:
    for line in fh:
        if line.strip():
            s, name = line.split()
            smiles[name] = s
energy = {}
with open(os.path.join(REPO, "docking_P0", "docking_replicates_P0_summary.csv")) as fh:
    for row in csv.DictReader(fh):
        energy[row["ligand"]] = float(row["mean_kcal_per_mol"])

ORDER = ["Pazopanib_Ref", "Candidate_4", "Candidate_2",
         "Candidate_3", "Candidate_1", "Candidate_5"]
LABEL = {"Pazopanib_Ref": "Pazopanib"}
ORDER.sort(key=lambda k: energy[k])          # most favourable first

# ---- panel B as SVG -> vector PDF ----------------------------------------
CELL_W, CELL_H, CAP_H, COLS, ROWS = 340, 260, 34, 3, 2
PAD_L, PAD_T = 8, 32
FONT = "Liberation Sans, Arial, Helvetica, sans-serif"


def draw(smi):
    m = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(CELL_W, CELL_H)
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 1.6
    o.minFontSize = 12
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    inner = d.GetDrawingText().split("<!-- END OF HEADER -->", 1)[1]
    inner = re.sub(r"<rect[^>]*>\s*</rect>", "", inner, count=1)
    return inner.replace("</svg>", "")


W = PAD_L * 2 + COLS * CELL_W
H = PAD_T + ROWS * (CELL_H + CAP_H) + 8
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}">',
       '<rect width="100%" height="100%" fill="white"/>',
       f'<text x="{PAD_L}" y="26" font-family="{FONT}" font-size="26" '
       f'font-weight="bold">(B)</text>']
for i, k in enumerate(ORDER):
    r, c = divmod(i, COLS)
    x, y = PAD_L + c * CELL_W, PAD_T + r * (CELL_H + CAP_H)
    svg.append(f'<g transform="translate({x},{y})">{draw(smiles[k])}</g>')
    svg.append(f'<text x="{x + CELL_W / 2}" y="{y + CELL_H + 20}" '
               f'text-anchor="middle" font-family="{FONT}" font-size="16">'
               f'{LABEL.get(k, k)} (−{abs(energy[k]):.2f})</text>')
svg.append("</svg>")
b_pdf = os.path.join(OUT, "_S6B.pdf")
cairosvg.svg2pdf(bytestring="\n".join(svg).encode(), write_to=b_pdf)

# ---- compose A (raster, native pixels) above B (vector) ------------------
imA = Image.open(FIG_A)
pB = PdfReader(b_pdf).pages[0]
bw, bh = float(pB.mediabox.width), float(pB.mediabox.height)
PAGE_W = bw                                  # width set by panel B
aw = PAGE_W - 16
ah = aw * imA.height / imA.width
GAP, TOP = 6, 4
PAGE_H = TOP + ah + GAP + bh

a_pdf = os.path.join(OUT, "_S6A.pdf")
cv = canvas.Canvas(a_pdf, pagesize=(PAGE_W, PAGE_H), initialFontName="LibSansBold")
cv.drawImage(ImageReader(imA), 8, PAGE_H - TOP - ah, width=aw, height=ah)
# Fig7A.png already carries its own "(A)" label
cv.save()

page = PdfReader(a_pdf).pages[0]
page.merge_transformed_page(pB, Transformation().translate(0, 0))
w = PdfWriter()
w.add_page(page)
with open(os.path.join(OUT, "Figure_S6.pdf"), "wb") as fh:
    w.write(fh)
for f in (a_pdf, b_pdf):
    os.remove(f)

print("panel A", imA.size, "-> %.0f ppi in the PDF" % (imA.width / (aw / 72)))
for k in ORDER:
    print(f"  {LABEL.get(k, k):<12} {energy[k]:.2f}  {smiles[k]}")
