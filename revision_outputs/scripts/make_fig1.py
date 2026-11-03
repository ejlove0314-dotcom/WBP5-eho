#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fig1.py — Figure 1 (Computational workflow)
=================================================
Corresponds one-to-one with Fig. 1 of the manuscript and with its legend.
Every element is a native PowerPoint shape, so the figure can be edited
directly in PowerPoint.

Output: Figure_1_workflow_20261103.pptx
        (PDF and PNG are produced with LibreOffice:
         soffice --headless --convert-to pdf)

-- revision of 2026-11-03 --------------------------------------------------
Before this revision the deposited script did not reproduce the submitted
Figure_1.pdf.  The wording, the order of the tool boxes and the layout all
differed, so the script was first brought into agreement with the submitted
figure and the corrections were then applied on top of it.  It now generates
Fig. 1 as the manuscript presents it.

  1. On the Finding line, "marks ... shared by" was replaced with the wording
     of the manuscript legend.  The revision withdrew both expressions -- the
     Discussion states that the data do not establish that WBP5 marks the
     program, and the title now reads "expressed across" rather than "shared
     by" -- but they survived in the figure alone.  The fibroblast and
     cancer-cell positive fractions (51.97% and 68.89%) are those reported in
     the Results.
  2. "76-residue core, 29-104" became "full-length model; core 29-104".
     Pocket detection and docking were run on the full-length model: the
     deposited pocket P_1 contains residues 27 and 28, which lie outside the
     structured core (29-104), so the core alone cannot have been the input.
  3. The citation was corrected from [11] to [14].
  4. "Drugs-lib; 4,574 approved drugs" was restored to the MTiOpenScreen box
     and "5 approved drugs + 5 de novo candidates" to the AutoDock Vina box.
     Both appear in the submitted figure and were absent from the script.
  5. The tool boxes were returned to two rows of three, in the submitted
     order; the script had drawn them as a single row of six.
  6. The slide height was changed from 7.5 to 7.62 in, which matches the
     548.64 pt page of the submitted PDF and keeps the closing note from being
     clipped once the second row is restored.

-- revision of 2026-10-22 (earlier history) --------------------------------
  1. The Phase I conclusion "not enriched in malignant cells", which was
     incorrect, became "not malignant-cell-specific" together with the
     statement that expression is largely absent from lymphoid cells.  That
     wording was itself replaced by item 1 above.
  2. The MTiOpenScreen box was added; it had been missing.
  3. "Exploratory structural ligandability" became "structural and chemical".
  4. "do NOT validate" became "do not validate".
  5. Hyphens became en dashes (29-104, Phases II-III).
  6. "Credibility suite" became "Robustness checks", removing the duplicated
     cross-cohort entry.
  7. The residue range 38-58 was added to the P_0 box.
"""
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

C = lambda h: RGBColor.from_string(h)

BLUE_FILL, BLUE_LINE, BLUE_TEXT = C("EAF3FB"), C("2E75B6"), C("1F4E79")
PUR_HEAD, PUR_LINE, PUR_FILL, PUR_TEXT = C("7030A0"), C("7030A0"), C("E9E1F4"), C("4E2583")
PUR_SUBLN = C("B9A5D9")
GRN_HEAD, GRN_LINE, GRN_FILL, GRN_TEXT = C("2E7D32"), C("2E7D32"), C("E4F0DC"), C("2F5B22")
GRN_SUBLN = C("A9D08E")
WHITE = C("FFFFFF")

FONT = "Calibri"


def add(slide, shape, x, y, w, h, fill, line, lw=1.0, adj=None):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill is None:
        s.fill.background()
    else:
        s.fill.solid()
        s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line
        s.line.width = Pt(lw)
    s.shadow.inherit = False
    if adj is not None and s.adjustments:
        s.adjustments[0] = adj
    s.text_frame.word_wrap = True
    return s


def text(shape, lines, size=10, color=None, bold=False, italic=False,
         align=PP_ALIGN.CENTER, space=0):
    """lines: a str, or a list of (text, bold, italic, size) tuples."""
    tf = shape.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Inches(0.06)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    if isinstance(lines, str):
        lines = [lines]
    for i, ln in enumerate(lines):
        if isinstance(ln, str):
            ln = (ln, bold, italic, size)
        t, b, it, sz = ln
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space)
        r = p.add_run()
        r.text = t
        f = r.font
        f.name, f.size, f.bold, f.italic = FONT, Pt(sz), b, it
        f.color.rgb = color if color else C("000000")
    return shape


def main():
    prs = Presentation()
    # 7.62 in matches the submitted Figure_1.pdf, whose page is 548.64 pt tall;
    # the saved script used 7.5 in, which clipped the closing note once the
    # tool boxes were restored to two rows.
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.62)
    sl = prs.slides.add_slide(prs.slide_layouts[6])

    L, W = 0.35, 12.63           # left margin, overall width
    R = L + W

    # -- 1. the premise box across the top ---------------------------
    b = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L, 0.22, W, 0.52,
            BLUE_FILL, BLUE_LINE, 1.25, adj=0.22)
    text(b, "WBP5 (TCEAL9): prognostic candidate carried forward from our prior "
            "bulk-tissue study [14] — not an unbiased discovery screen",
         size=12, color=BLUE_TEXT, bold=True)

    # -- 2. the Phase I container ------------------------------------
    p1_y, p1_h = 0.95, 3.56
    add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L, p1_y, W, p1_h, None, PUR_LINE, 1.5, adj=0.05)

    hd = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L + 0.08, p1_y + 0.07, W - 0.16, 0.42,
             PUR_HEAD, PUR_HEAD, 1.0, adj=0.30)
    text(hd, "Phase I — Transcriptomic target nomination  (main text)",
         size=12.5, color=WHITE, bold=True)

    cw, gap = 5.97, 0.25         # width of each of the two columns, gap
    cx1, cx2 = L + 0.25, L + 0.25 + cw + gap

    r1y, r1h = 1.58, 0.74
    s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, cx1, r1y, cw, r1h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.14)
    text(s, [("Discovery cohort: GSE103322 (HNSCC; 5,902 cells)", True, False, 10.5),
             ("Exploratory cohort: GSE198315 (HPV-negative OSCC)", True, False, 10.5)],
         color=PUR_TEXT, space=2)

    s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, cx2, r1y, cw, r1h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.14)
    text(s, [("Within-patient co-expression", True, False, 10.5),
             ("+ compartment-resolved analysis", True, False, 10.5)],
         color=PUR_TEXT, space=2)

    r2y, r2h = 2.42, 0.82
    s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, cx1, r2y, cw, r2h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.13)
    text(s, [("Program-level p-EMT enrichment", True, False, 10.5),
             ("depth-adjusted; expression-matched null gene sets", False, False, 9.5)],
         color=PUR_TEXT, space=2)

    s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, cx2, r2y, cw, r2h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.13)
    text(s, [("Robustness checks", True, False, 10.5),
             ("leave-one-patient-out; gene-family specificity; Xq22 module "
              "adjustment; negative and positive controls", False, False, 9.5)],
         color=PUR_TEXT, space=2)

    tw = 6.55
    s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L + (W - tw) / 2, 3.33, tw, 0.46,
            PUR_FILL, PUR_SUBLN, 1.0, adj=0.22)
    text(s, "TCGA HNSC bulk cohort (n = 520; purity-adjusted) — complementary",
         size=10.5, color=PUR_TEXT, bold=True)

    s = add(sl, MSO_SHAPE.RECTANGLE, L + 0.25, 3.86, W - 0.50, 0.56, None, None)
    text(s, ["Finding: WBP5 is expressed across the malignant and mesenchymal "
             "compartments and is largely absent from lymphoid cells "
             "(51.97% of fibroblasts and 68.89% of cancer cells positive);",
             "its association with the p-EMT program is detectable at the level "
             "of the program rather than at the level of individual markers"],
         size=10, color=PUR_TEXT, italic=True, space=1)

    # -- 3. the Phases II-III container ------------------------------
    p2_y, p2_h = 4.70, 2.82
    add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L, p2_y, W, p2_h, None, GRN_LINE, 1.5, adj=0.07)

    hd = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE, L + 0.08, p2_y + 0.07, W - 0.16, 0.42,
             GRN_HEAD, GRN_HEAD, 1.0, adj=0.30)
    text(hd, "Phases II–III — Exploratory structural and chemical "
             "ligandability  (Supplementary Information)",
         size=12.5, color=WHITE, bold=True)

    boxes = [
        [("AlphaFold2 model", True, False, 10.5),
         ("AF-Q9UHQ7-F1-model_v6", False, False, 9.0),
         ("full-length model; core 29–104", False, False, 9.0)],
        [("DoGSiteScorer", True, False, 10.5),
         ("candidate pocket P₀", False, False, 9.0),
         ("residues 38–58; Drug Score 0.84", False, False, 9.0)],
        [("MTiOpenScreen", True, False, 10.5),
         ("library-scale calibration", False, False, 9.0),
         ("Drugs-lib; 4,574 approved drugs", False, False, 9.0)],
        [("REINVENT4", True, False, 10.5),
         ("reinforcement-learning", False, False, 9.0),
         ("de novo generation", False, False, 9.0)],
        [("AutoDock Vina", True, False, 10.5),
         ("docking of 5 approved drugs", False, False, 9.0),
         ("+ 5 de novo candidates", False, False, 9.0)],
        [("SwissADME / ProTox", True, False, 10.5),
         ("drug-likeness and", False, False, 9.0),
         ("predicted toxicity", False, False, 9.0)],
    ]
    # The submitted figure lays the six tools out as two rows of three.  The
    # saved script drew a single row of six, which is one of the ways it failed
    # to reproduce the figure it was deposited as the source of.
    cols, g = 3, 0.15
    bw = (W - 0.50 - g * (cols - 1)) / cols
    bh, vg, by0 = 0.80, 0.12, 5.16
    for i, lines in enumerate(boxes):
        r, c = divmod(i, cols)
        s = add(sl, MSO_SHAPE.ROUNDED_RECTANGLE,
                L + 0.25 + c * (bw + g), by0 + r * (bh + vg), bw, bh,
                GRN_FILL, GRN_SUBLN, 1.0, adj=0.10)
        text(s, lines, color=GRN_TEXT, space=1)

    s = add(sl, MSO_SHAPE.RECTANGLE, L + 0.25, 6.98, W - 0.50, 0.50, None, None)
    text(s, ["Computationally prioritized hypotheses: these results do not validate "
             "the transcriptomic target nomination,",
             "and experimental validation is required"],
         size=10, color=GRN_TEXT, italic=True, space=1)

    # -- 4. the arrows between the phases ----------------------------
    for y0, y1 in ((0.74, 0.95), (4.51, 4.70)):
        a = add(sl, MSO_SHAPE.DOWN_ARROW, L + W / 2 - 0.11, y0, 0.22, y1 - y0,
                C("A6A6A6"), None)
        a.line.fill.background()

    out = "Figure_1_workflow_20261103.pptx"
    prs.save(out)
    print("saved ->", out)


if __name__ == "__main__":
    main()
