"""Draw the three pipeline overview diagrams (PPT style) as SVG."""
import sys
from xml.sax.saxutils import escape

TEAL, ORANGE, GRAY = "#2a9d8f", "#ee8a2a", "#6b7280"
INK, SUB, LINE = "#1f2937", "#4b5563", "#475569"
FONT = "Helvetica, Arial, sans-serif"


class Svg:
    def __init__(self, w, h):
        self.w, self.h, self.items = w, h, []

    def add(self, s):
        self.items.append(s)

    def text(self, x, y, s, size=12, weight="normal", color=INK, anchor="start", mono=False):
        fam = "Menlo, Consolas, monospace" if mono else FONT
        self.add(f'<text x="{x}" y="{y}" font-family="{fam}" font-size="{size}" font-weight="{weight}" '
                 f'fill="{color}" text-anchor="{anchor}">{escape(s)}</text>')

    def box(self, cx, y, w, num, title, files, kind="step", h=None, dashed=False):
        color = {"step": TEAL, "out": ORANGE, "ref": GRAY}[kind]
        files = [files] if isinstance(files, str) else files
        h = h or 44 + 15 * len(files)
        x = cx - w / 2
        dash = ' stroke-dasharray="5 4"' if dashed else ""
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#ffffff" stroke="{color}" stroke-width="1.6"{dash}/>')
        self.add(f'<rect x="{x}" y="{y}" width="6" height="{h}" fill="{color}"/>')
        numcol = {"step": "#1f7a6f", "out": "#c96a12", "ref": GRAY}[kind]
        self.text(x + 16, y + 16, num, 11, "bold", numcol)
        self.text(x + 16, y + 34, title, 14, "bold", INK)
        for i, f in enumerate(files):
            self.text(x + 16, y + 51 + 15 * i, f, 11, "normal", numcol if kind != "ref" else SUB)
        return {"cx": cx, "top": y, "bot": y + h, "left": x, "right": x + w, "mid": y + h / 2}

    def pill(self, cx, cy, label, w=None):
        w = w or 16 + 7.0 * len(label)
        self.add(f'<rect x="{cx - w / 2}" y="{cy - 11}" width="{w}" height="22" rx="11" fill="#eef2f7" stroke="#cbd5e1"/>')
        self.text(cx, cy + 4, label, 11.5, "bold", "#334155", "middle")

    def path(self, pts, arrow=True, dashed=False, color=LINE):
        d = "M " + " L ".join(f"{x} {y}" for x, y in pts)
        dash = ' stroke-dasharray="5 4"' if dashed else ""
        mk = ' marker-end="url(#arr)"' if arrow else ""
        self.add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.5"{dash}{mk}/>')

    def legend(self, x, y, entries):
        self.add(f'<rect x="{x}" y="{y}" width="170" height="{16 + 20 * len(entries)}" fill="#ffffff" stroke="#e5e7eb"/>')
        for i, (color, label, dashed) in enumerate(entries):
            yy = y + 12 + 20 * i
            dash = ' stroke-dasharray="3 2"' if dashed else ""
            self.add(f'<rect x="{x + 12}" y="{yy}" width="14" height="14" fill="#ffffff" stroke="{color}" stroke-width="1.5"{dash}/>')
            self.add(f'<rect x="{x + 12}" y="{yy}" width="4" height="14" fill="{color}"/>')
            self.text(x + 34, yy + 11, label, 11.5, "normal", SUB)

    def save(self, path):
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" viewBox="0 0 {self.w} {self.h}">'
                f'<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
                f'<path d="M0 0 L10 5 L0 10 z" fill="{LINE}"/></marker></defs>'
                f'<rect x="0" y="0" width="{self.w}" height="{self.h}" rx="10" fill="#ffffff" stroke="#d0d7de"/>')
        with open(path, "w") as f:
            f.write(head + "\n".join(self.items) + "</svg>\n")


LEGEND = [(TEAL, "Pipeline step", False), (ORANGE, "Final output", False), (GRAY, "Input / reference", False)]


def wgs(path):
    s = Svg(1100, 860)
    s.legend(910, 20, LEGEND)
    C = 550
    b0 = s.box(C, 24, 310, "00", "HiFi BAM with kinetics", "{sample_id}.hifi_reads.bam")
    s.pill(C, 112, "jasmine")
    s.path([(C, b0["bot"]), (C, 101)], arrow=False); s.path([(C, 123), (C, 140)])
    b1 = s.box(C, 142, 310, "01", "Methyl-tagged BAM", "{sample_id}.methyl.bam")
    alt = s.box(905, 142, 300, "input", "HiFi BAM with 5mC tags", "start_step=pbmm2 skips Jasmine", "ref")
    s.path([(alt["left"], alt["mid"] + 22), (C + 40, alt["mid"] + 22), (C + 40, 248)], dashed=True)
    s.pill(C, 230, "pbmm2")
    ref = s.box(175, 208, 290, "ref", "Reference genome", "GRCh38_no_alt_analysis_set", "ref")
    s.path([(ref["right"], 230), (C - 30, 230)])
    s.path([(C, b1["bot"]), (C, 219)], arrow=False); s.path([(C, 241), (C, 258)])
    b2 = s.box(C, 260, 310, "02", "Aligned BAM", "{sample_id}.aligned.sorted.bam")
    sex = s.box(905, 260, 300, "sex", "Samplesheet sex: XX / XY", "haploid chrX / chrY calling", "ref")
    cols = [190, 550, 910]
    s.path([(C, b2["bot"]), (C, 362)], arrow=False)
    s.path([(cols[0], 362), (cols[2], 362)], arrow=False)
    s.path([(sex["left"] + 30, sex["bot"]), (sex["left"] + 30, 362)], arrow=False, dashed=True, color=GRAY)
    for x, tool in zip(cols, ["DeepVariant", "Sawfish", "TRGT"]):
        s.pill(x if x != C else x - 95, 343, tool)
        s.path([(x, 362), (x, 398)])
    b3 = s.box(cols[0], 400, 300, "03", "VCF · Small variants", "{sample_id}.deepvariant.vcf.gz")
    s.box(cols[1], 400, 300, "04", "VCF · Structural variants", "{sample_id}.sv.vcf.gz")
    s.box(cols[2], 400, 300, "05", "VCF · Tandem repeats", "{sample_id}.trgt.vcf.gz")
    yb = b3["bot"]
    for x in cols:
        s.path([(x, yb), (x, 500)], arrow=False)
    four = [160, 420, 680, 940]
    s.path([(four[0], 500), (four[3], 500)], arrow=False)
    s.pill(C, 500, "HiPhase")
    for x in four:
        s.path([(x, 500), (x, 538)])
    a = s.box(four[0], 540, 240, "06a", "Phased BAM", "{sample_id}.haplotagged.bam", "out")
    s.box(four[1], 540, 240, "06b", "Phased SNV / indel", "{sample_id}.snv.phased.vcf.gz", "out")
    s.box(four[2], 540, 240, "06c", "Phased SV", "{sample_id}.sv.phased.vcf.gz", "out")
    s.box(four[3], 540, 240, "06d", "Phased tandem repeats", "{sample_id}.tr.phased.vcf.gz", "out")
    s.path([(four[0], a["bot"]), (four[0], 660)])
    s.path([(four[0], 630), (four[1], 630), (four[1], 660)])
    s.pill(four[0], 645, "pb-CpG-tools")
    s.pill(four[1], 645, "MethBat")
    s.box(four[0], 662, 240, "07", "CpG methylation", ["{sample_id}.combined.bed.gz", "+ hap1 / hap2"], "out")
    s.box(four[1], 662, 240, "07b", "Methylation pileup", ["{sample_id}.5mC / 5hmC / 6mA.bed.gz", "cpg_island_profile.tsv"], "out")
    s.path([(four[3], 599), (four[3], 660)])
    s.box(four[3], 662, 240, "QC", "QC summary", "qc_summary.csv", "out")
    s.text(30, 780, "start_step can also be deepvariant (aligned BAM) or hiphase (existing VCFs).", 12, "normal", SUB)
    s.text(30, 800, "07 and 07b run side by side on the same haplotagged BAM.", 12, "normal", SUB)
    s.h = 825
    s.save(path)


def rna(path):
    s = Svg(1100, 830)
    s.legend(910, 20, LEGEND)
    C = 550
    b0 = s.box(C, 24, 310, "00", "FLNC BAM", "{sample_id} FLNC BAM, one per sample")
    s.pill(C, 112, "pbmm2 (ISOSEQ preset)")
    ref = s.box(175, 90, 290, "ref", "hg38 + GENCODE v45", "genome FASTA + annotation GTF", "ref")
    s.path([(ref["right"], 112), (C - 85, 112)])
    s.path([(C, b0["bot"]), (C, 101)], arrow=False); s.path([(C, 123), (C, 140)])
    b1 = s.box(C, 142, 310, "01", "Aligned BAM", "alignment/pbmm2/ (shared by all callers)")
    cols = [190, 550, 910]
    s.path([(C, b1["bot"]), (C, 245)], arrow=False)
    s.path([(cols[0], 245), (cols[2], 245)], arrow=False)
    for x, tool in zip(cols, ["IsoQuant", "isocall", "FLAIR"]):
        s.pill(x if x != C else x - 70, 226, tool)
        s.path([(x, 245), (x, 278)])
    s.box(cols[0], 280, 300, "02a", "IsoQuant transcripts", "isoquant/discovery/")
    s.box(cols[1], 280, 300, "02b", "isocall transcripts", "isocall/{profile,discovery}/")
    s.box(cols[2], 280, 300, "02c", "FLAIR transcripts", "flair/discovery/")
    # cohort: only with two or more samples
    co = s.box(730, 400, 250, "only with 2+ samples", "Cohort transcriptome", ["cohort/isocall, cohort/flair", "joint merge across samples"], "step", dashed=True)
    s.path([(cols[1], 339), (cols[1], 380), (co["left"] + 40, 380), (co["left"] + 40, 398)], dashed=True)
    s.path([(cols[2], 339), (cols[2], 380), (co["right"] - 40, 380), (co["right"] - 40, 398)], dashed=True)
    y_m = 545
    s.path([(cols[0], 339), (cols[0], y_m)], arrow=False)
    s.path([(cols[1] - 120, 339), (cols[1] - 120, y_m)], arrow=False)
    s.path([(cols[2] + 120, 339), (cols[2] + 120, y_m)], arrow=False)
    s.path([(co["cx"], co["bot"]), (co["cx"], y_m)], arrow=False)
    s.path([(cols[0], y_m), (cols[2] + 120, y_m)], arrow=False)
    s.pill(C, y_m, "SQANTI3 QC + ORF (TD2 / PSAURON)", 270)
    outs = [190, 550, 910]
    for x in outs:
        s.path([(x if x != C else C, y_m + 11 if x == C else y_m), (x, 598)])
    s.box(outs[0], 600, 300, "03", "Sample results", ["samples/<id>/<method>/sqanti3/", "classification, ORF / protein FASTA"], "out")
    s.box(outs[1], 600, 300, "04", "Cohort results", ["cohort/<method>/sqanti3/", "joint transcriptome + counts"], "out")
    s.box(outs[2], 600, 300, "05", "Reports and QC summary", ["HTML report per method", "qc_summary/qc_summary.csv"], "out")
    s.text(30, 735, "Choose callers with --methods isoquant,isocall,flair. Cohort steps run only when the samplesheet has 2 or more samples", 12, "normal", SUB)
    s.text(30, 755, "(--run_cohort false to skip). IsoQuant always runs per sample.", 12, "normal", SUB)
    s.h = 780
    s.save(path)


def scrna(path):
    s = Svg(1100, 1120)
    s.legend(910, 20, LEGEND)
    C = 470
    W = 320
    y = 24
    b = s.box(C, y, W, "00", "FASTQ directory", "sample_R1 / R2 / I1 / I2 fastq.gz")
    ref = s.box(165, 98, 270, "ref", "10x reference", ["user-supplied transcriptome", "refdata-gex-GRCh38-2024-A"], "ref")
    steps = [
        ("cellranger", "01", "Cell Ranger count", ["raw_feature_bc_matrix.h5"]),
        ("CellBender", "02", "Remove background", ["{sample}.cellbender_cleaned.h5"]),
        ("Scanpy", "03", "Scanpy QC for SOLO", ["cell QC, MT filter off by default; all genes kept", "full.h5ad + hvg2000.h5ad + QC CSV"]),
        ("SOLO", "04", "SOLO doublet scoring", ["{sample}.solo/soloed.h5ad (scores corrected)"]),
        ("filter_singlets", "05", "Filter SOLO singlets", ["annotated full + singlet-only H5AD"]),
        ("merge_qc", "06", "QC summary merge", ["QC metrics + SOLO summary"]),
    ]
    prev = b
    boxes = {}
    for tool, num, title, files in steps:
        py = prev["bot"] + 26
        s.path([(C, prev["bot"]), (C, py - 11)], arrow=False)
        s.pill(C, py, tool)
        s.path([(C, py + 11), (C, py + 34)])
        nb = s.box(C, py + 36, W, num, title, files)
        boxes[num] = (py, nb)
        prev = nb
    py01 = boxes["01"][0]
    s.path([(ref["right"], py01), (C - 45, py01)])
    # start_step entry points
    py02, py03 = boxes["02"][0], boxes["03"][0]
    e1 = s.box(870, py02 - 30, 300, "start_step=cellbender", "Raw matrix H5", "raw_feature_bc_matrix.h5", "ref")
    s.path([(e1["left"], py02), (C + 50, py02)], dashed=True)
    e2 = s.box(870, py03 - 30, 300, "start_step=qc", "CellBender H5", "cellbender_cleaned.h5", "ref")
    s.path([(e2["left"], py03), (C + 35, py03)], dashed=True)
    # final outputs
    n5, n6 = boxes["05"][1], boxes["06"][1]
    o1 = s.box(870, n5["top"], 300, "final", "Final cell-level H5ADs", ["{sample}_qc_solo_annotated_full.h5ad", "{sample}_qc_solo_singlets_full.h5ad"], "out")
    s.path([(n5["right"], n5["mid"]), (o1["left"], n5["mid"])])
    o2 = s.box(870, n6["top"], 300, "final", "Final QC table", "scrnaseq_qc_summary.csv", "out")
    s.path([(n6["right"], n6["mid"]), (o2["left"], n6["mid"])])
    # bottom panels
    yp = n6["bot"] + 40
    panels = [
        (GRAY, "Required inputs", ["10x FASTQ directory", "10x Cell Ranger reference", "samplesheet CSV (sample_id, fastqs,", "  sample, transcriptome, start_step)"]),
        (TEAL, "Main processing", ["Cell Ranger count", "CellBender ambient RNA correction", "Scanpy QC and HVG selection", "SOLO doublet detection and filtering"]),
        (ORANGE, "Key final outputs", ["SOLO-annotated full H5AD", "SOLO singlet-only full H5AD", "Integrated QC summary CSV"]),
    ]
    for i, (col, title, lines) in enumerate(panels):
        x = 30 + i * 352
        s.add(f'<rect x="{x}" y="{yp}" width="335" height="128" fill="#ffffff" stroke="#e5e7eb"/>')
        s.add(f'<rect x="{x}" y="{yp}" width="5" height="128" fill="{col}"/>')
        s.text(x + 20, yp + 26, title, 14, "bold", INK)
        for j, line in enumerate(lines):
            s.text(x + 20, yp + 50 + 19 * j, ("• " if not line.startswith("  ") else "   ") + line.strip(), 12.5, "normal", SUB)
    s.h = yp + 150
    s.save(path)


if __name__ == "__main__":
    out = sys.argv[1]
    wgs(f"{out}/pblongreadwgs_pipeline.svg")
    rna(f"{out}/pblongreadrna_pipeline.svg")
    scrna(f"{out}/scrna_pipeline.svg")
