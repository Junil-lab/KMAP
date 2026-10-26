#!/usr/bin/env python3
"""Make synthetic FLNC PacBio reads and a small matching genome/annotation."""
import csv
from pathlib import Path
import random
import sys
import pysam

root = Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261021)
codons = [a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT' if a+b+c not in ('TAA','TAG','TGA')]
sequence = list(''.join(rng.choices('ACGT', k=100000)))
models, genes = [], []
for gene in range(12):
    start = 1000 + gene * 7000
    exons = [(start + i*600, start+i*600+300) for i in range(5)]
    for index, (left, right) in enumerate(exons):
        exon = ''.join(rng.choices(codons, k=100))
        if index == 0:
            exon = 'ATG' + exon[3:]
        if index == 4:
            exon = exon[:-3] + 'TAA'
        sequence[left:right] = exon
        sequence[left-2:left] = 'AG'
        sequence[right:right+2] = 'GT'
    genes.append((gene, exons))
    for label, order, count, annotated in [('A', [0,1,2,4], 12+gene, True),
                                          ('B', [0,1,3,4], 8+gene%4, True),
                                          ('NOVEL', [0,2,4], 6+gene%3, False)]:
        models.append((gene, label, [exons[i] for i in order], count, annotated))
genome = ''.join(sequence)
fasta = root/'genome.fa'
with fasta.open('w') as out:
    out.write('>chrSmoke\n')
    for i in range(0, len(genome), 80): out.write(genome[i:i+80]+'\n')
pysam.faidx(str(fasta))
with (root/'annotation.gtf').open('w') as out:
    def feature(kind, left, right, attributes):
        out.write(f'chrSmoke\tsynthetic\t{kind}\t{left+1}\t{right}\t.\t+\t.\t{attributes}\n')
    for gene, exons in genes:
        attrs = f'gene_id "G{gene}"; gene_name "Gene{gene}";'
        feature('gene',exons[0][0],exons[-1][1],attrs)
        for g, label, tx_exons, count, annotated in models:
            if g != gene or not annotated: continue
            tx = attrs + f' transcript_id "G{gene}_{label}";'
            feature('transcript',tx_exons[0][0],tx_exons[-1][1],tx)
            for i,(left,right) in enumerate(tx_exons):
                feature('exon',left,right,tx+f' exon_number "{i+1}";')
header = {'HD':{'VN':'1.6','SO':'unknown','pb':'3.0.1'},
          'RG':[{'ID':'12345678','PL':'PACBIO','PU':'mSmoke','SM':'Smoke',
                 'DS':'READTYPE=CCS;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-789-300;BASECALLERVERSION=5.0.0'}]}
with pysam.AlignmentFile(str(root/'flnc.bam'), 'wb', header=header) as out:
    zm=1
    for gene,label,exons,count,annotated in models:
        read=''.join(genome[a:b] for a,b in exons)
        for _ in range(count):
            rec=pysam.AlignedSegment()
            rec.query_name=f'mSmoke/{zm}/ccs'
            rec.query_sequence=read
            rec.flag=4
            rec.query_qualities=pysam.qualitystring_to_array('I'*len(read))
            rec.set_tag('RG','12345678')
            rec.set_tag('zm',zm)
            rec.set_tag('rq',0.999)
            rec.set_tag('np',12)
            out.write(rec)
            zm+=1
(root/'polyA_motifs.txt').write_text('AATAAA\nATTAAA\n')
with (root/'cage.bed').open('w') as cage, (root/'polya.bed').open('w') as polya:
    for gene, exons in genes:
        cage.write(f'chrSmoke\t{exons[0][0]-10}\t{exons[0][0]+10}\tCAGE{gene}\t1\t+\n')
        polya.write(f'chrSmoke\t{exons[-1][1]-10}\t{exons[-1][1]+10}\tpolyA{gene}\t1\t+\n')
with (root/'samples.csv').open('w',newline='') as out:
    writer=csv.writer(out)
    writer.writerow(['sample_id','flnc_bam','tissue'])
    # Test CSV quoting, relative BAM paths, and shell-safe metadata too.
    writer.writerow(['SmokeA','flnc.bam', 'synthetic, tissue'])
    writer.writerow(['SmokeB','flnc.bam', 'synthetic "quoted" tissue'])
print(f'{root}: {zm-1} synthetic reads/sample, {len(models)} transcript models')
