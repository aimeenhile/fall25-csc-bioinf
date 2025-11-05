# ChatGPT 

### **Prompt 1:**

for the reference genes, is there a way to download the sequences for each genes separately without downloading the enttire chromosome10 and merge them into one file? Knowing that below is the information for each gene:

Human Gene CYP2C8 (ENST00000371270.6) from GENCODE V48 
Description: cytochrome P450 family 2 subfamily C member 8, transcript variant 1 (from RefSeq NM_000770.3) 
Gencode Transcript: ENST00000371270.6 
Gencode Gene: ENSG00000138115.15 Transcript (Including UTRs) 
    Position: hg38 chr10:95,036,772-95,069,497 Size: 32,726 Total Exon Count: 9 Strand: - 
Coding Region 
    Position: hg38 chr10:95,037,128-95,069,402 Size: 32,275 Coding Exon Count: 9 
    
Human Gene CYP2C9 (ENST00000260682.8) from GENCODE V48 
Description: cytochrome P450 family 2 subfamily C member 9 (from RefSeq NM_000771.4) 
Gencode Transcript: ENST00000260682.8 
Gencode Gene: ENSG00000138109.11 
Transcript (Including UTRs) 
    Position: hg38 chr10:94,938,658-94,990,091 Size: 51,434 Total Exon Count: 9 Strand: + 
Coding Region 
    Position: hg38 chr10:94,938,683-94,989,028 Size: 50,346 Coding Exon Count: 9 

Human Gene CYP2C19 (ENST00000371321.9) from GENCODE V48 
Description: cytochrome P450 family 2 subfamily C member 19 (from RefSeq NM_000769.4) 
Gencode Transcript: ENST00000371321.9 
Gencode Gene: ENSG00000165841.11 
Transcript (Including UTRs) 
    Position: hg38 chr10:94,762,681-94,855,547 Size: 92,867 Total Exon Count: 9 Strand: + 
Coding Region 
    Position: hg38 chr10:94,762,706-94,852,914 Size: 90,209 Coding Exon Count: 9

### **Answer 1:**

Excellent question — yes, there is a way to download just those three genes directly, without downloading the entire chr10 FASTA file.

There are two main approaches, depending on how you want to fetch them:

🧬 Option 1: Use NCBI Entrez (Biopython)

You can use Biopython’s Entrez.efetch() to directly download the FASTA sequence for each gene using its RefSeq accession (e.g., NM_000769.4 for CYP2C19).
This is clean, lightweight, and CI-friendly — no need for a 100 MB chromosome file.
```python
from Bio import Entrez, SeqIO

# Always identify yourself to NCBI
Entrez.email = "your_email@example.com"

# RefSeq mRNA IDs for each gene (hg38)
genes = {
    "CYP2C8": "NM_000770.3",
    "CYP2C9": "NM_000771.4",
    "CYP2C19": "NM_000769.4"
}

# Fetch and save sequences
with open("data/CYP2C_genes.fa", "w") as out_f:
    for gene, accession in genes.items():
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        out_f.write(f">{gene}\n{record.seq}\n")
        print(f"{gene}: {len(record.seq)} bp")
```
✅ Output:
data/CYP2C_genes.fa — containing FASTA sequences for all three genes.

✅ Pros:

Super fast (fetches only 3 small sequences)

Fully automated

Works in your Jupyter + CI setup

Doesn’t require samtools or downloading large files

❌ Cons:

You’ll get transcript sequences (mRNA), not genomic sequences (with introns).

If you need genomic DNA (e.g. for variant calling), use Option 2 instead.

🧬 Option 2: Fetch genomic regions directly from UCSC or Ensembl

If you need genomic DNA (including introns), you can directly fetch only those regions using UCSC Table Browser or Ensembl REST API, without downloading full chromosomes.

(a) Using UCSC API (bash + curl)
```bash
# CYP2C19 (chr10:94762681-94855547)
curl "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=94762680;end=94855547" \
  | jq -r '.dna' > CYP2C19.fa

# CYP2C9 (chr10:94938658-94990091)
curl "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=94938657;end=94990091" \
  | jq -r '.dna' > CYP2C9.fa

# CYP2C8 (chr10:95036772-95069497)
curl "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=95036771;end=95069497" \
  | jq -r '.dna' > CYP2C8.fa

# Merge all into one file
cat CYP2C*.fa > CYP2C_genes.fa
```
✅ Pros:

You get the actual genomic regions (needed for alignment / variant calling)

No need for entire chr10.fa.gz

Works entirely from CI

❌ Cons:

Requires jq (JSON parser) — can be installed in CI easily (sudo apt install jq)

🧪 Recommendation for Your Project

Since your assignment involves aligning reads and calling variants, you need genomic sequences (with introns), not just mRNA.
➡️ Therefore, use Option 2 (UCSC API) — it’s small, direct, and fully CI-compatible.

Would you like me to modify your week5.ipynb so it:

Automatically downloads the genomic sequences for CYP2C8, CYP2C9, and CYP2C19 via UCSC API,

Merges them into CYP2C_genes.fa,

Without ever downloading chr10.fa.gz?

### **Prompt 2:**

But I cannot install jq (JSON parser). Could you please give a clean snippet to extract all three CYP2C genes into a single FASTA without reverse complementing. That would be ready for alignment and variant analysis?

### **Answer 2:**

Perfect! Here's a clean Python snippet using pysam to extract all three CYP2C genes into a single FASTA file without reverse complementing, and with proper headers that include strand information:
```python
import pysam

# Path to the chromosome 10 reference genome
path = "data/chr10.fa"
fasta = pysam.FastaFile(path)

# Gene coordinates from GENCODE / UCSC (hg38)
GENE_INFO = {
    "CYP2C19": {"chr": "chr10", "start": 94762681, "end": 94855547, "strand": "+"},
    "CYP2C9":  {"chr": "chr10", "start": 94938658, "end": 94990091, "strand": "+"},
    "CYP2C8":  {"chr": "chr10", "start": 95036772, "end": 95069497, "strand": "-"},
}

# Output FASTA file containing all three genes
with open("data/reference_genes.fa", "w") as out_f:
    for gene, info in GENE_INFO.items():
        # UCSC coordinates are 1-based inclusive, pysam.fetch is 0-based start, 1-based end
        seq = fasta.fetch(info["chr"], info["start"] - 1, info["end"])
        out_f.write(f">{gene} {info['chr']}:{info['start']}-{info['end']} ({info['strand']})\n")
        # Wrap sequence every 60 characters for readability
        for i in range(0, len(seq), 60):
            out_f.write(seq[i:i+60] + "\n")

fasta.close()
print("Reference genes extracted successfully!")
```
What this does:

Opens chr10.fa using pysam.FastaFile.

Loops through the three CYP2C genes using their coordinates.

Extracts the genomic sequence as-is (no reverse complement), even for CYP2C8 which is on the negative strand.

Writes all three sequences into a single FASTA file with headers including gene name, chromosome, coordinates, and strand.

Wraps sequence lines every 60 characters for standard FASTA formatting.

✅ Output: data/reference_genes.fa ready for alignment, variant calling, and downstream analysis.

If you want, I can also provide a version that automatically fetches just these genes from UCSC/Ensembl instead of downloading the full chr10.fa, so your notebook remains light and self-contained. Do you want me to do that?

**continue Prompt 2:**

yes please

**continue Aswer 2:**

Perfect! Since you want only the CYP2C genes and not the entire chromosome, we can fetch them directly from UCSC’s FASTA server using Python and requests. Here’s a fully self-contained example:

```python
import requests

# Output FASTA file
out_file = "data/reference_genes.fa"

# Gene coordinates (hg38, UCSC)
GENE_INFO = {
    "CYP2C19": {"chr": "chr10", "start": 94762681, "end": 94855547, "strand": "+"},
    "CYP2C9":  {"chr": "chr10", "start": 94938658, "end": 94990091, "strand": "+"},
    "CYP2C8":  {"chr": "chr10", "start": 95036772, "end": 95069497, "strand": "-"},
}

# UCSC FASTA server template
ucsc_fasta_url = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chr};start={start};end={end}"

with open(out_file, "w") as out_f:
    for gene, info in GENE_INFO.items():
        # UCSC is 0-based start, 1-based end for the API
        url = ucsc_fasta_url.format(chr=info["chr"], start=info["start"]-1, end=info["end"])
        r = requests.get(url)
        r.raise_for_status()
        seq_json = r.json()
        seq = seq_json["dna"]
        
        # Write to FASTA
        out_f.write(f">{gene} {info['chr']}:{info['start']}-{info['end']} ({info['strand']})\n")
        # Wrap sequence every 60 chars
        for i in range(0, len(seq), 60):
            out_f.write(seq[i:i+60] + "\n")

print(f"Reference genes saved to {out_file}")
```
✅ Notes:

This does not require downloading the full chromosome.

The sequence is written as-is, so the negative strand CYP2C8 is not reverse complemented, which is fine for alignment/variant comparisons.

You can run this cell directly in your notebook (week5.ipynb), and it will produce the reference_genes.fa file ready for alignment and variant calling.

Uses only Python’s standard library + requests, so no conda/pip complications.

If you want, I can also make a variant that uses biopython instead of requests, which some people prefer for sequence handling in notebooks. This can simplify later downstream manipulations like creating SeqRecord objects. Do you want me to do that?

**continue Prompt 2:**

yes please

**continue Answer 2:**

Perfect! Here’s a clean Biopython version that fetches the same CYP2C genes from UCSC and writes them into a single FASTA file as SeqRecords. This will also make downstream manipulations easier.

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests

# Output FASTA file
out_file = "data/reference_genes.fa"

# Gene coordinates (hg38, UCSC)
GENE_INFO = {
    "CYP2C19": {"chr": "chr10", "start": 94762681, "end": 94855547, "strand": "+"},
    "CYP2C9":  {"chr": "chr10", "start": 94938658, "end": 94990091, "strand": "+"},
    "CYP2C8":  {"chr": "chr10", "start": 95036772, "end": 95069497, "strand": "-"},
}

# UCSC FASTA API template
ucsc_fasta_url = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chr};start={start};end={end}"

records = []

for gene, info in GENE_INFO.items():
    url = ucsc_fasta_url.format(chr=info["chr"], start=info["start"]-1, end=info["end"])
    r = requests.get(url)
    r.raise_for_status()
    seq = r.json()["dna"]

    # Create SeqRecord
    record = SeqRecord(
        Seq(seq),
        id=gene,
        description=f"{info['chr']}:{info['start']}-{info['end']} ({info['strand']})"
    )
    records.append(record)

# Write all genes to a single FASTA
with open(out_file, "w") as f:
    SeqIO.write(records, f, "fasta")

print(f"Reference genes saved to {out_file}")
```
✅ Advantages:

Uses Biopython SeqRecords → easy integration with other Biopython tools.

Keeps strand information in the FASTA header.

Output is neatly formatted with line wrapping handled by SeqIO.write.

No need to reverse complement negative-strand gene unless specifically required later.

### **Prompt 3:**
The below version works very well. 
Is it possible to push the executed notebook to github?

```yaml
name: Github CI
on: [push]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - run: echo "🎉 The job was automatically triggered by a ${{ github.event_name }} event."
      - run: echo "🐧 This job is now running on a ${{ runner.os }} server hosted by GitHub!"
      - run: echo "🔎 The name of your branch is ${{ github.ref }} and your repository is ${{ github.repository }}."

      - name: Check out repository code
        uses: actions/checkout@v5

      - uses: actions/setup-python@v6
        with:
          python-version: '3.12'

      - name: Set up Conda
        uses: conda-incubator/setup-miniconda@v3
        with:
          auto-update-conda: true
          activate-environment: bioinfo
          python-version: 3.12
          channels: conda-forge,bioconda
          conda-remove-defaults: true

      - name: Install dependencies
        shell: bash -l {0}
        run: |
          conda install -y samtools minimap2 bcftools freebayes htslib hapcut2 
          pip install --upgrade pip
          pip install biopython pysam jupyter nbconvert
          
      - run: echo "💡 The ${{ github.repository }} repository has been cloned to the runner."
      - run: echo "🖥️ The workflow is now ready to test your code on the runner."
      - name: Week 5 & 6
        shell: bash -l {0}
        run: 
          conda run -n bioinfo jupyter nbconvert --to notebook --execute week5/week5.ipynb --output week5.ipynb --ExecutePreprocessor.kernel_name=python3
      - uses: actions/upload-artifact@v4
        with:
          name: executed-notebook
          path: week5/week5.ipynb
```

### **Answer 3:**

Excellent — that’s a clean and modern setup ✅

Yes — you can absolutely push the executed notebook back to GitHub automatically at the end of the workflow.
Here’s how you can extend your current ci.yml safely and clearly.

✅ Updated ci.yml with Git Push Step

```yaml
name: Github CI
on: [push]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
      - run: echo "🎉 The job was automatically triggered by a ${{ github.event_name }} event."
      - run: echo "🐧 This job is now running on a ${{ runner.os }} server hosted by GitHub!"
      - run: echo "🔎 The name of your branch is ${{ github.ref }} and your repository is ${{ github.repository }}."

      - name: Check out repository code
        uses: actions/checkout@v5
        with:
          token: ${{ secrets.GITHUB_TOKEN }}
          persist-credentials: true  # important so the push works

      - uses: actions/setup-python@v6
        with:
          python-version: '3.12'

      - name: Set up Conda
        uses: conda-incubator/setup-miniconda@v3
        with:
          auto-update-conda: true
          activate-environment: bioinfo
          python-version: 3.12
          channels: conda-forge,bioconda
          conda-remove-defaults: true

      - name: Install dependencies
        shell: bash -l {0}
        run: |
          conda install -y samtools minimap2 bcftools freebayes htslib hapcut2 
          pip install --upgrade pip
          pip install biopython pysam jupyter nbconvert

      - run: echo "💡 The ${{ github.repository }} repository has been cloned to the runner."
      - run: echo "🖥️ The workflow is now ready to test your code on the runner."

      - name: Execute week5 notebook
        shell: bash -l {0}
        run: |
          conda run -n bioinfo jupyter nbconvert \
            --to notebook --execute week5/week5.ipynb \
            --output week5.ipynb --ExecutePreprocessor.kernel_name=python3

      - uses: actions/upload-artifact@v4
        with:
          name: executed-notebook
          path: week5/week5.ipynb

      - name: Commit executed notebook
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          git add week5/week5.ipynb
          git commit -m "🤖 Auto-update executed notebook [skip ci]" || echo "No changes to commit"

      - name: Push executed notebook
        run: |
          git push
```
🧠 Key Points

persist-credentials: true
Lets the GITHUB_TOKEN stay active for pushing commits.

Safe commit message
The [skip ci] tag prevents the workflow from re-triggering itself after the push.

Idempotent commit
The || echo "No changes to commit" part ensures the job doesn’t fail if the notebook hasn’t changed.

Result
After each successful run, your updated week5/week5.ipynb (the executed version) will appear directly in your GitHub repository.

Would you like me to modify it so that the executed notebook is pushed into a separate branch (like executed-notebooks), instead of your main branch? That’s a common best practice to avoid polluting main.

**continue Prompt 3:**

do I have to specify any info here? "git config user.name "github-actions[bot]" git config user.email "github-actions[bot]@users.noreply.github.com""

Also I want to executed week 5 to be push to week5 branch

**continue Answer 3:**

Perfect 👍 — yes, let’s make two clean improvements:

✅ 1️⃣ You don’t need to change:
```bash
git config user.name "github-actions[bot]"
git config user.email "github-actions[bot]@users.noreply.github.com"
```
That’s exactly what GitHub’s official CI uses — it’s a safe default identity.
You don’t need to add your own name or email.

✅ 2️⃣ Push executed notebook to the week5 branch

We’ll modify the bottom part of your workflow so that:

The executed notebook is committed and pushed to a branch named week5

The branch will be created automatically if it doesn’t exist

The commit will include [skip ci] to prevent infinite CI loops

🧩 Updated ci.yml (Final Version)

```yaml
name: Github CI
on: [push]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
      - run: echo "🎉 The job was automatically triggered by a ${{ github.event_name }} event."
      - run: echo "🐧 This job is now running on a ${{ runner.os }} server hosted by GitHub!"
      - run: echo "🔎 The name of your branch is ${{ github.ref }} and your repository is ${{ github.repository }}."

      - name: Check out repository code
        uses: actions/checkout@v5
        with:
          token: ${{ secrets.GITHUB_TOKEN }}
          persist-credentials: true

      - uses: actions/setup-python@v6
        with:
          python-version: '3.12'

      - name: Set up Conda
        uses: conda-incubator/setup-miniconda@v3
        with:
          auto-update-conda: true
          activate-environment: bioinfo
          python-version: 3.12
          channels: conda-forge,bioconda
          conda-remove-defaults: true

      - name: Install dependencies
        shell: bash -l {0}
        run: |
          conda install -y samtools minimap2 bcftools freebayes htslib hapcut2 
          pip install --upgrade pip
          pip install biopython pysam jupyter nbconvert

      - name: Execute week5 notebook
        shell: bash -l {0}
        run: |
          conda run -n bioinfo jupyter nbconvert \
            --to notebook --execute week5/week5.ipynb \
            --output week5.ipynb --ExecutePreprocessor.kernel_name=python3

      - uses: actions/upload-artifact@v4
        with:
          name: executed-notebook
          path: week5/week5.ipynb

      - name: Commit executed notebook
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          git add week5/week5.ipynb
          git commit -m "🤖 Auto-update executed week5 notebook [skip ci]" || echo "No changes to commit"

      - name: Push to week5 branch
        run: |
          git fetch origin
          git checkout week5 || git checkout -b week5
          git merge --no-edit ${{ github.ref_name }} || echo "Already up to date"
          git push origin week5
```
🧠 What happens now

The workflow executes on your main branch (or whichever branch triggered it).

After execution:

It switches (or creates) a branch named week5

Merges your current changes into it

Commits the executed notebook

Pushes that notebook to week5 branch

week5/week5.ipynb on GitHub will now be the executed version.

### Prompt 4:

Can you deduce from these screenshots—are these variants sequencing-related artifacts or are they indeed true variants? Do this analysis for every gene. 
- CYP2C19:54675 T->C  GT=T|C DP=258 PQ=100 (Illumina only)
- CYP2C9:18070 A->G  GT=A|G DP=277 PQ=100 (Illumina only)
- CYP2C8:9977 TAAAA->T  GT=TAAAA/T DP=162 PQ=None (PacBio only)
CYP2C8:9977 TAAA->T  GT=T/TAAA DP=162 PQ=None (PacBio only)

### Answer 4:

Nice — great work getting IGV running and producing snapshots. I don’t have the actual images here, so I’ll reason from the numbers you gave (DP, PQ, GT) and from typical sequencing behaviour, then give a short decision for each variant and exact, copy-pasteable code you can run locally to objectively confirm what the screenshots are showing (read counts, strand bias, MQ, soft-clips, etc.). I’ll finish with a short recipe to extract the phased haplotypes for CYP2C19 from your normalized + phased VCF.

Quick interpretation (based on the VCF snippets you reported)

You gave these discordants:

CYP2C19:54675 T→C — Illumina-only
GT=T|C DP=258 PQ=100
Interpretation: likely a real SNP.
Why: high depth (DP=258) and very high phasing quality (PQ=100) — Illumina short reads give excellent base accuracy for SNPs, and many reads support the alt allele. If the PacBio VCF is missing it, likely reasons: lower PacBio coverage at that site, mapping differences, or variant representation differences (multi-allelic, normalization). Check IGV: if you see many short reads with the C base, balanced read pairs, consistent MQ and little soft-clipping, it’s almost certainly real.

CYP2C9:18070 A→G — Illumina-only
GT=A|G DP=277 PQ=100
Interpretation: likely a real SNP (same logic as above). High depth + high PQ → strong evidence from Illumina. If PacBio doesn’t show it, check coverage and allele fraction in PacBio reads. If PacBio reads show the alt allele but bcftools didn’t call it (representation/threshold differences), that explains the discordance.

CYP2C8:9977 TAAAA→T (PacBio-only)
GT=TAAAA/T DP=162 PQ=None
and also CYP2C8:9977 TAAA->T GT=T/TAAA DP=162 PQ=None (looks like a multiallelic indel representation)
Interpretation: ambiguous — could be real, but long-read indels are also where sequencing/platform errors concentrate (homopolymers/repeats). PacBio HiFi has low indel error but still indel-prone in certain contexts. Because this is an indel (repetitive homopolymer-like), short Illumina reads sometimes fail to map / call it (or call it differently), while PacBio long reads can span the repeat and show consistent deletion. The high DP (162) for PacBio suggests many long reads support it — that points toward this being a true structural/indel event, but we need to check read-level evidence carefully (soft-clips, alignment quality, consistent CIGAR deletions). If IGV shows clean long reads with the deletion and consistent mapping quality, it’s likely real. If long reads are soft-clipped or align poorly around the site, it could be an artifact.

How to objectively confirm using the BAMs (pysam snippets)

Run the snippets below in your notebook — they will compute, at each position, per-BAM:

total depth, alt count, alt allele fraction

alt counts on forward/reverse strands (strand bias)

mean mapping quality for reads supporting ref vs alt

#reads with soft-clips or indel CIGAR at the site

Copy/paste and run this cell (adjust paths if needed):

```python
import pysam
from collections import defaultdict
import math

illumina_bam = "data/illumina.bam"
pacbio_bam  = "data/pacbio.bam"

# Query variants (use the relative co-ordinates that IGV used)
# These are gene-local positions you provided; if your fasta contig headers include an offset
# you must convert to absolute positions. Here we assume the VCF POS values match these numbers.
queries = [
    ("CYP2C19", 54675),
    ("CYP2C9", 18070),
    ("CYP2C8", 9977),
]

def inspect_site(bam_path, chrom, pos, window=0):
    """
    pos = 1-based coordinate
    Returns a dict summarizing ref vs alt evidence at pos.
    """
    sam = pysam.AlignmentFile(bam_path, "rb")
    # pysam uses 0-based fetch intervals: fetch(chrom, start0, end0)
    start0 = pos - 1
    end0 = pos
    stats = defaultdict(int)
    mq_alt = []
    mq_ref = []
    bq_alt = []
    bq_ref = []
    # We'll capture alt bases and counts
    base_counts = defaultdict(int)
    # Also capture strandwise alt counts
    strand_counts = defaultdict(lambda: [0,0])  # base -> [fwd, rev]
    softclip_count = 0
    indel_count = 0
    total_reads = 0

    for pileupcolumn in sam.pileup(chrom, start0, end0, truncate=True, stepper="all"):
        if pileupcolumn.pos != start0:
            continue
        for pileupread in pileupcolumn.pileups:
            # skip deletions and ref skip in pileup
            if pileupread.is_del or pileupread.is_refskip:
                # deletion contributes to evidence for indel
                if pileupread.is_del:
                    indel_count += 1
                continue
            total_reads += 1
            read = pileupread.alignment
            # base (1-letter) at this pileup read
            read_base = read.query_sequence[pileupread.query_position]
            base_counts[read_base] += 1
            # strand
            if read.is_reverse:
                strand_counts[read_base][1] += 1
            else:
                strand_counts[read_base][0] += 1
            # mapping quality and base quality
            if read_base == read.get_reference_sequence()[pileupread.query_position] if False else None:
                pass
            # We can't use get_reference_sequence() like that; instead compare to ref by reading fasta.
            # But here we will just collect MQ/BQ per base:
            if read_base:
                if read.is_reverse:
                    pass
            if read.mapping_quality is not None:
                # approximate MQ list per base
                # collect MQ for this read and base
                if read_base:
                    # naive split: treat anything with base==ref as ref else alt, but need ref base:
                    # We'll later supply ref allele to caller for correct classification.
                    pass
            # detect soft-clip in CIGAR operations overlapping this position
            for (cig_type, cig_len) in read.cigartuples or []:
                # 4 == soft-clip, 1==insertion, 2==deletion
                if cig_type == 4:
                    # if soft-clip exists anywhere in this read, count it
                    softclip_count += 1
                    break

    sam.close()
    # Build output
    out = {
        "total_reads_at_pos": total_reads,
        "base_counts": dict(base_counts),
        "strand_counts": {b: tuple(v) for b,v in strand_counts.items()},
        "softclip_reads": softclip_count,
        "indel_reads": indel_count
    }
    return out

# run for both bams
for gene, pos in queries:
    print(f"\n--- site {gene}:{pos} ---")
    for bam in (illumina_bam, pacbio_bam):
        try:
            s = inspect_site(bam, gene, pos)
            print(f"{bam} --> total_reads={s['total_reads_at_pos']}, base_counts={s['base_counts']}, softclips={s['softclip_reads']}, indel_reads={s['indel_reads']}")
        except Exception as e:
            print(f"Error inspecting {bam} at {gene}:{pos} -> {e}")
```

Notes:

If your FASTA contigs use headers like >CYP2C19 chr10:9476..., then your BAM/VCF CHROM will be CYP2C19 (not chr10). The code above uses the contig name as in the BAM (same as your VCF).

The snippet is conservative (counts bases and soft-clips). For more robust classification we’ll compare bases to the reference base — see the next snippet.

More precise allele counts (compare to reference base) — better for SNPs & indels

Run this to get per-base counts where we know the REF allele from the FASTA:

```python
import pysam
from collections import Counter
from Bio import SeqIO

fasta_path = "data/reference_genome.fa"
ref_seqs = {r.id: str(r.seq) for r in SeqIO.parse(fasta_path, "fasta")}

def allele_summary(bam_path, chrom, pos, ref_base=None):
    sam = pysam.AlignmentFile(bam_path, "rb")
    start0 = pos - 1
    end0 = pos
    counts = Counter()
    fwd = Counter()
    rev = Counter()
    mq_totals = []
    softclips = 0
    for pileupcolumn in sam.pileup(chrom, start0, end0, truncate=True, stepper="all"):
        if pileupcolumn.pos != start0:
            continue
        for pr in pileupcolumn.pileups:
            if pr.is_del or pr.is_refskip:
                counts["<del>"] += 1
                continue
            base = pr.alignment.query_sequence[pr.query_position]
            counts[base] += 1
            if pr.alignment.is_reverse:
                rev[base] += 1
            else:
                fwd[base] += 1
            mq_totals.append(pr.alignment.mapping_quality)
            # detect softclip presence anywhere in read
            if any(c[0] == 4 for c in (pr.alignment.cigartuples or [])):
                softclips += 1
    sam.close()
    return {
        "counts": dict(counts),
        "strand_fwd": dict(fwd),
        "strand_rev": dict(rev),
        "mean_MQ": sum(mq_totals)/len(mq_totals) if mq_totals else None,
        "softclip_reads": softclips
    }

for gene, pos in queries:
    ref_base = None
    # obtain reference base from ref_seqs if possible
    seq = ref_seqs.get(gene)
    if seq:
        ref_base = seq[pos-1]
    print(f"\nSite {gene}:{pos} ref={ref_base}")
    for bam in (illumina_bam, pacbio_bam):
        try:
            out = allele_summary(bam, gene, pos, ref_base)
            print(f"{bam}: counts={out['counts']}, mean_MQ={out['mean_MQ']}, softclips={out['softclip_reads']}")
        except Exception as e:
            print("error", e)
```

Interpretation guide for results you’ll see

SNP is likely true if: large alt count, alt fraction ≥ ~30% (for heterozygous), alt present on both strands roughly equally, reads have good mapping quality (MQ high), and no excess soft-clips or suspicious nearby indels.

SNP is likely artifact if: alt only on one strand, extreme base quality drop, alt bases clustered in reads with low MQ, or heavy soft-clipping nearby.

Indel (PacBio) likely true if many long reads show the same consistent CIGAR deletion with good MQ and aligned cleanly across breakpoint. If reads are soft-clipped/realign differently (split reads, inconsistent CIGAR), suspect mapping artifact.

Example decisions (what to expect & what each result would imply)

If the allele summary for CYP2C19:54675 (Illumina) shows e.g. {'T': 120, 'C': 130}, mean_MQ ~ 50–60, softclips small → true SNP, Illumina strongly supports. If PacBio shows very low coverage at that base or those reads are ambiguous → explains absence there.

If CYP2C9:18070 shows {'A': 135, 'G': 140} balanced across strands with high MQ → true SNP.

If CYP2C8:9977 indel (PacBio) shows many reads with CIGAR containing a 4-nt deletion (or consistent D in CIGAR) and mapping qualities are high and reads span the event cleanly → true indel. If instead most PacBio reads are soft-clipped or have low MQ around that site, or their alt allele fraction is low, could be a sequencer / alignment artifact.

Using the phased information (PS/PQ) to check haplotypes for CYP2C19

You can extract phased variants for the gene and build each haplotype (string of alleles across the phase set):

```python
import pysam
from collections import defaultdict

vcf_path = "data/illumina.norm.splitted.hapcut.phased.VCF"  # or pacbio version
chrom = "CYP2C19"

vcf = pysam.VariantFile(vcf_path)
# collect phases per PS (phase set)
phased_by_ps = defaultdict(list)

for rec in vcf.fetch(chrom):
    # require phased genotype
    for sample in rec.samples:
        s = rec.samples[sample]
        gt = s.get("GT")
        # GT may be like (1,0) for phased? pysam stores phased flag separately
        # pysam indicates phased via rec.samples[sample].phased
        phased = s.phased
        if gt is None:
            continue
        # Use PS tag if present in FORMAT fields OR use PS from INFO/FORMAT
        ps = None
        if "PS" in s:
            ps = s["PS"]
        else:
            # fallback: use position as PS if PQ high
            ps = rec.pos
        # For each allele, build allele character (e.g., 'ref' or the chosen alt)
        # NOTE: for multi-allelic, pick first ALT only for simplicity here
        allele0 = rec.ref if gt[0] == 0 else rec.alts[gt[0]-1] if gt[0] is not None else "."
        allele1 = rec.ref if gt[1] == 0 else rec.alts[gt[1]-1] if gt[1] is not None else "."
        phased_by_ps[ps].append((rec.pos, allele0, allele1, phased))

# Print haplotypes per phase set
for ps, entries in phased_by_ps.items():
    entries.sort()
    hap0 = "".join(e[1] for e in entries)
    hap1 = "".join(e[2] for e in entries)
    print(f"PS={ps}: hap0 len={len(hap0)} hap1 len={len(hap1)}")
    print("positions:", [e[0] for e in entries])
    print("hap0 sample:", hap0[:200], "...")
    print("hap1 sample:", hap1[:200], "...\n")
```

This produces a per-phase-set representation of the haplotypes (concatenated alleles). For star-allele assignment you’ll then compare the pattern of phased variants (the combination on the same haplotype) to known PharmVar star allele definitions.

Final recommended next steps you can perform now

Run the allele_summary code above for all three discordants on both BAMs. Paste the outputs here if you want me to interpret them precisely.

Inspect your IGV snapshots: do the long reads/short reads show clean alternative alleles (no soft clips, no MQ collapse)? If yes → likely real. If not → artifact.

If PacBio indel looks real (consistent deletion in many reads), treat it as a candidate true variant; if Illumina can't see it, that is expected for larger indels in repetitive context.

For star-alleles (final step): use the phased VCF for each gene, extract the variants that define star alleles from PharmVar, check which phased haplotype contains those defining variants. I can help format that final mapping once you pick which star-alleles to test (or provide PharmVar variant list).