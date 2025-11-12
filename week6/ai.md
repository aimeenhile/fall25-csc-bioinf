# Prompt 1:

I have finish setting up the data. Now, i need help with the next step. Please explain the code line by line for me if possible please.

2. Use Alevin-fry to align and quantify this data to the reference genome. See here for the installation steps.
Simpleaf is designed to simplify the alevin-fry interface for single-cell and nucleus raw data processing. It encapsulates the whole processing pipeline into two steps:

    simpleaf index indexes the provided reference or makes a splici reference (spliced transcripts + introns) and index it.

    simpleaf quant maps the sequencing reads against the indexed reference and quantifies the mapping records to generate a gene count matrix.

More advanced usages and options for mapping with simpleaf can be found here.

When running simpleaf index, if a genome FASTA file (-f) and a gene annotation GTF file(-g) are provided, it will generate a splici reference and index it; if only a transcriptome FASTA file is provided (--refseq), it will directly index it. Currently, we recommend the splici index.

```bash
# simpleaf needs the environment variable ALEVIN_FRY_HOME to store configuration and data.
# For example, the paths to the underlying programs it uses and the CB permit list
mkdir alevin_fry_home && export ALEVIN_FRY_HOME='alevin_fry_home'

# the simpleaf set-paths command finds the path to the required tools and writes a configuration JSON file in the ALEVIN_FRY_HOME folder.
simpleaf set-paths

# simpleaf index
# Usage: simpleaf index -o out_dir [-f genome_fasta -g gene_annotation_GTF|--refseq transcriptome_fasta] -r read_length -t number_of_threads
## The -r read_lengh is the number of sequencing cycles performed by the sequencer to generate biological reads (read2 in Illumina).
## Publicly available datasets usually have the read length in the description. Sometimes they are called the number of cycles.
simpleaf index \
-o simpleaf_index \
-f toy_human_ref/fasta/genome.fa \
-g toy_human_ref/genes/genes.gtf \
-r 90 \
-t 8
```

In the output directory simpleaf_index, the ref folder contains the splici reference; The index folder contains the salmon index built upon the splici reference.

The next step, simpleaf quant, consumes an index directory and the mapping record FASTQ files to generate a gene count matrix. This command encapsulates all the major steps discussed in this section, including mapping, cell barcode correction, and UMI resolution.

```bash
# Collecting sequencing read files
## The reads1 and reads2 variables are defined by finding the filenames with the pattern "_R1_" and "_R2_" from the toy_read_fastq directory.
reads1_pat="_R1_"
reads2_pat="_R2_"

## The read files must be sorted and separated by a comma.
### The find command finds the files in the fastq_dir with the name pattern
### The sort command sorts the file names
### The awk command and the paste command together convert the file names into a comma-separated string.
reads1="$(find -L ${fastq_dir} -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
reads2="$(find -L ${fastq_dir} -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"

# simpleaf quant
## Usage: simpleaf quant -c chemistry -t threads -1 reads1 -2 reads2 -i index -u [unspliced permit list] -r resolution -m t2g_3col -o output_dir
simpleaf quant \
-c 10xv3 -t 8 \
-1 $reads1 -2 $reads2 \
-i simpleaf_index/index \
-u -r cr-like \
-m simpleaf_index/index/t2g_3col.tsv \
-o simpleaf_quant
```

After running these commands, the resulting quantification information can be found in the simpleaf_quant/af_quant/alevin folder. Within this directory, there are three files: quants_mat.mtx, quants_mat_cols.txt, and quants_mat_rows.txt, which correspond, respectively, to the count matrix, the gene names for each column of this matrix, and the corrected, filtered cell barcodes for each row of this matrix. The tail lines of these files are shown below. Of note here is the fact that alevin-fry was run in the USA-mode (unspliced, spliced, and ambiguous mode), and so quantification was performed for both the spliced and unspliced status of each gene — the resulting quants_mat_cols.txt file will then have a number of rows equal to 3 times the number of annotated genes which correspond, to the names used for the spliced (S), unspliced (U), and splicing-ambiguous variants (A) of each gene.

```bash
# Each line in `quants_mat.mtx` represents
# a non-zero entry in the format row column entry
$ tail -3 simpleaf_quant/af_quant/alevin/quants_mat.mtx
138 58 1
139 9 1
139 37 1

# Each line in `quants_mat_cols.txt` is a splice status
# of a gene in the format (gene name)-(splice status)
$ tail -3 simpleaf_quant/af_quant/alevin/quants_mat_cols.txt
ENSG00000120705-A
ENSG00000198961-A
ENSG00000245526-A

# Each line in `quants_mat_rows.txt` is a corrected
# (and, potentially, filtered) cell barcode
$ tail -3 simpleaf_quant/af_quant/alevin/quants_mat_rows.txt
TTCGATTTCTGAATCG
TGCTCGTGTTCGAAGG
ACTGTGAAGAAATTGC
```

We can load the count matrix into Python as an AnnData object using the load_fry function from pyroe. A similar function, loadFry, has been implemented in the fishpond R package.

```python
import pyroe

quant_dir = 'simpleaf_quant/af_quant'
adata_sa = pyroe.load_fry(quant_dir)
```

The default behavior loads the X layer of the Anndata object as the sum of the spliced and ambiguous counts for each gene. However, recent work [Pool et al., 2022] and updated practices suggest that the inclusion of intronic counts, even in single-cell RNA-seq data, may increase sensitivity and benefit downstream analyses. While the best way to make use of this information is the subject of ongoing research, since alevin-fry automatically quantifies spliced, unspliced, and ambiguous reads in each sample, the count matrix containing the total counts for each gene can be simply obtained as follows:

```python
import pyroe

quant_dir = 'simpleaf_quant/af_quant'
adata_usa = pyroe.load_fry(quant_dir, output_format={'X' : ['U','S','A']})
```

3.8.3. The complete alevin-fry pipeline

Simpleaf makes it possible to process single-cell raw data in the “standard” way with a few commands. Next, we will show how to generate the identical quantification result by explicitly calling the pyroe, salmon, and alevin-fry commands. On top of the pedagogical value, knowing the exact command of each step will be helpful if only a part of the pipeline needs to be rerun or if some parameters not currently exposed by simpleaf need to be specified.

Please note that the commands in the Preparation section should be executed in advance. All the tools called in the following commands, pyroe, salmon, and alevin-fry, have already been installed when installing simpleaf.
3.8.3.1. Building the index

First, we process the genome FASTA file and gene annotation GTF file to obtain the splici index. The commands in the following code chunk are analogous to the simpleaf index command discussed above. This includes two steps:

    Building the splici reference (spliced transcripts + introns) by calling pyroe make-splici, using the genome and gene annotation file

    Indexing the splici reference by calling salmon index

```bash
# make splici reference
## Usage: pyroe make-splici genome_file gtf_file read_length out_dir
## The read_lengh is the number of sequencing cycles performed by the sequencer. Ask your technician if you are not sure about it.
## Publicly available datasets usually have the read length in the description.
pyroe make-splici \
${ref_dir}/fasta/genome.fa \
${ref_dir}/genes/genes.gtf \
90 \
splici_rl90_ref

# Index the reference
## Usage: salmon index -t extend_txome.fa -i idx_out_dir -p num_threads
## The $() expression runs the command inside and puts the output in place.
## Please ensure that only one file ends with ".fa" in the `splici_ref` folder.
salmon index \
-t $(ls splici_rl90_ref/*\.fa) \
-i salmon_index \
-p 8
```

The splici index can be found in the salmon_index directory.
3.8.3.2. Mapping and quantification

Next, we will map the sequencing reads recorded against the splici index by calling salmon alevin. This will produce an output folder called salmon_alevin that contains all the information we need to process the mapped reads using alevin-fry.

```bash
# Collect FASTQ files
## The filenames are sorted and separated by space.
reads1="$(find -L $fastq_dir -name "*$reads1_pat*" -type f | sort | awk '{$1=$1;print}' | paste -sd' ')"
reads2="$(find -L $fastq_dir -name "*$reads2_pat*" -type f | sort | awk '{$1=$1;print}' | paste -sd' ')"

# Mapping
## Usage: salmon alevin -i index_dir -l library_type -1 reads1_files -2 reads2_files -p num_threads -o output_dir
## The variable reads1 and reads2 defined above are passed in using ${}.
salmon alevin \
-i salmon_index \
-l ISR \
-1 ${reads1} \
-2 ${reads2} \
-p 8 \
-o salmon_alevin \
--chromiumV3 \
--sketch
```

Then, we execute the cell barcode correction and UMI resolution step using alevin-fry. This procedure involves three alevin-fry commands:

    The generate-permit-list command is used for cell barcode correction.

    The collate command filters out invalid mapping records, corrects cell barcodes and collates mapping records originating from the same corrected cell barcode.

    The quant command performs UMI resolution and quantification.

# Cell barcode correction
## Usage: alevin-fry generate-permit-list -u CB_permit_list -d expected_orientation -o gpl_out_dir
## Here, the reads that map to the reverse complement strand of transcripts are filtered out by specifying `-d fw`.
alevin-fry generate-permit-list \
-u 3M-february-2018.txt \
-d fw \
-i salmon_alevin \
-o alevin_fry_gpl

# Filter mapping information
## Usage: alevin-fry collate -i gpl_out_dir -r alevin_map_dir -t num_threads
alevin-fry collate \
-i alevin_fry_gpl \
-r salmon_alevin \
-t 8

# UMI resolution + quantification
## Usage: alevin-fry quant -r resolution -m txp_to_gene_mapping -i gpl_out_dir -o quant_out_dir -t num_threads
## The file ends with `3col.tsv` in the splici_ref folder will be passed to the -m argument.
## Please ensure that there is only one such file in the `splici_ref` folder.
alevin-fry quant -r cr-like \
-m $(ls splici_rl90_ref/*3col.tsv) \
-i alevin_fry_gpl \
-o alevin_fry_quant \
-t 8

After running these commands, the resulting quantification information can be found in alevin_fry_quant/alevin. Other relevant information concerning the mapping, CB correction, and UMI resolution steps can be found in the salmon_alevin, alevin_fry_gpl, and alevin_fry_quant folders, respectively.

In the example given here, we demonstrate using simpleaf and alevin-fry to process a 10x Chromium 3’ v3 dataset. Alevin-fry and simpleaf provide many other options for processing different single-cell protocols, including but not limited to Dropseq [Macosko et al., 2015], sci-RNA-seq3 [Cao et al., 2019] and other 10x Chromium platforms. A more comprehensive list and description of available options for different stages of processing can be found in the alevin-fry and simpleaf documentation. alevin-fry also provides a nextflow-based workflow, called quantaf, for conveniently processing many samples from a simply-defined sample sheet.

Of course, similar resources exist for many of the other raw data processing tools referenced and described throughout this section, including zUMIs [Parekh et al., 2018], alevin [Srivastava et al., 2019], kallisto|bustools [Melsted et al., 2021], STARsolo [Kaminow et al., 2021] and CellRanger. The scrnaseq pipeline from nf-core also provides a nextflow-based pipeline for processing single-cell RNA-se

3. Perform cell clustering (via Leiden modularity algorithm). Output the clustering plot.

4. Perform automatic cell annotation via CellTypist. Annotate the plot with the cell types.

# Answer 1:

Nice — you’ve got the whole Alevin-fry workflow laid out. I’ll explain each command / code block line-by-line (and give short notes on why it’s there and common pitfalls). I’ll follow the order in your message and mark small tips you’ll want in CI / Jupyter.

Explanation — step by step
1) Prepare Alevin-fry home (config + state)
```bash
# simpleaf needs the environment variable ALEVIN_FRY_HOME to store configuration and data.
# For example, the paths to the underlying programs it uses and the CB permit list
mkdir alevin_fry_home && export ALEVIN_FRY_HOME='alevin_fry_home'
```
mkdir alevin_fry_home — create a directory named alevin_fry_home.

&& — run the next command only if the mkdir succeeded.

export ALEVIN_FRY_HOME='alevin_fry_home' — set an environment variable ALEVIN_FRY_HOME to that directory for the current shell session.
Why: simpleaf uses this directory to store its JSON config and local metadata (e.g., paths to salmon/alevin-fry binaries, CB permit list).
Tip: in a script or CI, ensure this export runs in the same shell where commands run (or prefix later commands with ALEVIN_FRY_HOME=...).

2) Let simpleaf discover tool locations and write config

```bash
# the simpleaf set-paths command finds the path to the required tools and writes a configuration JSON file in the ALEVIN_FRY_HOME folder.
simpleaf set-paths
```
simpleaf set-paths — interrogates your PATH to find required underlying programs (salmon, alevin-fry, pyroe, etc.) and writes a config file (usually inside $ALEVIN_FRY_HOME) so simpleaf knows where to call them.
Why: avoids hard-coding binary locations later and allows simpleaf to use its bundled/installed tools.
Tip: if set-paths can’t find tools, you’ll see an error; ensure conda env or PATH is set correctly.

3) Build a splici index with simpleaf
```bash
simpleaf index \
  -o simpleaf_index \
  -f toy_human_ref/fasta/genome.fa \
  -g toy_human_ref/genes/genes.gtf \
  -r 90 \
  -t 8
```
Breakdown:

simpleaf index — build an index for mapping.

\ — line continuation in shell for readability.

-o simpleaf_index — output directory for index artifacts.

-f toy_human_ref/fasta/genome.fa — path to genome FASTA. When combined with -g (GTF), simpleaf will generate a splici reference (spliced transcripts + intronic sequences).

-g toy_human_ref/genes/genes.gtf — gene annotation GTF to derive transcripts and introns.

-r 90 — read length (number of cycles / bases in the biological read, typically read 2 for 10x). Important when building intronic extensions / handling fragment lengths.

-t 8 — number of CPU threads to use.

What it produces:

In simpleaf_index/ref — the splici FASTA (extended transcriptome containing spliced sequences + introns).

In simpleaf_index/index — a salmon index built on that splici FASTA.

Pitfalls / tips:

Use the splici index for many single-cell pipelines — it lets you capture intronic reads.

Make sure the GTF and FASTA use the same chromosome names and the GTF is compatible.

4) Collect FASTQ file lists (for paired reads)
```bash
# patterns
reads1_pat="_R1_"
reads2_pat="_R2_"

reads1="$(find -L ${fastq_dir} -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
reads2="$(find -L ${fastq_dir} -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
```
Breakdown for each command:

reads1_pat / reads2_pat — strings used to identify read1 / read2 files (common Illumina pattern: _R1_ and _R2_).

find -L ${fastq_dir} -name "*$reads1_pat*" -type f — find files under $fastq_dir matching * _R1_*. -L follows symlinks.

| sort — ensure deterministic order (important to pair files correctly).

| awk -v OFS=, '{$1=$1;print}' — format output using , as output field separator (awk idiom to normalize whitespace).

| paste -sd, — join all lines into a single line separated by commas.

reads1="$(...)" — capture the full comma-separated list into reads1 shell variable.

Why: simpleaf quant (and many mapping commands) accept a comma separated list of fastq files for -1 and -2 arguments. Sorting makes pairing deterministic.

CI/Jupyter note: in Jupyter ! cells, variable assignments using $(...) are executed in a subshell; if you need to use those variables across cells, export them or keep that logic in a single bash cell.

5) Run simpleaf quant (mapping + quantification)
```bash
simpleaf quant \
  -c 10xv3 -t 8 \
  -1 $reads1 -2 $reads2 \
  -i simpleaf_index/index \
  -u -r cr-like \
  -m simpleaf_index/index/t2g_3col.tsv \
  -o simpleaf_quant
```
Flags explained:

simpleaf quant — high-level command that runs mapping, barcode correction and UMI resolution end-to-end.

-c 10xv3 — chemistry: tells simpleaf/alevin how CB/UMI are laid out (10x Chromium v3 here).

-t 8 — threads.

-1 $reads1 -2 $reads2 — paired read files (comma-separated lists you built).

-i simpleaf_index/index — path to the salmon index created earlier.

-u — use a CB permit list? (Depending on simpleaf version, -u may indicate that a whitelist is provided or to use implicit permit list — check simpleaf --help; some variants accept -u <permit_file>).

-r cr-like — resolution mode for UMI resolution. cr-like or cr-like corresponds to 10x UMI strategies (various options: cr-like, cr-like-no-gene, sketch, etc.). This controls how to collapse UMIs into counts (alevin-fry resolution strategy).

-m simpleaf_index/index/t2g_3col.tsv — transcript→gene three-column mapping file (txp, gene, status) used to summarize transcript-level mappings to gene-level counts. t2g_3col.tsv is created when building splici.

-o simpleaf_quant — output directory.

What it does: maps reads (via salmon), corrects cell barcodes, performs UMI-resolution and writes quant outputs (in simpleaf_quant/af_quant/alevin/ with quants_mat.mtx etc).

Tip: confirm the -u semantics for your simpleaf version (sometimes it expects -u <permit_file>; in your earlier blocks you passed -u and separately -u 3M-february-2018.txt to alevin-fry directly).

6) Resulting quant files — what they are

quants_mat.mtx — sparse matrix in Matrix Market format (each line: row col value).

quants_mat_cols.txt — column labels for the matrix (here columns represent gene+splice status like ENSG...-S, ENSG...-U, ENSG...-A).

quants_mat_rows.txt — row labels (cell barcodes, corrected/filtered).

Note on USA mode: USA (Unspliced / Spliced / Ambiguous) means each gene has 3 entries (S/U/A) — so final columns = (genes × 3). pyroe.load_fry can combine these into useful X representations.

7) Load into Python (AnnData) with pyroe
```python
import pyroe

quant_dir = 'simpleaf_quant/af_quant'
adata_sa = pyroe.load_fry(quant_dir)
```
pyroe.load_fry(quant_dir) — reads the quants_mat.* files and builds an AnnData object.
Default: loads X as (by default) the sum of spliced+ambiguous? (your note says default loads X as sum of spliced + ambiguous — check local pyroe version for exact defaults).

adata_sa — an AnnData object with .X, .var (genes), .obs (cells) where .X is the expression matrix.

Alternate where you get full USA layers: `adata_usa = pyroe.load_fry(quant_dir, output_format={'X' : ['U','S','A']})`
output_format={'X' : ['U','S','A']} — tell loader to populate X using U/S/A separately (could be stored as separate layers or stacked, depending on pyroe implementation). This yields an object where you can choose to analyze spliced-only, or total (U+S+A), etc.

Why: having U/S/A allows RNA velocity-aware analyses or to retain intronic info for sensitivity.

8) (Pedagogical) Explicit commands that simpleaf hides

If you want to run the pipeline explicitly (useful for partial reruns, custom flags):

8a) Make splici reference with pyroe
```bash
pyroe make-splici \
  ${ref_dir}/fasta/genome.fa \
  ${ref_dir}/genes/genes.gtf \
  90 \
  splici_rl90_ref
```
pyroe make-splici <genome.fa> <genes.gtf> <read_length> <out_dir> — constructs splici FASTA (splice + intron sequences) suitable for indexing.

8b) Build salmon index
```bash
salmon index \
  -t $(ls splici_rl90_ref/*\.fa) \
  -i salmon_index \
  -p 8
```
$(ls splici_rl90_ref/*\.fa) — command substitution; ensures the .fa splici file is passed as -t argument. Assumes exactly one .fa in that folder.

-i salmon_index — output index directory.

-p 8 — threads.

8c) Map reads with salmon alevin
```bash
salmon alevin \
  -i salmon_index \
  -l ISR \
  -1 ${reads1} \
  -2 ${reads2} \
  -p 8 \
  -o salmon_alevin \
  --chromiumV3 \
  --sketch
```
salmon alevin — salmon’s single-cell subcommand.

-l ISR — library type (fragment orientation / layout). ISR commonly used for paired-end where second read is the biological read on reverse orientation — check your dataset.

--chromiumV3 — 10x Chromium v3 chemistry hint (affects how CB/UMI are parsed).

--sketch — build the lightweight mapping index used for sketching (speeds things). Check docs for exact behavior.

Output: salmon_alevin folder containing Alevin map data (barcode/UMI assignments + mapping records).

8d) Generate permit list (CB correction)
```bash
alevin-fry generate-permit-list \
  -u 3M-february-2018.txt \
  -d fw \
  -i salmon_alevin \
  -o alevin_fry_gpl
```
generate-permit-list — produce a list of candidate cell barcodes (permitted CBs) based on whitelist & mapping.

-u 3M-february-2018.txt — use this whitelist file (10x 3M whitelist).

-d fw — expected read orientation: fw indicates forward mapping orientation is expected; using -d can filter reads mapping to wrong strand.

-i salmon_alevin — where the Alevin mapping results are.

-o alevin_fry_gpl — output directory for permit list artifacts.

8e) Collate (filter + group mapping records by corrected CB)
```bash
alevin-fry collate \
  -i alevin_fry_gpl \
  -r salmon_alevin \
  -t 8
```
collate — organizes mapping records by corrected barcode (collapse records with same corrected CB), and filters invalid records; prepares data for quantification.

-r salmon_alevin — original mappings used as input.

-t 8 — threads.

8f) UMI resolution & quant
```bash
alevin-fry quant -r cr-like \
  -m $(ls splici_rl90_ref/*3col.tsv) \
  -i alevin_fry_gpl \
  -o alevin_fry_quant \
  -t 8
```
alevin-fry quant — take collated records and perform UMI-resolution (collapse UMIs according to the resolution algorithm) to produce counts.

-r cr-like — resolution strategy (10x-like).

-m <txp_to_gene_map> — 3-column transcript→gene map (splici created these).

-i alevin_fry_gpl — input from generate-permit-list / collate stage.

-o alevin_fry_quant — output dir with alevin/ containing quants_mat*.

-t 8 — threads.

Why do this explicitly?

It gives you control and transparency; simpleaf hides these steps behind simpleaf quant.

9) Downstream: clustering + annotation (brief practical notes)

You asked steps 3 and 4 too (Leiden clustering + CellTypist annotation). High-level pointers:

After pyroe.load_fry you’ll have an AnnData (adata). Typical single-cell flow (scanpy):
```python
import scanpy as sc
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=1.0)
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')
```
Leiden: sc.tl.leiden runs the Leiden modularity algorithm on the neighborhood graph (requires scanpy and leidenalg).

CellTypist: programmatic annotation:
```python
import celltypist
pred = celltypist.annotate(adata, model='Immune_All_Low.pkl')  # example
adata.obs['celltypist'] = pred.predicted_labels
sc.pl.umap(adata, color='celltypist', legend_loc='on data')
```
Choose the appropriate CellTypist model for your sample (immune-only vs pan-cell).

CellTypist returns predictions and confidence; you can use thresholding or majority-vote per cluster.

Practical tips / gotchas & CI/Jupyter notes

Environment persistence: export ALEVIN_FRY_HOME=... must be visible to commands — if using multiple notebook ! cells, either run the export and all commands in one bash cell, or prefix commands with ALEVIN_FRY_HOME=... simpleaf ....

Conda activation in notebooks: conda activate doesn’t persist across separate ! cells. Either (A) use a single bash cell for setup and run everything there, (B) rely on CI job steps to create/activate envs, or (C) prefix commands with conda run -n af <command>.

Single-file assumptions: commands like $(ls splici_rl90_ref/*\.fa) or $(ls splici_rl90_ref/*3col.tsv) assume exactly one matching file. If multiple files exist you’ll get errors — check with ls first.

Threads/Memory: salmon/alevin-fry are multi-threaded and memory heavy for large refs — you used -t 8, fine for small examples.

Whitelist semantics: make sure -u usage is correct for your simpleaf version (some commands expect -u <file>).

USA mode: if you use U/S/A, be deliberate about which counts you analyze (S only for mature mRNA, U for intronic, total for sensitivity).

File locations: outputs are in simpleaf_quant/af_quant/alevin/ or alevin_fry_quant/alevin/ depending on path; inspect metrics files if something looks off.

# Prompt 2:

my goal is to run this notebook automatically in github actions. How do i prepare alevin-fry home in jupyter notebook or actions ci to do that? Please note that I do not want to change the structure and the content of my current CI. I only want to add the alevin-fry home where possible.
```yaml
my .github/workflows/ci.yml:
name: Week 6 Deliverable
on: [push]
jobs:
  test:
    runs-on: ubuntu-latest

    permissions:
      contents: write

    steps:
      - run: echo "🎉 The job was automatically triggered by a ${{ github.event_name }} event."
      - run: echo "🐧 This job is now running on a ${{ runner.os }} server hosted by GitHub!"
      - run: echo "🔎 The name of your branch is ${{ github.ref }} and your repository is ${{ github.repository }}."

      - name: Check out repository code
        uses: actions/checkout@v5
        with:
          ref: ${{ github.head_ref }}

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
          conda install -c conda-forge openjdk=21
          conda install -y -c bioconda simpleaf biopython pysam samtools minimap2 bcftools htslib hapcut2 igv pandas
          pip install --upgrade pip
          pip install cyvcf2 igv-notebook jupyter nbconvert
          
      - run: echo "💡 The ${{ github.repository }} repository has been cloned to the runner."
      - run: echo "🖥️ The workflow is now ready to test your code on the runner."

      - name: Week 6
        shell: bash -l {0}
        run: 
          conda run -n bioinfo jupyter nbconvert --to notebook --execute week6/week6.ipynb  
```

# Answer 2: 

Perfect — this is a great and realistic automation question 🔬

You already have an excellent CI skeleton (it sets up Conda, installs simpleaf, runs Jupyter automatically, etc.).
To make alevin-fry and simpleaf work non-interactively inside GitHub Actions or Jupyter, you just need to ensure one thing:

✅ The environment variable ALEVIN_FRY_HOME exists and points to a writable directory before you call simpleaf.

Below I’ll show you exactly where and how to set that up — both in your notebook and in your CI YAML — without breaking or restructuring anything.

🧩 1. What ALEVIN_FRY_HOME is

simpleaf needs a directory where it can store:

its config file (config.json),

cached binary paths to salmon, alevin-fry, and pyroe,

and sometimes the permit list or temp files.

If the variable isn’t defined, simpleaf set-paths will fail or hang asking for user input.

So you must create a directory (e.g., ~/.alevin_fry_home or alevin_fry_home) and export it before running simpleaf.

🧠 2. Minimal setup — inside your Jupyter Notebook

Add this at the very top of your week6.ipynb, before any call to simpleaf.
```python
# --- Alevin-fry setup ---
import os, subprocess

# Define the home directory for alevin-fry (can be inside the notebook workspace)
os.makedirs("alevin_fry_home", exist_ok=True)
os.environ["ALEVIN_FRY_HOME"] = os.path.abspath("alevin_fry_home")

# Optional: verify
print("ALEVIN_FRY_HOME =", os.environ["ALEVIN_FRY_HOME"])

# Run once to initialize the configuration (non-interactive)
subprocess.run(["simpleaf", "set-paths"], check=True)
```

✅ This ensures that even when CI executes Jupyter automatically,
ALEVIN_FRY_HOME exists and simpleaf knows where to find its tools.

⚙️ 3. Setting it up in your GitHub Actions (recommended safety)

You can add just one line to your existing CI YAML under the “Week 6” step, so that the environment is already prepared before the notebook runs.

Add this:
```yaml
      - name: Week 6
        shell: bash -l {0}
        env:
          ALEVIN_FRY_HOME: ${{ github.workspace }}/alevin_fry_home
        run: |
          mkdir -p $ALEVIN_FRY_HOME
          conda run -n bioinfo jupyter nbconvert --to notebook --execute week6/week6.ipynb
```
Explanation:

| Line                        | Meaning                                                                        |
| --------------------------- | ------------------------------------------------------------------------------ |
| `env:`                      | Sets environment variables visible to the notebook and `simpleaf`              |
| `ALEVIN_FRY_HOME:`          | Creates a per-run home directory for alevin-fry in your repository’s workspace |
| `mkdir -p $ALEVIN_FRY_HOME` | Makes sure the folder exists before executing the notebook                     |
| `conda run ...`             | Runs your existing notebook execution, unchanged                               |

✅ This keeps your CI structure intact.
✅ You do not need to modify Conda setup or install dependencies sections.

🧰 4. (Optional) Testing locally

You can test this exact setup in your terminal or inside a Jupyter cell:
```bash
export ALEVIN_FRY_HOME=$(pwd)/alevin_fry_home
mkdir -p $ALEVIN_FRY_HOME
simpleaf set-paths
```
Then check:
```bash
ls $ALEVIN_FRY_HOME
# should contain config.json and tool paths
```
✅ Summary

| Where                       | Add                                                             | Purpose                                               |
| --------------------------- | --------------------------------------------------------------- | ----------------------------------------------------- |
| **Notebook (top)**          | Python cell creating and exporting `ALEVIN_FRY_HOME`            | Works in local + CI execution                         |
| **CI YAML (`Week 6` step)** | `env: ALEVIN_FRY_HOME: ${{ github.workspace }}/alevin_fry_home` | Ensures environment var is defined when notebook runs |
| **No structural change**    | You keep all other steps exactly as before                      | Seamless addition                                     |

