# SPARK: simulated pre-mRNA and RNA kinetics
SPARK (Simulated Pre-mRNA and RNA Kinetics) is a comprehensive simulation framework for generating nascent RNA benchmarking datasets. SPARK uses real annotated gene sequences to simulate nascent RNA sequencing reads after performing in silico transcription, co-transcriptional RNA processing (splicing, cleavage and polyadenylation), nascent RNA enrichment, and library preparation steps.

In brief, these simulations start from annotated mRNA isoforms and assign transcription elongation rates to each nucleotide across a chosen number of genes. Simulated RNA Polymerase II (RNAPII) is then walked across each gene for a user-defined labeling time, producing nascent pre-mRNAs that can be co- or post-transcriptionally spliced, cleaved, and polyadenylated. User-defined simulation parameters can be tuned to vary across a range of experimental conditions, and SPARK can reproduce a wide variety of nascent RNA sequencing approaches. During in silico RNA synthesis, users can choose methods that either mimic metabolic labeling or labeling of RNAPII positions: to simulate metabolic labeling approaches like 4sU-seq or TT-seq, SPARK simulates the stochastic incorporation of nucleotide analogues (e.g., 4sU, 5-ethynyluridine (5eU)), while the incorporation of a "biotin-NTP" recreates the process of transcription termination inherent to run-on assays like PRO-seq. SPARK can also retain the chromatin-associated RNA fraction to mimic chromatin-associated RNA-seq (caRNA-seq).

Gene expression in SPARK is not assigned directly: each gene's read count follows from its simulated transcription kinetics (initiation rate, elongation rates, pausing and RNA processing) and from what the chosen assay captures, so the same gene has a different expression level in each assay, as in real data.

[![](https://github.com/ezecalvo/spark/raw/master/figures/SPARK_summary.png)](/ezecalvo/spark/blob/master/figures/SPARK_summary.png)
### Requirements (versions used for development)
- python (v3.12.9)
- bedtools (v2.31.1)
- R (v4.3.3)
- bash (v4.4.20)
- awk (v4.2.1)
#### Python Dependencies
- numpy (v2.2.3)
- pandas (v2.2.3)
- scipy (v1.15.2)
- tqdm (v4.67.1)
- pysam (v0.23.0)
#### R Dependencies
- data.table (v1.17.0)
- dplyr (v1.1.4)
- optparse (v1.7.5)
- R.utils (v2.13.0)
- tidyr (v1.3.1)
- magrittr (v2.0.3)
- scales (v1.3.0)

A spark.yml conda environment is provided and can be installed by running
```
conda env create -f spark.yml
```
`spark.yml` pins Linux (linux-64) builds. SPARK also runs on macOS; there, create an environment with the packages listed above.
## Table of Contents


[Quick start with presets](#quick-start-with-presets)


[Overview of SPARK (full pipeline)](#overview-of-spark)


[1. Gene Selection](#step-1-gene-selection)


[2. Elongation Rates](#step-2-elongation-rates)


[3. Nascent RNA generation](#step-3-nascent-rna-generation)


[4. RNA processing (splicing, cleavage & polyadenylation)](#step-4-rna-processing)


[5. Enrichment](#step-5-enrichment)


- [Metabolic labeling RNA enrichment](#metabolic-labeling-rna-enrichment)
- [Transient transcriptome (TT)-seq](#transient-transcriptome-tt-seq)
- [Precision Nuclear Run-on (PRO)-seq](#precision-nuclear-run-on-pro-seq)
- [Native Elongating Transcript (NET)-seq](#native-elongating-transcript-net-seq)
- [Chromatin-associated RNA](#chromatin-associated-rna)
- [Sequencing errors and background molecules](#sequencing-errors-and-background-molecules)


[6. Read Generation](#step-6-read-generation)


- [Sequencing depth and gene expression](#sequencing-depth-and-gene-expression)


[Ground truth](#ground-truth)


- [Consolidated per-gene ground truth](#consolidated-per-gene-ground-truth)
- [Elongation rates & pausing](#elongation-rates--pausing)
- [Splicing ground truth](#splicing-ground-truth)
- [Nucleoside analog incorporation & substitution](#nucleoside-analog-incorporation--substitution)
- [Fragmentation and size selection](#fragmentation-and-size-selection)
- [For read mapping](#for-read-mapping)


## Quick start with presets


SPARK has a set of presets that load a predefined combination of parameters that recreate common nascent RNA sequencing methods. A preset only sets the parameters that define that assay (e.g., experiment type, labeling base, incorporation probability, background fraction), and any argument value can still be overridden by passing the corresponding flag.
Load a preset with `--preset` and list all available presets and their values with `--show_presets`:
```
spark.py --show_presets
```


### Available presets
| Preset | Assay it mimics | Experiment type | Labeling | Conversion | nt incorporation prob. | Background |
| --- | --- | --- | --- | --- | --- | --- |
| `4sUseq` | 4sU-seq RNA pulldown | metaboliclabeling | T | none | 0.005,0.023 | 0.15 |
| `DRB4sUseq` | 4sUDRB-seq | metaboliclabeling + `--drb` | T | none | 0.005,0.023 | 0.15 |
| `BrUseq` | BrU-seq / Bru-seq | metaboliclabeling | T | none | 0.001,0.005 | 0.15 |
| `TTseq` | Transient transcriptome (TT)-seq | ttseq | T | none | 0.005,0.023 | 0.10 |
| `TT-TimelapseSeq` | TT-TimeLapse-seq (T→C recoding) | ttseq | T | C | 0.005,0.023 | 0.05 |
| `eSLAM-seq` | eSLAM-seq / enriched TimeLapse-seq (T→C recoding) | metaboliclabeling | T | C | 0.005,0.023 | 0.10 |
| `proseq` | Precision Nuclear Run-on (PRO)-seq | proseq | A,C,G,T | — | 0.005,0.023 (default) | 0.02 |
| `NETseq` | Native Elongating Transcript (NET)-seq | netseq | — | — | — | 0.10 |
| `chromatinassociated` | Chromatin-associated RNA-seq | chromatin_associated | — | — | — | 0.20 |
| `polyA` | Steady-state polyA RNA-seq | metaboliclabeling | — | none | — | 1 (all reads) |

The Background column is the fraction of each gene's reads that come from background (non-nascent) molecules (see [Sequencing errors and background molecules](#sequencing-errors-and-background-molecules)).

Presets can also be customized by additional flags on top of it. For example, to run the NETseq preset but with a deeper sequencing library:
```
spark.py --preset NETseq -o ./out_dir/ --n 1000 \
--gtf GRCh38.95.gtf --genome_fasta GRCh38.fa \
--seq_type PE --experiment_time 30 --seq_depth 40000000 --seed 1
```
`--seq_depth` is the total number of reads (read pairs for paired-end) in the final library. The final simulated libraries will be located in the `out_dir/final_libraries/` directory. One fastq file is generated per read (e.g., R1 and R2 for paired-end). The naming convention is `{seq_tech}_{experiment_type}_{seq_type}_{date}_{read}.fastq.gz`. For example, for a paired-end, short-read metabolic labeling experiment simulated on November 3, 2025, these are the output files:


```
shortread_metaboliclabeling_PE_20251103_R1.fastq.gz
shortread_metaboliclabeling_PE_20251103_R2.fastq.gz
```


The remainder of this README describes the full pipeline and every parameter in detail, including how to build custom conditions that go beyond the presets.
## Overview of SPARK
The SPARK pipeline was designed to be run either as a single continuous pipeline or in individual modes to select genes to simulate, assign transcription elongation rates, simulate RNA biogenesis, and generate reads from those simulated molecules. Here, we describe the usage of the fullpipeline mode and below, we provide sections that discuss each of the other modes individually and alternative parameter usage.
```
usage: spark.py


general arguments:
  -h, --help            show this help message and exit
  --mode {fullpipeline,gene_selection,rates_per_region,mRNAgeneration,readgeneration}
                        pipeline step to run (default: fullpipeline)
  -o O                  output directory (default: ./)
  --seed SEED           random seed for reproducibility (default: None)
  --experiment_type {metaboliclabeling,proseq,ttseq,netseq,chromatin_associated}
                        type of experiment to simulate (default: metaboliclabeling)
  --threads THREADS     number of threads for parallel processing (default: 1)
  --remove_temp         delete intermediate files as soon as no later step needs them, keeping only the final
                        libraries, ground truth files and the log. Files still needed by a later --mode step are
                        kept until that step runs (default: False)


presets:
  --preset {default,DRB4sUseq,4sUseq,TTseq,BrUseq,NETseq,proseq,chromatinassociated,TT-TimelapseSeq,eSLAM-seq,polyA}
                        load a predefined set of parameters (default: default)
  --show_presets        display all available presets and their default values, then exit


gene selection:
  --gtf GTF             genome annotation file
  --genome_fasta GENOME_FASTA
                        genome fasta file
  --protein_coding_only
                        only use protein-coding genes (default: False)
  --n_gene_clusters N_GENE_CLUSTERS
                        number of gene clusters to simulate (default: 1)
  --n N                 number of genes to simulate (set to 0 to simulate all genes in the GTF) (default: 0)


rates per region:
  --region_size_range REGION_SIZE_RANGE
                        range of elongation rate region sizes to split the gene into (default: 500,5000)
  --num_regions NUM_REGIONS
                        number of elongation rate regions to split the gene into
  --elong_rate_range ELONG_RATE_RANGE
                        range of elongation rates in nt/min (default: 500,5000)
  --num_pauses NUM_PAUSES
                        number of pausing events across the gene (default: 0)
  --pause_time PAUSE_TIME
                        length of the pausing event in minutes (default: 0.1,0.5)
  --flat_rates          all nucleotides have the same elongation rate (default: False)
  --promoter_pause_position PROMOTER_PAUSE_POSITION
                        range of positions for TSS pausing (default: 30,50)
  --promoter_pause_duration PROMOTER_PAUSE_DURATION
                        length of promoter pausing in minutes (default: 0,0)
  --initi_rate INITI_RATE
                        range of seconds between initiations, sampled log-uniformly per gene. Sets each gene's
                        expression (default: 5,200)


pre-mRNA generation:
  --experiment_time EXPERIMENT_TIME
                        experiment duration in minutes (default: 15)
  --seq_err SEQ_ERR     sequencing error rate (default: 0.0001,0.0002)
  --nt_inc_prob NT_INC_PROB
                        range of nucleotide incorporation probability (default: 0.005,0.023)
  --subs_prob SUBS_PROB
                        range of nucleotide substitution probability if simulating nucleotide recoding (default: 0.95,1)
  --labeling_base LABELING_BASE
                        base analog(s) incorporated into the nascent RNA. A single base or a comma-separated
                        list, e.g. "T" or "A,T" (default: T)
  --sub_base SUB_BASE   identity of the base(s) after conversion, matched positionally to --labeling_base, e.g.
                        "C" or "G,C" (A->G and T->C). A single value is applied to every labeling base. Omit
                        for incorporation without recoding (default: None)
  --bkg_molecules BKG_MOLECULES
                        fraction of each gene's reads that come from background (non-nascent) molecules (default: 0)
  --drb                 simulate a DRB treatment experiment for transcription synchronization (default: False)
  --mRNAcount MRNACOUNT
                        size of the molecule pool simulated per gene, including background if specified. Sets
                        simulation size, not expression (default: 5000)


pre-mRNA processing:
  --cleavage_half_life CLEAVAGE_HALF_LIFE
                        range of cleavage half-lives in min, randomly chosen per gene (default: 0.1,10.0)
  --chromatin_release_half_life CHROMATIN_RELEASE_HALF_LIFE
                        min,max chromatin release half-life in minutes (default: 2.0,5.0)
  --pas_elong_rate PAS_ELONG_RATE
                        range of elongation rates in kb/min after passing the polyA site (default: 0.5,2.0)
  --readthrough_deg_rate READTHROUGH_DEG_RATE
                        degradation rate of read-through transcripts after cleavage in kb/min (default: 1.56)
  --nocleavage          do not simulate cleavage (default: False)
  --immediate_cleavage  simulate immediate cleavage at the polyA site without a half-life delay (default: False)
  --intron_half_life INTRON_HALF_LIFE
                        range of intron half-life in minutes used to simulate splicing (default: 1,10)
  --nosplicing          do not simulate splicing (default: False)
  --polyA_length POLYA_LENGTH
                        range of polyA tail lengths to add to cleaved mRNAs (base range, +/- 10% jitter) (default: 50,150)
  --mRNA_coordinates    output a TSV with the coordinates of the surviving features (exons + retained introns) (default: False)


sequencing strategy:
  --seq_tech {longread,shortread}
                        sequencing technology (default: shortread)
  --seq_type {RNA,cDNA,SE,PE}
                        sequencing type (default: SE)
  --seq_depth SEQ_DEPTH
                        total library sequencing depth (default: 20000000)
  --insert_size INSERT_SIZE
                        length of insert (default: 200,300)
  --read_length READ_LENGTH
                        length of each read (default: 100)
  --s {rf,fr,unstranded}
                        library strandedness (default: rf)
  --fragments           export a ground truth for fragmentation & size selection (default: False)
  --sizeselectiontype {none,hardcut,probabilistic}
                        type of size selection. 'hardcut' strictly filters fragments outside the --insert_size
                        limits, while 'probabilistic' applies a double-sided sigmoid retention curve around
                        those limits (default: probabilistic)
  --no_fragmentation    if specified, do not fragment (only valid for metaboliclabeling and ttseq) (default: False)
```
Snippet to simulate 150nt rf-firststrand paired-end reads from 1000 genes from a TT-seq library:
```
spark.py --mode fullpipeline -o ./out_dir/ --experiment_type ttseq --n 1000 --gtf GRCh38.95.gtf --genome_fasta GRCh38.fa --s rf --seq_type PE --experiment_time 5 --read_length 150 --seed 1
```
SPARK can be run in full or as specific modules to change specific parameters of the simulation while keeping others constant. Users can choose to input a seed to replicate the randomly assigned values (e.g., elongation rates, intron splicing rates, RNA fragmentation, etc) across runs. With the same `--seed`, every per-gene parameter (selected genes, elongation rates and pauses, initiation rate, cleavage, splicing and chromatin release half-lives) is identical across runs and experiment types. This is the way to simulate several assays (e.g., TT-seq and NET-seq) on the same genes and kinetics: run each assay with the same `--seed` into its own output directory. When running modes separately, pass the same simulation parameters (e.g., `--experiment_time`, `--mRNAcount`, `--bkg_molecules`, labeling options) in every call: if read generation needs to simulate a gene again (see [Sequencing depth and gene expression](#sequencing-depth-and-gene-expression)), it reuses that gene's stored per-gene parameters but takes the other settings from the current call.
### Reusing Simulated Data Across Runs


SPARK can share simulated molecules between runs. To do this, simply use the same output directory (`-o`) for each subsequent run. Reuse only works between experiment types that simulate the same molecules: `metaboliclabeling` and `ttseq` share them, so a metabolic labeling run can be re-sequenced as TT-seq, or with a different sequencing technology, read type, depth or insert size. NET-seq and PRO-seq simulate a different molecule population (elongating polymerases), and chromatin-associated runs filter molecules during generation, so these need their own simulation (same `--seed`, separate output directory).


```
# Run 1: Generates simulated mRNAs for a metabolic labeling experiment
spark.py --mode fullpipeline -o ./out_dir/ \
 --experiment_type metaboliclabeling --n 1000 \
 --gtf GRCh38.95.gtf --genome_fasta GRCh38.fa \
 --s rf --seq_type PE --experiment_time 15 --seed 1
```


Next, to generate a TT-seq library from the exact same simulated molecules, point to the same output directory and run only read generation, with the same simulation parameters:


```
# Run 2: Reuses the molecules from Run 1 to generate a TT-seq library
spark.py --mode readgeneration -o ./out_dir/ \
 --experiment_type ttseq --experiment_time 15 --seed 1 \
 --s rf --seq_type PE
```
### Removing intermediate files
By default, SPARK keeps every intermediate file it creates in the output directory (per-gene sequences, pre-mRNA and mRNA tables, per-gene reads, etc.). This is what allows reusing simulated data across runs, but it can take a lot of disk space for large simulations. Adding `--remove_temp` deletes each intermediate file as soon as no later step needs it, which also keeps disk usage low while the pipeline is running. At the end of a full run, the output directory only contains:
- `final_libraries/`
- `ground_truth_per_gene.tsv`
- `rate_per_gene/` (elongation rates & pausing ground truth)
- `intron_half_lives/` (when splicing is simulated)
- `ground_truth_pre_size_selection/` and `ground_truth_after_size_selection/` (when `--fragments` is used)
- `spark_pipeline.log`

When running individual modes, files needed by a later mode are kept until that mode runs (e.g., `--mode mRNAgeneration --remove_temp` keeps the simulated mRNAs for `--mode readgeneration`). Because intermediate files are deleted, the output of a run with `--remove_temp` cannot be reused by later runs as described above.
### Step 1: Gene Selection
**Clustering by gene features (`--mode gene_selection`)**
SPARK is designed to select genes from which to simulate nascent RNA reads that represent a diversity of genomic features. Users specify the number of genes from which they would like to simulate reads and provide: (1) a gtf with gene annotations and (2) a reference fasta file. For each gene, the longest isoform (upstream-most transcript start site and downstream-most transcript end site) is selected, and its features are used for classification. Genes in the gtf are first categorized by gene length, total transcript length(s), number of exons, mean exon lengths, number of introns, mean intron lengths, first intron length, lengths of 5′ and 3′ untranslated regions (UTRs), and exonic and intronic nucleotide composition. Hierarchical clustering is used to identify clusters of transcripts with similar characteristics.
The number of genes to simulate can be established using `--n` (set `--n 0` to simulate every gene in the GTF), and the number of clusters to divide genes into can be changed with `--n_gene_clusters`. An equal number of genes is selected from each cluster for downstream simulations, ensuring a balanced and unbiased set of genes reflecting genome diversity. SPARK can also simulate genes for any species with annotated genomes and selected genes from a pre-filtered gtf file.


To simulate multiple isoforms per gene, users can provide a gtf where each isoform to simulate has a different gene id.
[![](https://github.com/ezecalvo/spark/raw/master/figures/gene_clustering.png)](/ezecalvo/spark/blob/master/figures/gene_clustering.png)
### Step 2: Elongation Rates
**Assigning an RNA Pol II elongation rate per nucleotide (`--mode rates_per_region`)**
RNAPII elongation rates are likely not static across a gene. Additionally, any changes in elongation rates likely occur continuously, rather than in a step-wise fashion. To simulate the variability in elongation rates across the gene sequence, this mode will assign an elongation rate (in the range specified by `--elong_rate_range`) to each nucleotide of the gene. Elongation rates within this range are assigned to regions whose length is chosen from the range specified by `--region_size_range`. Within each region, the elongation rate (in nt/min) changes linearly from the rate at the first nucleotide to the rate at the last nucleotide in that region. For example, in a 10nt region that starts at 1 nt/min and ends at 2 nt/min, nucleotide #1 = 1 nt/min, #2 = 1.11, #3 = 1.22, ... until #10 at 2 nt/min.
[![](https://github.com/ezecalvo/spark/raw/master/figures/rates_per_position.png)](/ezecalvo/spark/blob/master/figures/rates_per_position.png)
To simulate step-wise changes in elongation rates, users can specify `--flat_rates` to ensure all the nucleotides in a region have the same elongation rate (each region still gets its own rate). To simulate a single elongation region across the entire gene, users can specify the number of regions with the same elongation rate per nucleotide using `--num_regions`. Specifying both `--flat_rates` and `--num_regions 1` will create a single constant elongation rate across the entire gene.
[![](https://github.com/ezecalvo/spark/raw/master/figures/rates_flat_genelevel.png)](/ezecalvo/spark/blob/master/figures/rates_flat_genelevel.png)
During the simulated experiment, each gene initiates new transcripts as a Poisson process. The mean interval between initiations is drawn per gene, log-uniformly, from the range set by `--initi_rate` (in seconds, default 5–200 s, a 40-fold range). The initiation rate is what sets each gene's expression: the number of molecules a gene produces, and therefore its reads in every assay, scales with it.


Finally, users can choose to include RNAPII pausing events throughout the gene. Pause sites can be placed randomly across the gene. The number and duration of those events can be modulated with `--num_pauses` and `--pause_time`, respectively. TSS pausing can also be added using `--promoter_pause_position` and `--promoter_pause_duration`.
[![](https://github.com/ezecalvo/spark/raw/master/figures/rates_with_pause.png)](/ezecalvo/spark/blob/master/figures/rates_with_pause.png)
#### Output files
Two output files per gene will be created in a new directory named *rate_per_gene* to record the ground truth of elongation rates. The *RatesandTraversalTimes* files contain the elongation rate per nucleotide, while the *VariableElongationRateRegions* files contain the elongation rate variability information per region (see the [ground truth](#ground-truth) section for details of file format).
### Step 3: Nascent RNA generation
**pre-mRNA generation (`--mode mRNAgeneration`)**
SPARK first simulates full-length pre-mRNA isoforms based on the elongation rates generated by `--mode rates_per_region`. Using these rates, transcription is simulated to create nascent pre-mRNA transcripts. The length of the experiment is specified by `--experiment_time`, which is used to calculate a total traveled distance by each RNAPII considering the cumulative elongation rates across the gene. These molecules can be labeled with a nucleotide analog or not, depending on the specified `--experiment_type`.

Each gene is simulated as a population of independent gene copies (analogous to the same locus across a population of cells). At the start of the experiment, each copy carries the number of RNAPII expected at steady state, *N*<sub>RNAP</sub> = *i* × *T*, where *i* is the gene's initiation rate and *T* is the time needed to cross the gene (the sum of the traversal times of its nucleotides, including pauses). RNAPII are placed along the gene in proportion to the time they spend in each region, and new transcripts initiate during the experiment at rate *i*. `--mRNAcount` sets the size of the molecule pool simulated per gene (`--mRNAcount` × `--experiment_time` / 10, ×10 for NET-seq and PRO-seq): SPARK simulates as many gene copies as needed to fill the pool and records how many copies the pool represents. `--mRNAcount` therefore controls simulation size and runtime, not expression (see [Sequencing depth and gene expression](#sequencing-depth-and-gene-expression)).

For nucleoside-analog dependent experiment types (`--experiment_type metaboliclabeling`, `--experiment_type proseq`, or `--experiment_type ttseq`), SPARK will simulate temporally-resolved nascent RNA sequencing using metabolic labeling. There are several parameters that determine the labeling probability, detection, and efficiency. During this "time", the probability that a labeled nucleotide specified with `--labeling_base` is *in silico* incorporated into nascent RNA molecules is specified by `--nt_inc_prob`. Once incorporated, the probability that these labeled nucleotides undergo nucleoside re-encoding is specified by `--subs_prob`. The nature of this re-encoding is specified by `--sub_base`. For example, iodoacetamide or TFEA-induced re-encoding of 4sU results in a T-to-C conversion, so `--labeling_base T --sub_base C` would simulate eSLAM-seq or enriched TimeLapse-seq experiments. Several bases can be labeled at once by giving `--labeling_base` a comma-separated list (e.g., `A,C,G,T`). `--sub_base` is then matched positionally (e.g., `--labeling_base A,T --sub_base G,C` simulates A-to-G and T-to-C conversions), and a single `--sub_base` value is applied to every labeling base.
Two transcriptional conditions can be simulated using `--drb`:
(1) RNA enrichment after transcription synchronization with DRB (*i.e.* 4sUDRB-seq), such that all RNAPIIs initiate at the transcription start site at the start of the labeling window (including `--drb`).
(2) A native context assuming steady-state transcription, in which RNAPII are already distributed across the gene body when labeling begins (default, not including `--drb`).
### Step 4: RNA processing
After nascent pre-mRNAs are generated, SPARK simulates the co-transcriptional RNA processing events of splicing, 3′-end cleavage, and polyadenylation. These steps are applied to the population of nascent molecules before enrichment.
#### Splicing
SPARK simulates co-transcriptional splicing as a first-order kinetic process per intron (`spark_spliceosome.py`). Each intron is assigned a half-life drawn from the range set by `--intron_half_life`. For every molecule, the time available for splicing an intron is the interval between when its 3′ splice site finished being transcribed and the molecule's current position in time; an intron is removed with probability `1 - exp(-t / mean_life)`, where `mean_life = half_life / ln(2)`. As a result, longer-lived introns and recently transcribed introns are more likely to be retained, recapitulating the co-transcriptional nature of splicing. Splicing can be turned off entirely with `--nosplicing`, in which case transcripts are treated as continuous, unspliced isoforms.
[![](https://github.com/ezecalvo/spark/raw/master/figures/splicing.png)]
#### Cleavage and polyadenylation
Once RNAPII transcribes past the polyA site (PAS), the transcript can be cleaved and polyadenylated. Cleavage is modeled as a first-order process with a per-gene half-life drawn from `--cleavage_half_life`. Upon cleavage, the molecule is split into two products: a *cleaved body* (the mature transcript upstream of the PAS, which receives a polyA tail) and a *read-through tail* (the downstream RNA still attached to the elongating polymerase). After the PAS, elongation of the read-through species proceeds at the rate set by `--pas_elong_rate`, and the uncleaved/read-through RNA is degraded at the rate set by `--readthrough_deg_rate`. PolyA tail lengths are drawn from `--polyA_length` (with ±10% jitter). `--immediate_cleavage` cleaves at the PAS with no half-life delay (every molecule that reaches the PAS is cleaved instantly), and `--nocleavage` disables cleavage entirely so molecules remain uncleaved. These cleavage states (`_uncleaved`, `_cleaved_body`, `_readthrough_tail`) are tracked per molecule and are what the chromatin-associated enrichment step (below) acts on.
[![](https://github.com/ezecalvo/spark/raw/master/figures/cleavage.png)]
### Step 5: Enrichment
Following the synthesis and processing of RNA molecules, SPARK can generate reads from one of several nascent RNA experimental approaches: metabolic labeling with full-molecule nascent RNA enrichment, metabolic labeling with fragmented nascent RNA enrichment, RNAPII run-on, nascent RNA associated with elongating RNAPIIs, or the chromatin-associated RNA fraction.
#### Metabolic labeling RNA enrichment
The metabolic labeling enrichment mode (`--experiment_type metaboliclabeling`) simulates a 4sU-seq style nascent RNA pulldown experiment by conditioning on full-length molecules that have at least one nucleotide analog incorporated.
#### Transient transcriptome (TT)-seq
The TT-seq experimental mode (`--experiment_type ttseq`) simulates a TT-seq style experiment by performing *in silico* fragmentation (see below) directly after full-length RNA generation and labeling but before nascent RNA enrichment. After this fragmentation, SPARK enriches for nascent RNA fragments (rather than nascent RNA molecules). Fragments without a nucleotide analog in their sequenced ends (the first and last `--read_length` nucleotides of the fragment) are discarded from the population of fragments before proceeding to library preparation.
#### Precision Nuclear Run-on (PRO)-seq
The PRO-seq experimental mode (`--experiment_type proseq`) simulates a run-on assay by simulating the incorporation of a biotinylated nucleotide during the *in silico* transcription process. The most upstream biotin-NTP incorporation event then defines the 3′ terminus of the nascent transcript by truncating the molecule at that nucleotide, and molecules without a biotin-NTP are discarded. From each truncated molecule, SPARK keeps a single 3′-end fragment ending at the biotin-NTP, with a length drawn uniformly from the `--insert_size` range, for subsequent library preparation simulation. The `proseq` preset allows biotin-NTP incorporation at all four bases (`--labeling_base A,C,G,T`), so the run-on can stop at any nucleotide.
#### Native Elongating Transcript (NET)-seq
The NET-seq experimental mode (`--experiment_type netseq`) simulates RNAPII pulldown-based nascent RNA sequencing assays (e.g., NET-seq and mNET-seq) by identifying all molecules still engaged with an elongating RNAPII at the end of the simulated experiment time. This includes RNAPII downstream of the PAS whose read-through RNA has not yet been degraded up to the polymerase. The simulation mimics RNAPII pulldown followed by MNase digestion. For each engaged molecule, SPARK extracts the 3′-terminal fragment that represents the nascent RNA physically protected by the RNAPII complex by selecting the last N nt of the molecule, where N is drawn from a uniform distribution between 35–100 nucleotides. Only these 3′ fragments are retained for subsequent library preparation steps. NET-seq does not use nucleotide analogs, so no analog incorporation is simulated in this mode.


#### Chromatin-associated RNA
The chromatin-associated experimental mode (`--experiment_type chromatin_associated`) simulates a caRNA-seq style experiment by enriching for the RNA still tethered to chromatin (`chromatin_isolation.py`). Molecules that have not yet been cleaved (`_uncleaved`) and read-through tails (`_readthrough_tail`) are assumed to remain chromatin-associated and are retained in full. Cleaved bodies (`_cleaved_body`), by contrast, are released from chromatin over time: each is retained with probability `exp(-t / mean_life)`, where `t` is the time since cleavage and the chromatin release half-life is drawn per gene from `--chromatin_release_half_life`. The resulting chromatin-associated population is therefore enriched for actively transcribing and recently cleaved RNA, while mature, released transcripts are progressively depleted.
[![](https://github.com/ezecalvo/spark/raw/master/figures/experiment_types.png)](/ezecalvo/spark/blob/master/figures/experiment_types.png)
#### Sequencing errors and background molecules
SPARK also allows for the introduction of two background signal sources:

**Errors** (substitutions) are introduced by technical steps in the experiment (e.g., base calling from sequencing). These errors can be incorporated by specifying `--seq_err`, which sets the probability that there will be a random substitution at any base in the molecule.

**Background, non-nascent transcripts.** These are transcripts that are not transcribed during the labeling window, but may be present due to inefficient or non-specific enrichment. `--bkg_molecules` sets the fraction of each gene's reads that come from background molecules: for every nascent read, a gene receives `bkg / (1 − bkg)` background reads, so background is exactly that fraction of every gene's reads. Background molecules are full-length mature (spliced) transcripts, or unspliced transcripts with `--nosplicing`, and carry no nucleotide analogs. They go through the same library preparation as the nascent molecules of the chosen `--experiment_type`:
- For metabolic labeling, fragments of the full-length background molecules are kept even though they carry no analog.
- For TT-seq, unlabeled fragments of background molecules are kept after the fragmentation step.
- For NET-seq, each background molecule contributes a random fragment of 35–100nt.
- For PRO-seq, each background molecule contributes a 3′-end fragment, as nascent molecules do.
- For chromatin-associated, background molecules bypass the chromatin release step.
- For long reads, background molecules are sequenced whole.

Setting `--bkg_molecules 1` (as in the `polyA` preset) produces a library made entirely of background (mature) molecules, approximating a steady-state polyA RNA-seq library. Reads are then split across genes in proportion to each gene's initiation rate and the number of fragments each mature molecule yields (proportional to initiation rate × mature transcript length for short reads, and to initiation rate for long reads).
### Step 6: Read Generation
**Read generation (`--mode readgeneration`)**
SPARK can generate either short-read or long-read sequencing reads from the RNA molecules produced by the enrichment step. The choice of sequencing technology is specified by `--seq_tech`, which has options for *shortread* (e.g. Illumina) or *longread* (e.g., Oxford Nanopore Technologies).

For short-read simulations, fragmentation and size selection are modeled after the simulations described in [Pai *et al.*, eLife 2017](https://elifesciences.org/articles/32537). Fragment insert sizes are specified by `--insert_size`. Read length is specified by `--read_length`. The simulation of single-end (*SE*) or paired-end (*PE*) is specified by `--seq_type`. Finally, strandedness is specified by `--s`, choosing from reverse-forward (*rf*), forward-reverse (*fr*), or *unstranded*. The default is *rf*, which is typical for dUTP-based strand-specific RNA-seq protocols.
For long-read data, SPARK defaults to Oxford Nanopore strand specifications: for directRNA, transcripts are sequenced in the 3′ to 5′ direction, while cDNA reads may originate from either end. This behavior is specified by `--seq_type`, which can be set to "RNA" or "cDNA" for long reads.
#### Sequencing depth and gene expression
SPARK does not assign expression values directly. Each gene's read count follows from its simulated kinetics:
1. After enrichment and library preparation, SPARK counts the sequenceable units in each gene's molecule pool: size-selected fragments (metabolic labeling, TT-seq, chromatin-associated), 35–100nt footprints (NET-seq), 3′ run-on fragments (PRO-seq), or whole molecules (long reads).
2. Dividing by the number of gene copies the pool represents gives the signal that one cell produces for that gene in the chosen assay.
3. `--seq_depth` reads are split across genes in proportion to that per-cell signal (after setting aside the background fraction), so the final library contains exactly `--seq_depth` reads.
4. If a gene is allocated more reads than its pool can supply, SPARK simulates that gene again with proportionally more gene copies before writing reads, so reads are never resampled. The gene keeps the per-gene parameters drawn the first time (initiation rate, post-PAS elongation rate, cleavage and chromatin release half-lives), which SPARK stores in `rate_per_gene/<gene>_parameters.tsv`. At realistic depths this affects only the most highly expressed genes.

As a result, expression is controlled by `--initi_rate` and follows each assay's kinetics. For example, TT-seq reads scale with the initiation rate × gene length (RNA synthesized during labeling), while NET-seq reads scale with the number of engaged RNAPII (initiation rate × time to cross the gene), so a slowly elongating gene has more NET-seq reads relative to TT-seq than a fast one. The RPKM and TPM computed from the reads each gene received are reported in the [ground truth](#consolidated-per-gene-ground-truth).
### Ground truth
SPARK can output multiple types of ground truth that can be used to validate a diversity of benchmarking metrics.
#### Consolidated per-gene ground truth
At the end of a run, SPARK consolidates the per-gene ground truth metrics produced by each module into a single `out_dir/ground_truth_per_gene.tsv` file. Depending on the simulated conditions, this table can include:

| Column | Description |
| --- | --- |
| `initiation_rate` | mean interval between initiation events, in seconds |
| `tss_pause_position`, `tss_pause_duration_min` | TSS pause position and duration |
| `min_elongation_rate`, `max_elongation_rate`, `mean_elongation_rate` | elongation rates across the gene (nt/min); the mean is the arithmetic mean over nucleotides |
| `num_pauses`, `pause_duration_min`, `pause_positions` | gene body pauses |
| `cleavage_half_life` | cleavage half-life (min) |
| `PAS_elongation_rate` | elongation rate after the PAS (kb/min) |
| `rpkm`, `tpm` | expression computed from the reads each gene received, using its exonic length |
| `mean_seq_err` | mean sequencing error rate |
| `chromatin_release_halflife_min` | chromatin release half-life (chromatin-associated runs) |

For metrics that depend on the time polymerases spend on a gene (e.g., elongation velocity from TT-seq/NET-seq), the matching average rate is gene length divided by the sum of `time_for_this_nt` in the *RatesandTraversalTimes* file, rather than `mean_elongation_rate`. The two differ when rates vary along the gene.






#### Elongation rates & pausing
1. Regions over which elongation rates are varied, including the position of pausing events, which are stored in the `out_dir/rate_per_gene/*_VariableElongationRateRegions.gtf` files containing the elongation rate information per region:
| Column Name | Description |
| --- | --- |
| chromosome | gene chromosome |
| absolute_start | absolute start position of the region |
| absolute_end | absolute end position of the region |
| region_start_coord | relative start position of the region |
| region_end_coord | relative end position of the region |
| strand | strand |
| rate_initial | elongation rate in bp per minute for the first nucleotide in the region |
| rate_final | elongation rate in bp per minute for the last nucleotide in the region |
| rate_change_per_nt | change in rate for each individual nucleotide in the region |
| time_to_traverse | time in minutes to traverse the entire region |
| sequence | region sequence |
2. Per nucleotide elongation rates, including the position of pausing events, which is stored in the `out_dir/rate_per_gene/*_RatesandTraversalTimes.gtf` files:
| Column Name | Description |
| --- | --- |
| chromosome | gene chromosome |
| absolute_position | absolute position of each nucleotide |
| strand | strand |
| region_number | gene region the nucleotide belongs to |
| nucleotide_coord | relative position of the nucleotide within the gene |
| time_for_this_nt | time in minutes required to traverse the nucleotide |
| rate_for_this_nt | elongation rate in bp per minute for the nucleotide |
| rate_change_per_nt | elongation rate change between nucleotide n and n+1 |
#### Splicing ground truth
When splicing is simulated, SPARK records the per-intron half-lives used for each gene (`out_dir/intron_half_lives/`) together with a validation file comparing the input half-life to the half-life re-estimated from the simulated molecules. If `--mRNA_coordinates` is set, a per-molecule TSV of surviving features (exons + retained introns) is also written, capturing exactly which introns were removed from each molecule.

#### Nucleoside analog incorporation & substitution
For metabolic labeling and chromatin-associated short reads, the read names contain the number and position of incorporated and converted nucleotides (if any) in the read. For example:
```
@XWXMPWMRCUZ_uncleavedtzkbq_12:108628977-108628878:-_ninc2:27,70_nsubs0:
Molecule ID: XWXMPWMRCUZ
Cleavage state of the molecule: uncleaved
Read suffix: tzkbq
Number of modified nucleotides in the read: 2
Position of the modified nucleotides in the read: 27,70
Number of substitutions: 0
Position of the substitutions in the read: NA
```


```
@2LOX927M95I_cleaved_bodyfdjie_12:108626879-108626780:-_ninc3:17,21,90_nsubs1:21
Molecule ID: 2LOX927M95I
Cleavage state of the molecule: cleaved_body
Read suffix: fdjie
Number of modified nucleotides in the read: 3
Position of the modified nucleotides in the read: 17,21,90
Number of substitutions: 1
Position of the substitutions in the read: 21
```


Each read name begins with an 11-character capitalized molecule ID, followed by the molecule's cleavage state (`_uncleaved`, `_cleaved_body` or `_readthrough_tail`); reads coming from the same molecule share this prefix. A 5-letter lowercase suffix then differentiates individual reads. For read pairs, both reads will have identical IDs.

TT-seq, NET-seq and PRO-seq read names have the same structure without the incorporation and substitution fields (e.g., `@XWXMPWMRCUZ_uncleavedtzkbq_12:108628977-108628878:-`). Long-read names contain the molecule ID followed by the gene ID.

Background molecules, coming from unspecific enrichment steps, will contain 'BG' in their names but no nucleoside analog incorporation or substitution due to them being background:
```
@69UOVR0SU3O_BG_zvhzi_12:108625373-108625274:-_ninc0:_nsubs0:
```
#### Fragmentation and size selection
For experiments that involve a size selection step (e.g., `--experiment_type metaboliclabeling` or `--experiment_type ttseq`), you can save the ground truth per-base coverage data from before and after this step. To enable this, add the `--fragments` argument to your command.
When this option is used, SPARK will create two directories in your output folder, containing files with the per-base coverage data:
`out_dir/ground_truth_pre_size_selection/` This directory contains the coverage data after fragmentation but before the size selection is applied.
`out_dir/ground_truth_after_size_selection/` This directory contains the final coverage data after the size selection has been applied.
For metabolic labeling, both files describe each gene's full simulated fragment pool. For TT-seq, they are sampled to the number of reads the gene received.
#### For read mapping
To use SPARK to benchmark read alignment under different experiment simulation conditions, the read names contain the genomic coordinates from which each read was simulated:
```
@CN1369NOND9_uncleavedqjfsv_12:108625124-108625025:-_ninc5:52,50,37,21,94_nsubs0:
Molecule ID: CN1369NOND9
Chromosome: 12
Start: 108625124
End: 108625025
Strand: -
Number of modified nucleotides in the read: 5
Position of the modified nucleotides in the read: 52,50,37,21,94
Number of substitutions: 0
Position of the substitutions in the read: NA
```
For TT-seq, NET-seq and PRO-seq reads that fall outside the annotated gene (e.g., in the read-through region downstream of the PAS), the coordinates are reported as `NA`. For metabolic labeling and chromatin-associated reads, they are extrapolated from the gene's coordinates.




