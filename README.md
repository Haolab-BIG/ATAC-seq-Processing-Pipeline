# ATAC-seq Pipeline

The ATAC-seq analysis pipeline processes raw FASTQ data through a series of steps including adapter trimming, quality control, genome mapping, peak calling,transcription factor occupancy prediction,annotate and compare peak data. Using Singularity for reproducibility and supporting batch analysis of multiple samples.

## Workflow Diagram

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/ATAC-seq_pipeline.png)

## Requirements

1. **Recommended System Configuration**:

 * 8-core CPU
 * 80 GB RAM

2. **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

 * **Step 1: Install System Dependencies**

```bash
# Update package lists and install dependencies
sudo apt-get update
sudo apt-get install -y \
build-essential \
libseccomp-dev \
libfuse3-dev \
pkg-config \
squashfs-tools \
cryptsetup \
curl wget git
```

 * **Step 2: Install Go Language**

```bash
# Download and install Go
wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
rm go1.21.3.linux-amd64.tar.gz
 
# Configure Go environment variables and apply them
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc
```

 * **Step 3: Download, Build, and Install Singularity**

```bash
# Note: The script navigates to /mnt/share/software. 
# You can change this to your preferred directory for source code.
cd /mnt/share/software
 
# Download the Singularity CE source code
wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz
 
# Extract the archive and clean up
tar -xvzf singularity-ce-4.0.1.tar.gz
rm singularity-ce-4.0.1.tar.gz
cd singularity-ce-4.0.1
 
# Configure the build
./mconfig
 
# Build Singularity (this can be time-consuming)
cd builddir
make
 
# Install Singularity to the system
sudo make install
```

 * **Step 4: Verify the Installation**

```bash
# Check the installed version
singularity --version

# Display help information
singularity -h
```
3. **snakemake**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. 

```bash
pip install snakemake
```

4. **Pipeline Files**:

 * `ATAC-seq.smk`
 * `atacseq.sif` (The Singularity container)

5. **Reference Data**: A directory containing all necessary reference files (e.g., bwa-mem2 indices and GTF annotation, etc.).

**Note on Sequencing Type:**
This pipeline supports paired-end (PE) data. The example below shows the format for paired-end data.

### 1. Prepare the Reference Genome

The pipeline requires several pre-built reference files. Below are the steps to generate them for the human hg38 genome using GENCODE annotations.

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p data
cd data
```

#### Common Reference Files

We require the following base files:

**Download Reference Files:**

```bash
# Download Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz

# Download GTF Annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz

# Unzip the files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v48.primary_assembly.annotation.gtf.gz
```

##### Build bwa-mem2 Indices:

```bash
# Build the main Genome Index
mkdir -p GRCh38.primary_assembly_genome.bwa-mem2_index
singularity exec atacseq.sif bwa-mem2 index -p GRCh38.primary_assembly_genome.bwa-mem2_index/GRCh38.primary_assembly_genome.bwa-mem2_index GRCh38.primary_assembly.genome.fa
```
### 2. Prepare the formated TF motifs data

The data run by this pipeline is backup of the jaspar, hocomoco, and CIS-BP databases from the meme database.The specific processing method is as follows

#### Create a dedicated directory for the MEME TF motifs data

```bash
mkdir -p motifs/individual
cd motifs
```

#### Download the MEME TF motifs data #### 

```bash
wget https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.27.tgz
tar -zxvf motif_databases.12.27.tgz
mv H12CORE_meme_format.meme  Homo_sapiens.meme  JASPAR2024_CORE_vertebrates_non-redundant_v2.meme individual
```

#### Format the MEME TF motifs data

```bash
singularity exec atacseq.sif TOBIAS FormatMotifs --input motifs/individual --format pfm --task join --output motifs/all_motifs.txt
```

### 3. Prepare the test fastq data

The data run by this pipeline is from SRR13171471,SRR13171472,SRR13171473 and SRR13171474 in the SRA database.The specific processing method is as follows

#### Create a dedicated directory for the test sra data:

```bash
mkdir -p data/samples
cd data/samples
```
#### Download the test sra data
```bash
prefetch SRR13171471
prefetch SRR13171472
prefetch SRR13171473
prefetch SRR13171474
```
#### Convert sra data to fastq data
```bash
fastq-dump --split-files SRR13171471\SRR13171471.sra
fastq-dump --split-files SRR13171472\SRR13171472.sra
fastq-dump --split-files SRR13171473\SRR13171473.sra
fastq-dump --split-files SRR13171474\SRR13171474.sra
```

#### Randomly sample fastq data:

```bash
singularity exec atacseq.sif seqtk sample -s100 SRR13171471_1.fastq 53086 > SRR13171471_R1.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171471_2.fastq 53086 > SRR13171471_R2.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171472_1.fastq 53086 > SRR13171472_R1.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171472_2.fastq 53086 > SRR13171472_R2.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171473_1.fastq 53086 > SRR13171473_R1.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171473_2.fastq 53086 > SRR13171473_R2.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171474_1.fastq 53086 > SRR13171474_R1.fastq
singularity exec atacseq.sif seqtk sample -s100 SRR13171474_2.fastq 53086 > SRR13171474_R2.fastq
pigz -p 8 SRR13171471_R1.fastq
pigz -p 8 SRR13171471_R2.fastq
pigz -p 8 SRR13171472_R1.fastq
pigz -p 8 SRR13171472_R2.fastq
pigz -p 8 SRR13171473_R1.fastq
pigz -p 8 SRR13171473_R2.fastq
pigz -p 8 SRR13171474_R1.fastq
pigz -p 8 SRR13171474_R2.fastq
```
### 4. Check the ATAC-seq snakemake workflow

The specific snakemake workflow is as follows

#### Dry-run:

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -np \
  -s ATAC-seq.smk \
  --use-singularity \
  --singularity-args "--bind /mnt/liuq/test/singularity/:/root/"
```
#### Draw a detailed DAG picture

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -s ATAC-seq.smk \
  --use-singularity \
  --singularity-args "--bind /mnt/liuq/test/singularity/:/root/" \
  --dag  | \
dot -Tsvg > dag.svg
```

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/dag.png)
### 5. Check the current working directory

The initial working structure for snakemake in /mnt/liuq/test/singularity/:
## Input Structure and Interpretation

Before the pipeline, the input directory contain several files and directories. Below is a detailed explanation of what each file you should supply and how it can be used.
```bash
├── atacseq.sif
├── ATAC-seq.smk
├── config.yaml
├── dag.svg
├── data
│   ├── gencode.v48.primary_assembly.annotation.gtf
│   ├── GRCh38.primary_assembly_genome.bwa-mem2_index
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.0123
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.amb
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.ann
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.bwt.2bit.64
│   │   └── GRCh38.primary_assembly_genome.bwa-mem2_index.pac
│   ├── GRCh38.primary_assembly.genome.fa
│   └── samples
│       ├── SRR13171471_R1.fastq.gz
│       ├── SRR13171471_R2.fastq.gz
│    ├── SRR13171472_R1.fastq.gz
│       ├── SRR13171472_R2.fastq.gz
│       ├── SRR13171473_R1.fastq.gz
│       ├── SRR13171473_R2.fastq.gz
│       ├── SRR13171474_R1.fastq.gz
│       └── SRR13171474_R2.fastq.gz
├── motifs
│   └── all_motifs.txt
├── samplesheet.txt
└── scripts
    ├── CHIPPEAKANNO_ATACseq.R
    └── CSAW_ATACseq.R
```
- **ATAC-seq.smk** — The main Snakemake workflow script.

- **config.yaml** — Configuration file containing paths, parameters, and sample information.
  ⚠️ Must be located in the same directory as `ATAC-seq.smk`.
  
- **atacseq.sif** — Singularity container image with all required software and dependencies pre-installed.

- **dag.svg**— Detailed pipeline DAG picture.

- **scripts** — Additional scripts required for program execution.

- **motifs** — TF motif file formatted by TOBIAS.

- **samplesheet.txt** — Experimental design table,describe the group of each sample.

  ```bash
  namecondition
  SRR13171471 control
  SRR13171472 control
  SRR13171473 test
  SRR13171474 test
  ```

- **data** — Reference genome index for bwa-mem2,samples which were sequenced and other required files which are listed as follows.

#### Per-Sample Files (`data/samples/`)

- **`genome.fa`**
  - **Content**: This is the genome fasta file.
  - **Application**: It is used for pause_index rule in snakemake workflow.
- **`genome.bwa-mem2_index`**
  - **Content**:  This is the whole genome index created by bwa-mem2,The building method is as above.
  - **Application**:  It's the primary reference for whole genome read alignment.
- **`annotation.gtf`**
  - **Content**: This is the genome annotation gtf file.
  - **Application**: It is used for annotate in snakemake workflow.
# Part III Running

After adding the corresponding parameters in config.yaml,you can execute the pipeline using a single command, the only important parameter is thread (**-j**).

**Example code**

* **Step 1: Edit `config.yaml`**

```shell
samples:
	KYSE30:
		SRR13171473: 
			R1 : "data/samples/SRR13171473_R1.fastq.gz"
			R2 : "data/samples/SRR13171473_R2.fastq.gz"
		SRR13171474: 
			R1 : "data/samples/SRR13171474_R1.fastq.gz"
			R2 : "data/samples/SRR13171474_R2.fastq.gz"
	HET1A:
		SRR13171471: 
			R1 : "data/samples/SRR13171471_R1.fastq.gz"
			R2 : "data/samples/SRR13171471_R2.fastq.gz"
		SRR13171472: 
			R1 : "data/samples/SRR13171472_R1.fastq.gz"
			R2 : "data/samples/SRR13171472_R2.fastq.gz"
fn_gtf : "data/gencode.v48.primary_assembly.annotation.gtf"
fa_in : 'data/GRCh38.primary_assembly.genome.fa'
BWA-MEM2IDX : "data/GRCh38.primary_assembly_genome.bwa-mem2_index/GRCh38.primary_assembly_genome.bwa-mem2_index"
home : "/mnt/liuq/test/singularity"
container : "atacseq.sif"
all_motifs : "motifs/all_motifs.txt"
species : "hsa"
```

* **Step 2: run snakemake**

Here /mnt/liuq/test/singularity/ represents the root directory.

```bash
snakemake -s ATAC-seq.smk \
		  -j 30 \
  --use-singularity \
  --singularity-args "--bind /mnt/liuq/test/singularity/:/root/"
```

**Command Parameters**

**edit `config.yaml`**

- `samples`: Path to the FASTQ file. For paired-end (required)
- `BWA-MEM2IDX`: Path to the directory where bwa-mem2 reference build with prefix (required)
- `container`: Path to the singularity environment file (required)
- `fa_in` : Reference Genome FASTA (required)
- `fn_gtf` : Reference Genome GTF Annotation (required)
- `home`: The root directory where running program (required)
- `all_motifs`: TF motif file formatted by TOBIAS (required)
- `species`: species which was sequenced (required)

**run snakemake**

- `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
- `--singularity-args`: Allows passing additional arguments to the Singularity runtime (e.g., `--bind`, `--nv`, or custom options).
- `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel when executing workflow rules.
- `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references. The format is `/project_directory:/project_directory`. Multiple directories can be mounted by separating them with commas, for example: `/path1:/path1,/path2:/path2` (required)

Then delete the intermediate files and folders

```bash
rm -rf qc mapped filtered_bam shifted_bam name_sorted sorted peaks group bias_correction deepTools_qc tmp data/multiqc_report/*/multiqc_data
rm Filtered.results.bed Filtered.results.DOWN.bed Filtered.results.MIXED.bed Filtered.results.UP.bed Full.results.bed TMM_normalizedCounts.pdf DiffBinding_scores.txt DiffBinding_modelfit.pdf DiffBinding_analysis.Rdata DiffBinding_allregions.bed CSAW.session_info.txt
rm TFBS/bindetect_figures.pdf TFBS/*.html TFBS/bindetect_results.xlsx TFBS/bindetect_distances.txt
cd TFBS
find . -type d -exec rm -rf {} \;
```
# Part IV Output

After the pipeline completes, the output directory will contain several files and directories. Below is a detailed explanation of what each file is and how it can be used.
```bash
├── atacseq.sif
├── ATAC-seq.smk
├── bigwig
│   ├── SRR13171471.bw
│   ├── SRR13171472.bw
│   ├── SRR13171473.bw
│   └── SRR13171474.bw
├── config.yaml
├── dag.svg
├── data
│   ├── gencode.v48.primary_assembly.annotation.gtf
│   ├── GRCh38.primary_assembly_genome.bwa-mem2_index
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.0123
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.amb
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.ann
│   │   ├── GRCh38.primary_assembly_genome.bwa-mem2_index.bwt.2bit.64
│   │   └── GRCh38.primary_assembly_genome.bwa-mem2_index.pac
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── multiqc_report
│   │   ├── SRR13171471
│   │   │   └── multiqc_report.html
│   │   ├── SRR13171472
│   │   │   └── multiqc_report.html
│   │   ├── SRR13171473
│   │   │   └── multiqc_report.html
│   │   └── SRR13171474
│   │   └── multiqc_report.html
│   ├── samples
│   │   ├── SRR13171471_R1.fastq.gz
│   │   ├── SRR13171471_R2.fastq.gz
│   │   ├── SRR13171472_R1.fastq.gz
│   │   ├── SRR13171472_R2.fastq.gz
│   │   ├── SRR13171473_R1.fastq.gz
│   │   ├── SRR13171473_R2.fastq.gz
│   │   ├── SRR13171474_R1.fastq.gz
│   │   └── SRR13171474_R2.fastq.gz
│   ├── SRR13171471
│   │   ├── SRR13171471_R1_val_1_fastqc.html
│   │   ├── SRR13171471_R1_val_1_fastqc.zip
│   │   ├── SRR13171471_R2_val_2_fastqc.html
│   │   └── SRR13171471_R2_val_2_fastqc.zip
│   ├── SRR13171472
│   │   ├── SRR13171472_R1_val_1_fastqc.html
│   │   ├── SRR13171472_R1_val_1_fastqc.zip
│   │   ├── SRR13171472_R2_val_2_fastqc.html
│   │   └── SRR13171472_R2_val_2_fastqc.zip
│   ├── SRR13171473
│   │   ├── SRR13171473_R1_val_1_fastqc.html
│   │   ├── SRR13171473_R1_val_1_fastqc.zip
│   │   ├── SRR13171473_R2_val_2_fastqc.html
│   │   └── SRR13171473_R2_val_2_fastqc.zip
│   └── SRR13171474
│       ├── SRR13171474_R1_val_1_fastqc.html
│       ├── SRR13171474_R1_val_1_fastqc.zip
│       ├── SRR13171474_R2_val_2_fastqc.html
│       └── SRR13171474_R2_val_2_fastqc.zip
├── diff
│   ├── annotated_peak.csv
│   ├── DiffBinding_significant.bed
│   ├── enrichmentPlot_GO.pdf
│   └── enrichmentPlot_PATH.pdf
├── motifs
│   └── all_motifs.txt
├── samplesheet.txt
├── scripts
│   ├── CHIPPEAKANNO_ATACseq.R
│   └── CSAW_ATACseq.R
└── TFBS
    └── bindetect_results.txt
```

- **`multiqc_report.html`**
  Open multiqc_report.html in a web browser to explore all sections interactively.
  
- **Application**: This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .
  
- **General Statistics**: A combined table summarizing important metrics for each sample:
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/general_statistic.png)
- **FastQC**: Quality-control metrics on raw and trimmed reads, including  
  'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
  'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
  'Sequence Length Distribution', 'Sequence Duplication Levels',  
  'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.

- **Sequence Quality Histograms**: The mean quality value across each base position in the read. 
  
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/fastqc_per_base_sequence_quality_plot.png)

  - **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/fastqc_per_sequence_quality_scores_plot.png)

- **`fastqc.html(zip)`**
  
  - **Content**: Open fastqc.html in a web browser to explore all sections interactively.Similar to the above multiqc results.
  
  - **Application**:This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .Similar to the above multiqc results.
  
- **`sample.bw`**
  
  - **Content**: BigWig files that represent the ATAC-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format.
  
  - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows chromatin opening levels visually across chromosomes. Highly chromatin opening will appear as peaks.
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/track.png)

- **`DiffBinding_significant.bed`**
  
  - **Content**: Bed files that represent differential chromatin accessibility regions.
  
  - **Application**: It is used to find differential epigenetic regulation regions between the case and control sample groups.![](D:\RLK_database\project\数据来源\毕业\博士后\杂事\EasyOmics\ATAC_seq\picture\DiffBinding_significant.png)

- **`annotated_peak.csv`**
  
  - **Content**: CSV format annotation file about regions of differential chromatin accessibility region.The "start" and "end" columns represent the start and end positions of the differential chromatin accessibility regions, while the "start position" and "end position" columns represent the start and end sites of the genes closest to the differential chromatin accessibility regions.
  
  - **Application**: It is used to indicate which genes contain differential chromatin accessibility regions and the target genes of cis-regulatory elements located in chromatin accessibility regions.
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/annotated_peak.png)

- **`enrichmentPlot_GO.pdf`**
  
  - **Content**: GO annotations for the genes of annotated_peak.csv.
  
  - **Application**:  This annotation implies the mechanism by which differential chromatin accessibility regions influence the phenotype about the case group relative to the control group.
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/enrichmentPlot_GO.png)

  

- **`enrichmentPlot_PATH.pdf`**
  - **Content**:  Biological pathways annotations for the genes of annotated_peak.csv.
  
  - **Application**:  Similar to GO annotations,This annotation implies the mechanism by which differential chromatin accessibility regions influence the phenotype about the case group relative to the control group.
  

![](https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/picture/enrichmentPlot_PATH.png)

**`bindetect_results.txt`**

- **Content**:  Transcription factors at chromatin-accessibility sites in both the case and control sample groups, including differentially binded transcription factors at chromatin-accessibility sites between the two groups.

- **Application**:  This result reveals the transcription factors that cause the phenotypic differences between the case sample group and the control sample group, and these transcription factors can serve as targets for drug therapy.

```bash
output_prefixname motif_id  cluster total_tfbsKYSE30_footprints_mean_score KYSE30_footprints_bound HET1A_footprints_mean_score  HET1A_footprints_bound  KYSE30_footprints_HET1A_footprints_change KYSE30_footprints_HET1A_footprints_pvalue KYSE30_footprints_HET1A_footprints_highlighted
_AHR.H12CORE.0.P.B  AHR.H12CORE.0.P.B C_6 64795.36961  1 82612.82311  2 -1.67202  7.01272E-66  False
_AHRR.H12CORE.0.P.C AHRR.H12CORE.0.P.CC_ZNF287  5 56837.34503  1 56710.89476  1 -0.8216 9.31477E-32  False
_ALX1.H12CORE.0.SM.BALX1.H12CORE.0.SM.B  C_2 51932.19363  0 35721.39629  0 0.23852 1.91053E-02  False
_ALX3.H12CORE.0.SM.BALX3.H12CORE.0.SM.B  C_5 60738.00209  1 104037.89701 3 -2.4982 2.74018E-73  False
_ALX3.H12CORE.1.S.B ALX3.H12CORE.1.S.BC_2 69028.88359  1 121701.66706 1 -3.63492  1.15582E-60  False
```

## Video Tutorials

Watch this video tutorial to see a complete walk through of running the pipeline:

https://github.com/Haolab-BIG/ATAC-seq-Processing-Pipeline/raw/main/snakemake.mp4
