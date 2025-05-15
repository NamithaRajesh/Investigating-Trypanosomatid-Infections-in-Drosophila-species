Sure! Here's your **Commands_Run.md** file formatted to match the structure you prefer:

```md
# RNA-seq Analysis Project: Command Line Operations  

## Overview  
This document contains the command-line operations used in the RNA-seq analysis project, detailing steps from quality control to alignment and classification of unmapped reads.  

## Repository Structure  
- `project_details.md` – Background information about the project.  
- `methods.md` – Steps taken so far, including plots explaining the actions taken.  
- `commands.md` – This file containing all command-line operations used.  
- `references.md` – A list of tools, datasets, and publications referenced in the project.  

## Quality Control & Preprocessing  
RNA-seq data requires careful preprocessing to ensure high-quality results. The **SRA Toolkit** provides utilities for working with sequencing data.  

### **Key Tools in the SRA Toolkit**
- `fastq-dump` – Converts SRA files into FASTQ format.  
- `prefetch` – Downloads SRA files from NCBI.  
- `fasterq-dump` – Optimized version of `fastq-dump`.  
- `sam-dump` – Converts SRA files into SAM format.  

#### **Downloading & Setting Up the Toolkit**
```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:/path/to/sratoolkit/bin
cd sratoolkit.3.1.1-ubuntu64/
source ~/.bashrc
```

#### **Processing RNA-seq Sample (SRR5647735)**
```sh
fastq-dump --split-3 -X 10000 SRR5647735
fastp --cut_tail -i SRR5647735_1.fastq -I SRR5647735_2.fastq -o SRR5647735_1.trim.fq -O SRR5647735_2.trim.fq
```

---

## Trimming & Adapter Removal  
FastQC analysis revealed issues requiring trimming. **Trimmomatic** was used to remove low-quality sequences and adapters.  

#### **Creating Adapter File**
```sh
echo ">illumina" > adapter.fa
echo "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" >> adapter.fa
```

#### **Testing Different Trimming Settings**
```sh
java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 4 -phred33 \
  SRR5647735_1.fastq SRR5647735_2.fastq \
  SRR5647735_1.trim.fq SRR5647735_1.unpaired.fq \
  SRR5647735_2.trim.fq SRR5647735_2.unpaired.fq \
  ILLUMINACLIP:adapter.fa:2:30:10 \
  SLIDINGWINDOW:4:20 \
  TRAILING:20 \
  MINLEN:36
```
✅ **Best Setting:** Successfully removed overrepresented sequences and adapter content, while addressing sequence irregularities.  

---

## Alignment to Reference Genome  
Short reads were aligned to the **Drosophila melanogaster** genome using **BWA**, a tool optimized for read mapping.  

#### **Setting Up Environment & Indexing Reference Genome**
```sh
conda update -n base -c conda-forge conda
conda create -n bwa_env bwa -c bioconda
conda activate bwa_env

bwa index dmel-all-chromosome-r6.61.fasta
bwa mem dmel-all-chromosome-r6.61.fasta SRR5647735_1.trim.fq SRR5647735_2.trim.fq > output1.sam
```

#### **Processing SAM/BAM Files with Samtools**
```sh
conda install -c bioconda samtools
conda create -n sam_env samtools -c bioconda
conda activate sam_env

samtools view -bS output_4.sam > output_4.bam
samtools view -b -f 4 output_4.bam > unmapped_4.bam
samtools fastq unmapped_4.bam > unmapped_reads_4.fastq
```

---

## Identification of Unmapped Reads  
To classify unmapped reads, **Kraken2** was used for metagenomic analysis.  

#### **Installing & Building Kraken2 Database**
```sh
sudo apt install kraken2
sudo kraken2-build --standard --db /home/hlnrajes/Project/kraken/kraken_database
sudo kraken2-build --build --db /home/hlnrajes/Project/kraken/kraken_database --threads 8
```

#### **Running Kraken2 for Classification**
```sh
kraken2 --db /home/hlnrajes/Project/kraken/kraken_database --threads 4 \
        --report output_4.report --output output_4.kraken unmapped_reads_4.fastq
```

#### **Extracting Species-Level Classifications**
```sh
awk '$4 == "S"' output_4.report > species_classification.txt
```

✅ **Final Observations:**  
- Trimming improved read quality significantly.  
- Alignment successfully mapped reads to the reference genome.  
- Kraken2 classified species from unmapped reads.  

## Future Work  
- Refining the classification of unmapped reads.  
- Investigating additional analysis techniques.  
- Incorporating visualization methods for improved reporting.  

```
