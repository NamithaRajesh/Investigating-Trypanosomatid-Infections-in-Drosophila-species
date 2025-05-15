# Commands Run (Linux Environment)
Project Directory: `/home/hlnrajes/PROJECT/`

## Quality Control

### Key Tools in the SRA Toolkit
- `fastq-dump` – Converts SRA files into FASTQ format.  
- `prefetch` – Downloads SRA files from NCBI.  
- `fasterq-dump` – Optimized version of `fastq-dump`.  
- `sam-dump` – Converts SRA files into SAM format.  

#### Downloading & Setting Up the Toolkit
```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:/path/to/sratoolkit/bin
cd sratoolkit.3.1.1-ubuntu64/
source ~/.bashrc
```

#### Processing RNA-seq Sample (SRR5647735)
```sh
fastq-dump --split-3 -X 10000 SRR5647735
fastp --cut_tail -i SRR5647735_1.fastq -I SRR5647735_2.fastq -o SRR5647735_1.trim.fq -O SRR5647735_2.trim.fq
```

### FastQC Results: 
- **Pros**: Good quality reads, stable GC content.  
- **Cons**: Issues with Per Base Sequence Content at ends, overrepresented sequences, and adapter content → **Trimming Required**  

---

## Trimming & Adapter Removal  

### **Creating Adapter File**
```sh
echo ">illumina" > adapter.fa
echo "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" >> adapter.fa
```

### **Testing Different Trimmomatic Settings**  
Multiple settings were tested to optimize trimming.  

##### **Setting 1**
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

##### **Setting 2**
```sh
java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 4 -phred33 \
  SRR5647735_1.fastq SRR5647735_2.fastq \
  SRR5647735_1.trim.fq SRR5647735_1.unpaired.fq \
  SRR5647735_2.trim.fq SRR5647735_2.unpaired.fq \
  ILLUMINACLIP:adapter.fa:2:30:10 \
  SLIDINGWINDOW:4:20 \
  TRAILING:20 
```

##### **Setting 3**
```sh
java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 4 -phred33 \
  SRR5647735_1.fastq SRR5647735_2.fastq \
  SRR5647735_1.trim3.fq SRR5647735_1.unpaired3.fq \
  SRR5647735_2.trim3.fq SRR5647735_2.unpaired3.fq \
  TRAILING:20
```

##### **Setting 4 (Best Outcome)**
```sh
java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 4 -phred33 \
  SRR5647735_1.fastq SRR5647735_2.fastq \
  SRR5647735_1.trim4.fq SRR5647735_1.unpaired4.fq \
  SRR5647735_2.trim4.fq SRR5647735_2.unpaired4.fq \
  CROP:6 \
  HEADCROP:142
```

**Setting 4** showed the best outcome, resolving overrepresented sequences and adapter content, which were observed in other settings also **but Setting 4 also fixed the Per Base Sequence Content irregularities**, making it the final choice for trimming.

---

## Alignment to Reference Genome  

#### Setting Up Environment & Indexing Reference Genome
```sh
# Drosophila Genome collected from FlyBase database  
# URL: Index of genomes/Drosophila_melanogaster/current/fasta/  
# The dmel-all-chromosome-r6.63.fasta.gz file was downloaded.

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

#### **Installing & Building Kraken2 Database**
```sh
**Project Directory DIR:** `/home/hlnrajes/PROJECT/kraken/kraken_database`

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

**Final Observations:**  
- Trimming improved read quality significantly.  
- Alignment successfully mapped reads to the reference genome.  
- Kraken2 classified species from unmapped reads.  
