# Methods  

## Collecting Sample Reads  
To begin the analysis, I downloaded RNA-seq data from the **SRR5647735 sample**, which is part of the **Kapun et al. (2020)** study. This dataset was selected for its relatively high trypanosomatid presence.  
**Code for this step is available in `Commands_Run.md`.**  

## Quality Assessment of Raw Reads  
To evaluate sequencing quality, I ran **FastQC** on the unprocessed sample reads.  

### **Initial FastQC Results (Untrimmed Reads)**  
- **Read Lengths**: 35 - 151  
- **Total Sequences**: 34,117,123  
- **Cons**:  
  - Per Base Sequence Content (Pbsc) irregularities observed at the sequence ends  
  - Presence of **overrepresented sequences**  
  - **Adapter contamination detected**  
- **Pros**:  
  - High-quality reads  
  - Stable GC content  

## Need for Trimming  
Given the identified adapter content and sequence irregularities, **Trimmomatic** was used to refine the dataset. I tested **four different trimming settings**, followed by **FastQC analysis** to assess improvements.  
**Code for trimming steps is available in `Commands_Run.md`.**  

## Trimmomatic: Testing Four Settings  
Each trimming setting was evaluated for **read length, sequence integrity, and GC content stability**.  

### **FastQC Results After Trimming**  
#### **Setting 1**  
- **Read Lengths**: 36 - 151  
- **Total Sequences**: 9,040  
- **Cons**:  
  - **GC content graph inconsistent** (as seen in last meeting)  
  - **Pbsc irregularities improved but still present**  
- **Pros**:  
  - Overrepresented sequences **successfully removed**  
  - Adapter contamination **resolved**  

#### **Setting 2**  
- **Read Lengths**: 3 - 151  
- **Total Sequences**: 9,836  
- **Cons**:  
  - **GC content graph showed more peaks**, lacking smooth distribution  
  - **Pbsc issue persists**  
- **Pros**:  
  - Overrepresented sequences **removed**  
  - Adapter contamination **resolved**  

#### **Setting 3**  
- **Read Lengths**: 37 - 151  
- **Total Sequences**: 9,992  
- **Cons**:  
  - GC content graph **less erratic**, but Pbsc still present  
- **Pros**:  
  - Overrepresented sequences **removed**  
  - Adapter contamination **resolved**  

#### **Setting 4 (Final Selection)**  
- **Read Lengths**: 32 - 142  
- **Total Sequences**: 10,000  
- **Cons**:  
  - **Marginal duplicate sequence content observed**  
  - **Lowest average read length compared to other settings**  
- **Pros**:  
  - **Pbsc issue fully resolved**  
  - **Highest number of total sequences among all settings**  

### **Choice of Trimming Strategy**  
Setting **4 was chosen**, as it provided the **best resolution of Pbsc issues** while retaining a high number of quality reads.  

---

## Alignment to Drosophila Genome  
Using the **trimmed reads from Setting 4**, I aligned them against the *Drosophila melanogaster* genome.  

### **Reference Genome Selection**  
The reference genome was obtained from **FlyBase**:  
- **Source**: [Index of genomes/Drosophila_melanogaster/current/fasta/](https://flybase.org)  
- **Downloaded File**: `dmel-all-chromosome-r6.63.fasta.gz`  

### **BWA Alignment**  
The trimmed reads were aligned to the *Drosophila melanogaster* genome using **BWA-MEM**.  
**Code for alignment steps is available in `Commands_Run.md`.**  

---

## Extracting Unmapped Reads  
After alignment, **Samtools** was used to extract unmapped reads for further classification.  
**Code for extracting unmapped reads is available in `Commands_Run.md`.**  

---

## Taxonomic Classification with Kraken2  
The extracted unmapped reads were analyzed using **Kraken2**, a metagenomic classification tool.  
**Code for Kraken2 analysis is available in `Commands_Run.md`.**  

### **Extracting Species-Level Classification**  
The classified species data was saved in **species_classification.txt**, which was further analyzed in R. 
**Code for data analysis is available in `Commands_Run.md`.**  

---

## Visualization Using R  
The **species_classification.txt** file was loaded into **R**, where a graph was generated to depict the relative percentages of identified species.  

Further analysis will refine **graph formatting**, optimize taxonomic classification methods, and compare findings against published datasets.  

---
### 🎯 **Next Steps** 
- Interpreting Kraken output
- Expanding taxonomic classification  
- Refining visualization techniques  
- Investigating additional classification tools 
