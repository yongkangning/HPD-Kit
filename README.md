# HPD-Kit
HPD-Kit(Henbio Pathogen Detection Toolkit) is a robust, plug-and-play, and reproducible pathogen analysis tool. It is implemented using Singularity container technology and the Nextflow workflow framework.

## Installing Singularity on Your Linux System

To install Singularity, run the following command:

```bash
conda install singularity=3.8.6
```

## Downloading Image and Database Files

- Download the HPD-Kit-v2.0.img image and the accompanying database from the following link:
  (https://pan.baidu.com/s/1nH5Uv_BoeBpGalr7bllt0A?pwd=68cd)
- Unzip the database file:
  ```bash
  unzip virus.zip
  ```

## Viewing Help

To view the help information for HPD-Kit, execute the following command:

```bash
singularity exec HPD-Kit-v2.0.img nextflow run /opt/tools/HPD-Kit/main.nf --help
```

## Testing

### Built-in Dataset Test

HPD-Kit includes a built-in dataset for testing. To test if HPD-Kit is functioning correctly, run:

```bash
singularity exec HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --sample_id=SRR3214089 \
    --fastq1=/opt/tools/HPD-Kit/test_data/SRR3214089.fastq.gz \
    --pathogenTypes=virus \
    --outDir=result
```

You should find the analysis results in the `result/identified_pathogens` directory:
- `SRR3214089.pathogens.final.xls`
- `SRR3214089.virus.final.xls`

When `pathogenTypes` specifies only one type of pathogen, the contents of these two files are the same.

### Paired-end Sequencing Data

If you have a paired-end sequencing sample (e.g., `SRR8476435_1.fastq.gz`, `SRR8476435_2.fastq.gz`) stored in `/home/user/data`, and the HPD-Kit database and image have been downloaded and extracted to `/home/user/data/HPD-Kit`, you can analyze the sample for four types of pathogens (bacteria, virus, fungi, parasite) using the following command:

```bash
singularity exec /home/user/data/HPD-Kit/HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --sample_id=SRR8476435 \
    --fastq1=/home/user/data/SRR8476435_1.fastq.gz \
    --fastq2=/home/user/data/SRR8476435_2.fastq.gz \
    --pathogenTypes=bacteria,virus,fungi,parasite \
    --outDir=result
```

By default, dehosting is not performed. If you need to remove the host, you can use `bowtie2` or `bbduk`. For `bowtie2`, create an index file with:

```bash
bowtie2-build -f /home/user/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --threads 6 /home/user/data/bowtie2_index/Homo_sapiens
```

Then run:

```bash
singularity exec /home/user/data/HPD-Kit/HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --sample_id=SRR8476435 \
    --fastq1=/home/user/data/SRR8476435_1.fastq.gz \
    --fastq2=/home/user/data/SRR8476435_2.fastq.gz \
    --dehost_method=bowtie2 \
    --host_fasta=/home/user/data/bowtie2_index/Homo_sapiens \
    --pathogenTypes=bacteria,virus,fungi,parasite \
    --outDir=result
```

For `bbduk`, provide the host fasta genome file:

```bash
singularity exec /home/user/data/HPD-Kit/HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --sample_id=SRR8476435 \
    --fastq1=/home/user/data/SRR8476435_1.fastq.gz \
    --fastq2=/home/user/data/SRR8476435_2.fastq.gz \
    --dehost_method=bbduk \
    --host_fasta=/home/user/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --pathogenTypes=bacteria,virus,fungi,parasite \
    --outDir=result
```

### Batch Analysis

For multiple samples, create a `metadata.txt` describing the sample information:

```
SampleId,Read1,Read2,Host
SRR3214089,/data/SRR3214089.fastq.gz,,/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SRR3214092,/data/SRR3214092.fastq.gz,,/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SRR8476532,/data/SRR8476532_1.fastq.gz,/data/SRR8476532_2.fastq.gz,no
SRR8476535,/data/SRR8476535_1.fastq.gz,/data/SRR8476535_2.fastq.gz,no
```

Run the batch analysis with:

```bash
singularity exec /home/user/data/HPD-Kit/HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --metadataFile=/home/user/data/metadata.txt \
    --dehost_method=bbduk \
    --pathogenTypes=bacteria,virus,fungi,parasite \
    --outDir=result
```

### Comparison Analysis

If you have control samples, set up a comparison file `comparison_file.txt`:

```
SRR3214089     SRR3214092
SRR8476532     SRR8476535
```

Run the comparison analysis with:

```bash
singularity exec /home/user/data/HPD-Kit/HPD-Kit-v2.0.img \
  nextflow run /opt/tools/HPD-Kit/main.nf \
    --metadataFile=/home/user/data/metadata.txt \
    --comparisonFile=/home/user/data/comparison_file.txt \
    --dehost_method=bbduk \
    --pathogenTypes=bacteria,virus,fungi,parasite \
    --outDir=result
```

## Results Interpretation

Upon successful completion, you will find the following directories in the `result` folder:
- `identified_pathogens`: Final pathogen identification results.
- `kraken2`: Preliminary identification results from Kraken2.
- `bwt2`: BAM and BAI files from Bowtie2 alignments.
- `coverage_plot`: Coverage plot files (generated if `--plot_coverage=1` is set).

### Column Descriptions in `sampleId.pathogens.final.xls`

- `taxid`: Pathogen taxid.
- `scientific_name`: Scientific name of the pathogen.
- `assembly_accession`: Pathogen reference genome assembly version used.
- `accession_version`: Pathogen sequence names matched, up to 20 sequences displayed, separated by `|`.
- `pathogen_type`: Type of pathogen (bacteria, virus, fungi, parasite).
- `sample_total_reads`: Read count after host removal.
- `unique_reads`: Unique reads after Bowtie2 and BLAST alignment.
- `bwt2_reads_count`: Unique reads from Bowtie2 alignment.
- `k2_reads_count`: Read count identified by Kraken2.
- `k2_kmers`: Kmer count identified by Kraken2.
- `k2_unique_kmers`: Unique kmer count identified by Kraken2.
- `k2_relative_abundance`: Relative abundance of the pathogen in the sample as identified by Kraken2.
- `NPA`: Normalized Pathogen Abundance.
- `NPAS`: Normalized Pathogen Abundance Score.
- `NPA_adjusted`: Adjusted NPA when control samples are present.
- `NPAS_adjusted`: Adjusted NPAS when control samples are present.
- `base_coverage(%)`: Genomic base coverage percentage.
- `base_depth`: Sequencing depth.
- `genome_size`: Genome size of the pathogen.
- `sequences_count`: Number of sequences in the pathogen genome.
- `mapped_sequences_count`: Number of sequences mapped.
- `sequence_coverage(%)`: Genomic sequence coverage percentage.
- `avg_similarity(%)`: Average similarity percentage.
- `min_similarity(%)`: Minimum similarity percentage.
- `max_similarity(%)`: Maximum similarity percentage.
- `passed_filter`: Positive or negative judgment,

Contact Us

If you have any questions or require further assistance, feel free to contact us via email at zhangzn340@163.com. We are here to help!
