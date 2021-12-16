### Hypotheses under investigation: 
* I am first interested in the differences in alpha diversity between healthy and sick individuals. 

* I am asking the question: Is alpha diversity between healthy and sick sea stars constant throughout the sites sampled?

* Further, are there specific bacterial taxa that are present in only sick or only healthy sites/sea stars? How much of the variation is due to site differences?

### In order to test these questions, I will use the full dataset and perform the steps we used in class. I will calculate alpha diversity in the samples and compare that between sites.

### First, I had to activate qiime, so i used this code every time I started a new session
```{bash}
conda activate qiime2-2021.8
```

### I used a temporary directory while running this because the directory on the server is too small.
```{bash}
export TMPDIR="/data/project_data/16S/tmptmpdir"
echo $TMPDIR 
```

### In order to create the paired-end output file from the metadata, I need to specify the file type, provide the path to the metadata file, specify the format of the data, and create an output file.
```{bash}
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /data/project_data/16S/pyc_manifest \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux-paired-end_full.qza
```

### Next I visualized the quality of the data to determine how much to trim. I had to create a .qzv file in order to visualize it and view it in qiime.
```{bash}
qiime demux summarize \
  --i-data demux-paired-end.qza \         
  --o-visualization demux-pyc-sub.qzv      
```

### Using DADA2, I was able to 'denoise' the data to whatever parameters I thought best. 

#### On forward reads, I trimmed at base 12 and 278, and on reverse reads I trimmed at 0 and 256. 

### I chose these because I did not want to cut off too much of the data if I didn't have to, and I just cut off data that had a score much lower than 30.

```{bash}
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-n-threads 4 \
  --p-trim-left-f 12 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 278 \
  --p-trunc-len-r 256 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```

### I created visualization files based on the previous .qza files generated.
```{bash}
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file pyc_manifest

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

```

### I built a phylogenetic tree among identified taxa in order to perform diversity metrics.
```{bash}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

### After I made my taxa-bar-plots.qzv file, I was able to visualize the taxa on various different levels and comparing between different individuals!
```{bash}
qiime taxa barplot \
> --i-table table.qza \
> --i-taxonomy taxonomy.qza \
> --m-metadata-file /data/project_data/16S/pyc_manifest \
> --o-visualization taxa-bar-plots.qzv
```

#### Just looking at the genus level in healthy vs sick sea stars, I notice that in many of the sick individuals, bacteroidales is very common. This genus is not as common in the healthy individuals. 
#### Spirochaetaceae is very prevalent in the healthy sea stars, and not seen much at all in the sick.
#### Also, vibrionaceae is very common in the sick stars as well as at sick sites.

### I used ancom to test for differences in abundance and ran it in a screen.
```{bash}
screen
cd /data/project_data/16S/mprun/taxa
conda activate qiime2-2021.8
export TMPDIR="/data/project_data/16S/tmptmpdir"

qiime composition add-pseudocount \
  --i-table /data/project_data/16S/mprun/table.qza \
  --o-composition-table comp-table.qza
 
qiime composition ancom \
  --i-table comp-table.qza \
  --m-metadata-file /data/project_data/16S/pyc_manifest \
  --m-metadata-column site-animal-health \
  --o-visualization ancom-site-animal-health.qzv 
```