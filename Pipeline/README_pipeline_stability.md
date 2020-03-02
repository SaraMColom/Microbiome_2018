# Chapter 3 Microbiome project data analysis

[Table of Contents](#title1)

[Important files](#title1)

[Overview of analysis pipeline](#title1)

[Experimental design](#title1)

[Table sample number by species and treatment](#title1)

[Table number of maternal line per species](#table-number-of-maternal-line-per-species)

[Code for running analysis in Mothur](#code-for-running-analysis)

- [OTU based analysis](#otu-based-analysis)

- [Rarefaction curve](#title)

- [Beta Diversity](#beta-diversity)


# Important files

**FastaQ files:** All files ending in 'fastaq.gz' include the forward and reverse DBA reads of the 16SU from a subset samplig of soil taken from the rhizosphere of plants in the presence and absence of competition. In the file prefix, the number indicates field ID and the letter 'H' and 'P' indicate species. 

**Taxonomy file:** This file ('silva.full_v132.tax') was made following Pat Schloss' tutorial here: http://blog.mothur.org/2018/01/10/SILVAv132referencefiles/. This allows us to map the sequences to a specific taxonomy, 'Archaea', 'Bacteria' and 'Eukaryota'.  

<pre><code>
list.seqs(fasta=silva.nr_v132.pcr.align) # list sequences in template file
get.seqs(taxonomy=silva.nr_v132.tax, accnos=current) # select those reads from the taxonomy file
</code></pre>


**SILVA reference file**This file ('silva.nr_v132.pcr.align) was made following Pat Schloss' tutorial (linke above). This contains the reference alighnment data for bacteria, archaea and eukaryota. Note*, I ran the 'running pcr.seqs' function. *Note* I had an error using the reference silva fasta file I created because not all of the reference sequences had a taxanomy assaigned to it. Therefore, I removed the fasta sequences with no taxanomy and used this 'reduced' silva file for my analysis ('Silva.v4.Red.fasta'). To create this fasta file I ran some code in R. 

Specifics on how I obtained the SILVA and tax files are found here: https://github.com/SaraMColom/MicrobiomeProject/blob/master/README_SILVA.md.

# Overview of analysis pipeline
![](https://raw.githubusercontent.com/SaraMColom/Microbiome_2018/master/Pipeline/Chapter3Pipeline.jpg)

# Experimental design
We took samples from focal host plant, _I. purpurea_, in the presence and absence of competition from its sister species _I. hederacea_. We also sampled the competitor _I. hederacea_ plant. We sampled a total of 183 plants, 27 alone plants and a total of 32 unique competitive combination types. Note, we had 29 combination types with at least 2 biological replicates, and three combinations with a single representation. Each mate 24rnal line of _I. purpurea_ was replicated 513 times across competition treatment (except for a single maternal line), and 24 times within alone treatment.

Below are a few tables describing our sample sizes by species and treatment, and the number of maternal lines and combination types per species.

### Table sample number by species and treatment


| Species          | Treatment       | N      |
|--------------    |-------------    |----    |
| I. purpurea      | Alone           | 27     |
| I. purpurea      | Competition     | 78     |
| I. hederacea     | Competition     | 78     |


### Table number of maternal line per species

| Species          | Number of ML     |
|--------------    |--------------    |
| I. purpurea      |      10          |
| I. hederacea     |       5          |


# Code for running analysis

*Step 1:* Unzip the fasta.gz files in the directory, this step is in linux


```TestRunMothur scolom$ for i in `ls`; do gunzip $i; done```


# Mothur pipeline

The following code should run from within Mothur in the Flux system. Run the stability file ('stability.batch') to automate the step.

### Create stability file in mothur


```make.file(inputdir=., type=fastq, prefix=stability)```

Change the name of the file from stability.files to whatever suits your study

make.contigs(file=stability.files, processors=8)

This implementation of the command will remove any sequences with ambiguous bases and anything longer than 275 bp. 


```screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275)```

## Processing improved sequences

### Remove duplicates

```unique.seqs()```

This will generate a file called stability.trim.contigs.good.count_table.

```count.seqs(name=current, group=current)```

### Align our sequences to the reference alignment.


```align.seqs(fasta=current, reference=silva.nr_v132.pcr.align)```
```screen.seqs(fasta=current, count=current, start=1968, end=11550, maxhomop=8)```


### Filter the sequences to remove the overhangs at both ends.Remove gaps

`filter.seqs(fasta=current, vertical=T, trump=.)`

### Remove redundancy

```unique.seqs(fasta=current, count=current)```


### Denoise with 'pre.cluster'
This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence.


```pre.cluster(fasta=current, count=current, diffs=2)```

### Search & remove chimeras


```chimera.uchime(fasta=current, count=current, dereplicate=t)```


Also remove them from the fasta files.


```remove.seqs(fasta=current, accnos=current)```

### Remove unwanted DNA sequences

More cleanup, remove _undesirables_, e.g., Archaea, chlorpolast and mitochondria DNA.

First _classify_ the sequences using our reference file.


```classify.seqs(fasta=current, count=current, reference=, taxonomy=, cutoff=80)```

Now that its classified, remove the undesirables.


```remove.lineage(fasta=current, count=current, taxonomy=current,taxon=ChloroplastMitochondriaunknownArchaeaEukaryota)```


I did *not* run the code below because I did not have a Mock community.


```remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock) ```# Not run

## OTUS

In this approach, we use the taxonomic information to split the sequences into bins and then cluster within each bin.
The advantage of the cluster.split approach is that it should be faster,use less memory, and can be run on multiple processors.


```cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify,taxlevel=4, cutoff=0.15)```


Next we want to know how many sequences are in each OTU from each group and we can do this using the make.shared command.
Here we tell mothur that we're really only interested in the 0.03 cutoff level
For some analyses you may desire to bin your sequences in to phylotypes according to their taxonomic classification. We can do this using the phylotype command.
We also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the classify.otu command.
The cutoff numbering is a bit different for phylotype compared to cluster/cluster.split. Here you see 1 through 6 listed; these correspond to Genus through Kingdom levels, respectively.
So if you want the genuslevel shared file we'll run the'make.shared' function.

Lastly, We also want to know who these OTUs are and can run classify.otu on our phylotypes.



```make.shared(list=current, count=current, label=0.03)```
```classify.otu(list=current, count=current, taxonomy=current, label=0.03)```
```phylotype(taxonomy=current)```
```make.shared(list=current, count=current, label=1)```
```classify.otu(list=current, count=current, taxonomy=current, label=1)```


# OTU based analysis

## Rarefaction curve

Allows us to analyze the alpha diversity of our samples. The 'rarefaction.single' function in the mothur program can be used to do this. The 'summary.single' function can be used to sample the coverage, number of observed OTU's and inverse Simpson diversity estimate. 


```rarefaction.single(shared=stability.opti_mcc.shared, calc=sobs, freq=100)```
```summary.single(shared=stability.opti_mcc.shared, calc=nseqscoveragesobsinvsimpson,subsample=T)```

The R package vegan can also give us alpha diversity with the function _multipart()_

## Beta Diversity
Beta diversity can be estimated in mothur with the 'dist.shared' function and the code below can be implemented for this.


```dist.shared(shared=stability.opti_mcc.shared, calc=thetaycjclass, subsample=t)```


Similarly, the Beta Diversity can be estimated with the R package vegan using the function: _betadisper_ based on the BrayCurtis distances.

## Species richness
For species richness see the _summary.single()_ function.

# Runing job on Great Lakes

Log into great lakes thru terminal: 

`ssh -l UserName greatlakes.arc-ts.umich.edu`
`/scratch/lsa_root/lsa/UserName/microbiome`


# ANOVA in R

Once we have our indeces of microbial community variables (e.g., diversity and richness), we will compare how these variables vary by: (1) Species, (2) Treatment, (3) Maternal line.

We will use the _lm_ and _anova_ functions in base R to perform simple fixed linear models, and the _lmer_ function from the lmerTest package to perform more complex linear mixed models.


# Rscripts 


