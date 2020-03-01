# README for the SILVA v132 reference files

Jan 10, 2018 • PD Schloss • 9 min read 

_Modified on Sep 16, 2019, It took me a couple days to figure this out because some of the code initially provided did not grab the correct columns etc. Alterations to original docx are noted below._

The good people at SILVA have released a new version of the SILVA database. A little bit of tweaking is needed to get their files to be compatible with mothur. This README document describes the process that I used to generate the mothur-compatible reference files.

Curation of references
Getting the data in and out of the ARB database
This README file explains how we generated the silva reference files for use with mothur’s classify.seqs and align.seqs commands. I’ll assume that you have a functioning copy of arb installed on your computer. For this README we are using version 6.0. First we need to download the database and decompress it. From the command line we do the following:

<pre><code>
wget -N https://www.arb-silva.de/fileadmin/arb_web_db/release_132/ARB_files/SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz
gunzip SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz
arb SILVA_132_SSURef_NR99_13_12_17_opt.arb
</code></pre>

This will launch us into the arb environment with the ‘‘Ref NR 99’’ database opened. This database has 695,171 sequences within it that are not more than 99% similar to each other. The release notes for this database as well as the idea behind the non-redundant database are available from the silva website. Within arb do the following:

Click the search button
Set the first search field to ‘ARB_color’ and set it to 1. Click on the equal sign until it indicates not equal (this removes low quality reads and chimeras)
Click ‘Search’. This yielded 629,211 hits
Click the “Mark Listed Unmark Rest” button
Close the “Search and Query” box
Now click on File->export->export to external format
In this box the Export option should be set to marked, Filter to none, and Compression should be set to no.
In the field for Choose an output file name make sure the path has you in the correct working directory and enter silva.full_v132.fasta`.
Select a format: fasta_mothur.eft. This is a custom formatting file that I have created that includes the sequences accession number and it’s taxonomy across the top line. To create one for you will need to create fasta_mothur.eft in the /opt/local/share/arb/lib/export/ folder with the following:

<pre><code>
SUFFIX          fasta    
BEGIN    
>*(acc).*(name)\t*(align_ident_slv)\t*(tax_slv);    
*(|export_sequence)    
</code></pre>

Save this as silva.full_v132.fasta
You can now quit arb.

## Modification 
*Instead* of copoying the code above, I created the file mannually in the directory I am working in. Specifically, I opened a blank file in a text editor called 'fasta_mothur.eft'and saved it in the /opt/local/share/arb/lib/export/  folder. This was because I was unsure how to authenticate such a file. 

Screening the sequences
Now we need to screen the sequences for those that span the 27f and 1492r primer region, have 5 or fewer ambiguous base calls, and that are unique. We’ll also extract the taxonomic information from the header line. Run the following commands from a bash terminal:

<pre><code>
mothur "#screen.seqs(fasta=silva.full_v132.fasta, start=1044, end=43116, maxambig=5, processors=8);
        pcr.seqs(start=1044, end=43116, keepdots=T);
        degap.seqs();
        unique.seqs();"

#identify the unique sequences without regard to their alignment
grep ">" silva.full_v132.good.pcr.ng.unique.fasta | cut -f 1 | cut -c 2- > silva.full_v132.good.pcr.ng.unique.accnos

#get the unique sequences without regard to their alignment
mothur "#get.seqs(fasta=silva.full_v132.good.pcr.fasta, accnos=silva.full_v132.good.pcr.ng.unique.accnos)"

#generate alignment file
mv silva.full_v132.good.pcr.pick.fasta silva.nr_v132.align

#generate taxonomy file
grep '>' silva.nr_v132.align | cut -f1,3 | cut -f2 -d'>' > silva.nr_v132.full
</code></pre>

The mothur commands above do several things. First the screen.seqs command removes sequences that are not full length and have more than 5 ambiguous base calls. Note: this will remove a number of Archaea since the ARB RN reference database lets in shorter (>900 bp) archaeal 16S rRNA gene sequences. Second, pcr.seqs convert any base calls that occur before position 1044 and after 43116 to . to make them only span the region between the 27f and 1492r priming sites. Finally, it is possible that weird things happen in the alignments and so we unalign the sequences (degap.seqs) and identify the unique sequences (unique.seqs). We then convert the resulting fasta file into an accnos file so that we can go back into mothur and pull out the unique sequences from the aligned file (get.seqs).

## Modification!
*Instead* of runing the line above, 'grep '>' silva.nr_v132.align | cut -f1,3 | cut -f2 -d'>' > silva.nr_v132.full', I ran the following below:

grep ">" silva.nr_v132.align | awk '{print $1,"\t",$3,$4}' | cut -f2 -d ">" > silva.nr_v132.full

I ran the line below because the original code line provided did not copy over the column of interest containing the taxonomic labels referred to as 'taxlabels' in the R scripting. Note this introduces a space before the taxonomic lables and therefore additional external modification was done to remove the spacing, see below.

## External modification
I opened the 'silva.nr_v132.full' file in a text editor and did find and replace of spaces. This is very important, otherwise the taxonomic labels will not match with the names in the taxonomic text document later downloaded.

#Formatting the taxonomy files
Now we want to make sure the taxonomy file is properly formatted for use with mothur. First we want to grab the SILVA taxa mapping file by running the following in bash:

<pre><code>
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_132.txt
</code></pre>

Thanks to Eric Collins at the University of Alaska Fairbanks, we have some nice R code to map all of the taxa names to the six Linnean levels (kingdom, phylum, class, order, family, and genus). We’ll run the following code from within R:


<pre><code>
# Set directory first

map.in <- read.table("tax_slv_ssu_132.txt",header=F,sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")
map.in <- rbind(map.in, c("Bacteria;RsaHf231;", "phylum")) #wasn't in tax_slv_ssu_132.txt

#fix Escherichia nonsense
map.in$taxlevel[which(map.in$taxlabel=="Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;")] <- "genus"

taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus")
taxabb <- c("ro","do","mc","pk","ki","bk","ik","pp","ph","bp","ip","pc","cl","bc","ic","po","or","bo","pf","fa","bf","ge")
tax.mat <- matrix(data="",nrow=nrow(map.in),ncol=length(taxlevels))
tax.mat[,1] <- "root"
colnames(tax.mat) <- taxlevels

outlevels <- c("domain","phylum","class","order","family","genus")

for(i in 1:nrow(map.in)) {
	taxname <- unlist(strsplit(as.character(map.in[i,1]), split=';'))
	#print(taxname);

	while ( length(taxname) > 0) {
		#regex to look for exact match

		tax.exp <- paste(paste(taxname,collapse=";"),";",sep="")
		tax.match <- match(tax.exp,map.in$taxlabel)
		tax.mat[i,map.in[tax.match,2]] <- tail(taxname,1)
		taxname <- head(taxname,-1)
	}
}

for(i in 1:nrow(tax.mat)) {
	#this fills in the empty gaps by using the closest higher taxonomic level appended with an abbreviation for the current taxonomic level
	#if you don't want this behavior, cut it out
	for(j in 1:ncol(tax.mat)) {
		if(tax.mat[i,j] < 0) { tax.mat[i,j] <- paste(tmptax,taxabb[j],sep="_")}
		else { tmptax <- tax.mat[i,j]}
	}

	#this maps the new name to the input taxonomic levels
	map.in[i,"taxout"] <- paste(paste(tax.mat[i,outlevels],collapse=";"),";",sep="")
}

# replace spaces with underscores
map.in$taxout <- gsub(" ","_",map.in$taxout)

# bring in the old taxonomic levels from SILVA and remap them using the new levels
tax.in <- read.table("silva.nr_v132.full",header=F,stringsAsFactors=F,sep="\t")
colnames(tax.in) <- c("taxid","taxlabel")

# Following line corrects the Bacteria;Bacteroidetes;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Polaribacter;Polaribacter; problem
tax.in$taxlabel <- gsub("Polaribacter;Polaribacter;", "Polaribacter;", tax.in$taxlabel)
tax.in$taxlabel <- gsub(";[[:space:]]+$", ";", tax.in$taxlabel)

tax.in$id <- 1:nrow(tax.in)

tax.write <- merge(tax.in,map.in,all.x=T,sort=F)
tax.write <- tax.write[order(tax.write$id),]


#we want to see whether everything has 6 taxonomic level (kingdom to genus)
getDepth <- function(taxonString){
	initial <- nchar(taxonString)
	removed <- nchar(gsub(";", "", taxonString))
	return(initial-removed)
}

depth <- getDepth(tax.write$taxout)
summary(depth) #should all be 6 and there should be no NAs <- I got 4982 NA's, note this is due to some of the tax labels from the 'txt' file I downloaded above. I later do 'na.omit' to exclude those that do not match to the SILVA file I work with.
bacteria <- grepl("Bacteria;", tax.write$taxout)
archaea <- grepl("Archaea;", tax.write$taxout)
eukarya <- grepl("Eukaryota;", tax.write$taxout)

tax.write[depth > 6 & bacteria,] #good to go
tax.write[depth > 6 & archaea,]  #good to go
tax.write[depth > 6 & eukarya,]  #good to go

write.table(tax.write[,c("taxid","taxout")],file="silva.full_v132.tax",sep="\t",row.names=F,quote=F,col.names=F)
</code></pre>

#Building the SEED references
The first thing to note is that SILVA does not release their SEED; it is private. By screening through the ARB databases we can attempt to recreate it. Our previous publications show that classify.seqs with the recreated SEED does an excellent job of realigning sequences to look like they would if you used SINA and the true SEED. Now we want to try to figure out which sequences are part of the seed. Earlier, when we exported the sequences from ARB, we included the align_ident_slv field from the database in our output. Let’s generate an accnos file that contains the names of the sequences with 100% to the SEED database and then use mothur to generate SEED fasta and taxonomy files. While we’re at it we’ll also generate the nr_132 taxonomy file as well. The following code will be run from within a bash terminal:

<pre><code>
grep ">" silva.nr_v132.align | cut -f 1,2 | grep "\t100" | cut -f 1 | cut -c 2- > silva.seed_v132.accnos
mothur "#get.seqs(fasta=silva.nr_v132.align, taxonomy=silva.full_v132.tax, accnos=silva.seed_v132.accnos)"
mv silva.nr_v132.pick.align silva.seed_v132.align
mv silva.full_v132.pick.tax silva.seed_v132.tax

mothur "#get.seqs(taxonomy=silva.full_v132.tax, accnos=silva.full_v132.good.pcr.ng.unique.accnos)"
mv silva.full_v132.pick.tax silva.nr_v132.tax
Taxonomic representation
Let’s look to see how many different taxa we have for each taxonomic level within the silva.nr_v132.tax, silva.seed_v132.tax. To do this we’ll run the following in R:

getNumTaxaNames <- function(file, kingdom){
  taxonomy <- read.table(file=file, row.names=1)
  sub.tax <- as.character(taxonomy[grepl(kingdom, taxonomy[,1]),])

  phyla <- as.vector(levels(as.factor(gsub("[^;]*;([^;]*;).*", "\\1", sub.tax))))
  phyla <- sum(!grepl(kingdom, phyla))

  class <- as.vector(levels(as.factor(gsub("[^;]*;[^;]*;([^;]*;).*", "\\1", sub.tax))))
  class <- sum(!grepl(kingdom, class))

  order <- as.vector(levels(as.factor(gsub("[^;]*;[^;]*;[^;]*;([^;]*;).*", "\\1", sub.tax))))
  order <- sum(!grepl(kingdom, order))

  family <- as.vector(levels(as.factor(gsub("[^;]*;[^;]*;[^;]*;[^;]*;([^;]*;).*", "\\1", sub.tax))))
  family <- sum(!grepl(kingdom, family))

  genus <- as.vector(levels(as.factor(gsub("[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;([^;]*;).*", "\\1", sub.tax))))
  genus <- sum(!grepl(kingdom, genus))

  n.seqs <- length(sub.tax)
  return(c(phyla=phyla, class=class, order=order, family=family, genus=genus, n.seqs=n.seqs))
}

kingdoms <- c("Bacteria", "Archaea", "Eukaryota")
tax.levels <- c("phyla", "class", "order", "family", "genus", "n.seqs")

nr.file <- "silva.nr_v132.tax"
nr.matrix <- matrix(rep(0,18), nrow=3)
nr.matrix[1,] <- getNumTaxaNames(nr.file, kingdoms[1])
nr.matrix[2,] <- getNumTaxaNames(nr.file, kingdoms[2])
nr.matrix[3,] <- getNumTaxaNames(nr.file, kingdoms[3])
rownames(nr.matrix) <- kingdoms
colnames(nr.matrix) <- tax.levels
nr.matrix
#          phyla class order family genus n.seqs
#Bacteria     80   204   580   1052  3971 188247
#Archaea      11    30    52     85   210   4626
#Eukaryota    93   240   648    923  3018  20246

# My nr.matrix
#         phyla class order family genus n.seqs
#Bacteria     70   174   480    847  3036  54226
#Archaea      11    26    42     69   165   1842
#Eukaryota    71   176   489    710  2083   6140

seed.file <- "silva.seed_v132.tax"
seed.matrix <- matrix(rep(0,18), nrow=3)
seed.matrix[1,] <- getNumTaxaNames(seed.file, kingdoms[1])
seed.matrix[2,] <- getNumTaxaNames(seed.file, kingdoms[2])
seed.matrix[3,] <- getNumTaxaNames(seed.file, kingdoms[3])
rownames(seed.matrix) <- kingdoms
colnames(seed.matrix) <- tax.levels
seed.matrix
#          phyla class order family genus n.seqs
#Bacteria     50   110   301    530  1436   8517
#Archaea       7    15    26     39    62    147
#Eukaryota    41   100   287    478  1040   2516

# My seed.matrix 
#          phyla class order family genus n.seqs
#Bacteria     48   107   253    438   928   3099
#Archaea       6    15    24     33    51    109
#Eukaryota    34    68   138    158   233    332

seed.matrix / nr.matrix
#              phyla     class     order    family     genus     n.seqs
#Bacteria  0.6250000 0.5392157 0.5189655 0.5038023 0.3616218 0.04524375
#Archaea   0.6363636 0.5000000 0.5000000 0.4588235 0.2952381 0.03177691
#Eukaryota 0.4408602 0.4166667 0.4429012 0.5178765 0.3445991 0.12427146

# My seed.matrix/nr.matrix
#             phyla     class     order    family     genus     n.seqs
#Bacteria  0.6857143 0.6149425 0.5270833 0.5171192 0.3056653 0.05714971
#Archaea   0.5454545 0.5769231 0.5714286 0.4782609 0.3090909 0.05917481
#Eukaryota 0.4788732 0.3863636 0.2822086 0.2225352 0.1118579 0.05407166
</code></pre>

The Archaea take a beating and recall they lost a bunch of sequences in the initial steps since many of the arachaeal sequences in SILVA are between 900 and 1200 nt long. If you are interested in analyzing the Archaea and the Eukaryota, I would suggest duplicating my efforts here but modify the screen.seqs and pcr.seqs steps to target your region of interest.

Finally, we want to compress the resulting alignment and this README file into the full length and SEED archives using commands in the bash terminal:

<pre><code>
tar cvzf silva.nr_v132.tgz silva.nr_v132.tax silva.nr_v132.align README.md
tar cvzf silva.seed_v132.tgz silva.seed_v132.tax silva.seed_v132.align README.md
</code></pre>

#Application
So… which to use for what application? If you have the RAM, I’d suggest using silva.nr_v132.align in align.seqs. It took about 10 minutes to read in the database file and a minute or so to align a 1000 full-length sequences. Here is an example workflow for use within mothur that will get you the V4 region of the 16S rRNA gene:

<pre><code>
mothur "#pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F, processors=8);
        unique.seqs()"
</code></pre>

This will get you 139,321 unique sequences to then align against (meh.). Other tricks to consider would be to use get.lineage to pull out the reference sequences that are from the Bacteria, this will probably only reduce the size of the database by ~10%. You could also try using filter.seqs with vertical=T; however, that might be problematic if there are insertions in your sequences (can’t know a priori). It’s likely that you can just use the silva.seed_v132.align reference for aligning. For classifying sequences, I would strongly recommend using the silva.nr_v132.align and silva.nr_v132.tax references after running pcr.seqs on silva.nr_v132.align. I probably wouldn’t advise using unique.seqs on the output.

**Important** I used the 'pcr.seqs' command as suggested right above, on the 'silva.nr_v132.align' and created a new file called 'silva.nr_v132.align.pcr.align'. Then I made the latter a fasta file on a MAC runing the following code in the terminal:

<pre><code>
mv silva.nr_v132.pcr.align silva.v4.fasta
</code></pre>

When I tried out the 'stability.batch' file by Pat Schloss using the files I created from the procedure above, I encountered an error, 'xxx is in your taxonomy file and is not in your template file. Please correct.' I used a Mothur Forum to trouble shoot this linked here: https://github.com/mothur/mothur/issues/505. It appeared that the issue was because when I "ran pcr.seqs some of the references were scrapped. You can select the "good" sequences from the template using the get.seqs command." I followed the suggested codes below, which made a new file called 'silva.full_v132.pick.tax', and I used that one instead and it worked.

<pre><code>
list.seqs(fasta=silva.nr_v132.pcr.align) # list sequences in template file
get.seqs(taxonomy=silva.nr_v132.tax, accnos=current) # select those reads from the taxonomy file
</code></pre>



mothur > list.seqs(fasta=silva.nr_v132.pcr.align) - list sequences in template file
mothur > get.seqs(taxonomy=silva.nr_v132.tax, accnos=current) - select those reads from the taxonomy file

Legalese
If you are going to use the files generated in this README, you should be aware of SILVA’s dual use license. We’ll leave it to you to work out the details.

# Note the following:

## RawData folder
- Contains the files used and created throughout the process of creating the SILVA reference files.

## CleanData
- Contains the zipped files produced at the end of the entire process above.

## R script
- Contains the R code ultimately used to format the taxonomy files and to evaluate the taxonomic representation.
