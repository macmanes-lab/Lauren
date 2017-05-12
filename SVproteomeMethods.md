#Seminal Vessicle Proteome


**Source Data: Protein Prospector RAW spectra data (DATAset S)**

***This is the information about dataset-S***

 

	This dataset consisted of 6 seminal vesicles: 3 WET & 3 DRY

	WET: 1W, W4, W5

	DRY: 1D, D3, D8


1D and 1W each have one dataset apiece of spectra data

D3, D8, W4,& W5 have two separate datasets apiece of spectra data, because two of the gel bands (Bands 2&3) were analyzed separately from from the other bands, due to their high protein density.

Therefore, when the data was imported for these samples, we merged the Bands 2&3 data for each sample with the other bands data for each sample.

example: SD323 (corresponds to DatasetS, sample D3, bands 2&3) and SD3oth (coresponds to DatasetS, sample D3, other bands), and these two datasets for sample D3 are merged into SD3 (DatasetS, sample D3)

----


**1) I imported this dataset into R, (R version 3.3.1 (2016-06-21))**


	setwd("/Users/lauren/Desktop/proteome/SPECTRArawDatasetS")

	S1D<-read.csv("S1Draw.csv",header=T)
	S1W<-read.csv("S1Wraw.csv",header=T)
	SD323<-read.csv("SD323raw.csv",header=T)
		SD3oth<-read.csv("SD3othraw.csv",header=T)
	SD823<-read.csv("SD823raw.csv",header=T)
	SD8oth<-read.csv("SD8othraw.csv",header=T)
	SW423<-read.csv("SW423raw.csv",header=T)
	SW4oth<-read.csv("SW4othraw.csv",header=T)
	SW523<-read.csv("SW523raw.csv",header=T)
	SW5oth<-read.csv("SW5othraw.csv",header=T)

	SD3<-merge(SD323,SD3oth,all=TRUE) 
	SD3$D3norm<-rowSums(cbind(SD3$D3raw23,SD3$D3rawOth),na.rm=T)
	SD3<-SD3[,c(1,4)]
	head(SD3)

	SD8<-merge(SD823,SD8oth,all=TRUE) 
	SD8$D8norm<-rowSums(cbind(SD8$D8raw23,SD8$D8rawOth),na.rm=T)
	head(SD8)
	SD8<-SD8[,c(1,4)]
	head(SD8)

	SW4<-merge(SW423,SW4oth,all=TRUE) 
	SW4$W4norm<-rowSums(cbind(SW4$W4raw23,SW4$W4rawOth),na.rm=T)
	head(SW4)
	SW4<-SW4[,c(1,4)]
	head(SW4)

	SW5<-merge(SW523,SW5oth,all=TRUE) 
	SW5$W5norm<-rowSums(cbind(SW5$W5raw23,SW5$W5rawOth),na.rm=T)
	head(SW5)
	SW5<-SW5[,c(1,4)]
	head(SW5)


	Smerged<-Reduce(function(x, y) merge(x, y, all=TRUE), 
                list(S1D,S1W,SD3,SD8,SW4,SW5))
	head(Smerged)

	Smerged[is.na(Smerged)]<-0
	head(Smerged)



	write.csv(Smerged,file="Smerged.csv")

***There are 1142 different contig acc#***

*The naked contigIDs list is called ContigIDsDatasetS.csv*

It is also named as :  ContigIDsDatasetS

--


**2) I retrieved all the sequences (from the P. eremicus amino acid database fasta file) corresponding to the 1142 contig IDs from ContigIDsDatasetS to make an unannotated seminal vesicle proteome**


In order to do so, I also uploaded the fasta file of the peromyscus eremicus database (amino acid seqs): peromyscus_eremicus.fasta


I used a program that Matt MacManes wrote to pull out amino acid sequences from the peromyscus eremicus database file that correspond to the ContigIDsDatasetS

*The program is filter.py*


	"""
	%prog some.fasta wanted-list.txt
	"""
	from Bio import SeqIO
	import sys

	wanted = [line.strip() for line in open(sys.argv[2])]
	seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')
	SeqIO.write((seq for seq in seqiter if seq.id in wanted), 	sys.stdout, "fasta")


*Execution of the program:*

	python filter.py peromyscus_eremicus.fasta ContigIDsDatasetS >> proteome.fasta


***The results file (proteome.fasta) contians the Contig IDs and Sequences for the proteome, this has been downloaded and it is the unannotated seminal vesicle proteome, which is available on Uniprot : SVproteome.fasta.***


--

**3) Next I annotated the Contig IDs from the proteome with Mus musculus to make an annotated Contig ID list.**


First I made a Mus Musculus Database to do a BLASTp search of these sequences in the Mus Musculus PEP file 

The Ensembl Mus database was dowloaded at 8:45 am EST on 3/9/17:

Index of /pub/release-87/fasta/mus_musculus/pep/

ensembl Mus Musculus pep file downloaded March 9th (10.7 MB, last date modified 11/24/16)

ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz


	mkdir MUSp

	wget ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz

	gzip -d Mus_musculus.GRCm38.pep.all.fa.gz


	makeblastdb -in /mnt/data3/lauren/proteome/MUSp/Mus_musculus.GRCm38.pep.all.fa -out MUSpep -dbtype prot


This is a protein database, and I have a protein query for my proteome, so it is a BLASTp search:

	blastp -query /mnt/data3/lauren/proteome/proteome.fasta \
	-db /mnt/data3/lauren/proteome/MUSp/MUSpep \
	-max_target_seqs 1 \
	-outfmt '6 qseqid pident evalue stitle' \
	-evalue 1e-5 -num_threads 10 | tee proteomeHITS.txt


***This output file (proteomeHITS.txt) is the P. eremicus proteome contig list annotated with Mus Musculus protein accession numbers, I have downloaded it, and it is called SVproteomeBLAST.fasta***


--

**4) Next I did a gene ontology (GO) analysis in Panther. In order to do so I first retrieved the Mus musculus IDs corresponding to the annotated SV proteome sequences**

I have to modify the BLASTp results to make a matrix of all of the contigIDs corresponding to the Mus accession numbers, I will then download this matrix, and I just use all of the Accession numbers to do my panther search

	cat proteomeHITS.txt | wc -l 

1488 (this is the total number of sequencs in the annotated proteome)


	cat proteomeHITS.txt | awk '{print $1 "\t" $4}' > NEWproteomeHITS.txt

	cat NEWproteomeHITS.txt | sort -uk1,1 > SORTproteomeHITS.txt	

	cat SORTproteomeHITS.txt | wc -l


1141 (total number of lines in file)

SORTproteomeHITS.txt is file that only has one contigID per MUS Ensembl gene ID that has matched, but it has duplicate gene matches; 

Now we must remove duplicate gene matches:

	cat SORTproteomeHITS.txt | sort -uk2,2 > newSORTproteomeHITS.txt 

	cat newSORTproteomeHITS.txt | wc -l

1084 (total number of lines in file)


I am proceeding with the file newSORTproteomeHITS.txt, which I download.  

The first column is the Contig ID matches, and the second column is the contig ID matches.

***I retained only the second column, which is Mus Accession IDs, and this is called AccessionONLY.csv.  This file is available on GitHub***

***I uploaded AccessionONLY.csv to Panther (pantherdb.org) and conducted my search against Mus musculus ***

***I have uploaded all of the PANTHER results to GitHub***


***This is the R generated graph with the Biological Processes GO data from PANTHER:***


	data<-data.frame("Biological_Process"=c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                 "Genes"=c(396,339,115,100,68,53,43,43,25,11,5,2,1))
	barplot(data$Genes,names.arg=data$Biological_Process,col=rainbow(13),
        ylim=c(0,400),xlim=c(0,28),xlab="Biological Process",ylab="Genes")

	leg.txt<-c('1: metabolic process',     
           '2: cellular process', '3: cellular component organization', '4: 	localization', '5: biological regulation',
           '6: response to stimulus', '7: developmental process', '8: multicellular organismal process',
           '9: immune system process', '10: biological adhesion', '11: reproduction', '12: locomotion',
           '13: growth')
	legend(16,360,leg.txt, cex=0.75, text.width=11.4, y.intersp=1.5)

 