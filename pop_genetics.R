library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(popgenr)


data(snp)

plot(snp$type)

curve(x^2, 0, 1, xlab="Allele frequencies", 
      ylab="Genotype frequencies", col="green", lwd=2)
text(0.6, 0.2, "Homozygotes", col="green")
curve(2*x*(1-x), 0, 1, add=TRUE, xlab="Allele frequencies", 
      ylab="Genotype frequencies", col="blue", lwd=2)
text(0.25, 0.5, "Heterozygotes", col = "blue")
points(snp$p, snp$hom, pch = 19, col = "green")
points(snp$p, snp$het, pch = 19, col = "blue")


AA <- 34
AS <- 100
SS <- 0

n <- AA+AS+SS
p <- (SS + (AS/2))/n

EAA <- n*(1-p)^2
EAS <- n*2*p*(1-p)
ESS <- n*p^2

geno <- c(AA,AS,SS)
expe <- c(EAA,EAS,ESS)

chi2 <- sum((expe-geno)^2/expe)
print(paste("chi-square:",chi2))

pvalue <- pchisq(chi2, df=1, lower.tail=FALSE)
print(paste("P-value:",pvalue))

dat <- matrix(c(geno,expe), nrow = 2, byrow = T)
barplot(dat,beside=T, col=c("turquoise4", "sienna1"), names.arg=c("AA", "SA", "SS"))
legend(x="topright", legend=c("Observed","Expected"), pch=15, col=c("turquoise4","sienna1"))


#SELECTION: 1st part without and 2nd - under selection

#1 part
p <- 0.25 #Initial allele frequency
gen <- 100 #Number of generations
N <- 1000 #Population size

#Initialize a plot that we can add lines to later

plot(x=NULL, y=NULL, xlim=c(1,gen), ylim=c(0,1), xlab='Generations', 
     ylab='Allele frequency')

for(j in 1:(gen-1)){     
  #Draw the number of alleles 
  a <- rbinom(n=1,size=2*N,prob=p[j]) 
  f <- a/(2*N)  #Get the allele frequency 
  p <- c(p,f) #Save the frequency 
  }
  
lines(x=1:gen, y=p, lwd=2)


#2 part
p <- 0.25 #Initial allele frequency
gen <- 100 #Number of generations
N <- 1000 #Population size
s <- 0.1 #Selection coefficient

#Initialize a plot that we can add lines to later

plot(x=NULL, y=NULL, xlim=c(1,gen), ylim=c(0,1), xlab='Generations', 
     ylab='Allele frequency')

for(j in 1:(gen-1)){     
  #Draw the number of alleles 
  a <- rbinom(n=1,size=2*N,prob=p[j]) 
  f <- a/(2*N)  #Get the allele frequency 
  p <- c(p,(f*(1+s))/(f*s+1))
}

lines(x=1:gen, y=p, lwd=2)


# biomaRt usage and getting SNPs for psoriasis

library(biomaRt)
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snps <- getBM(attributes = c('refsnp_id','allele','phenotype_description','allele_1', 'associated_gene'), 
              filters = c('phenotype_description'), 
              values = "psoriasis", 
              mart = snpmart)

attrib <- listAttributes(snpmart)
filters <- listFilters(snpmart) 



#IRanges and Genomic Ranges

ir1 <- IRanges(start = c(1, 3, 7), width = 2)
ir2 <- IRanges(start = c(3, 5, 6), width = 2)
intersect(ir1, ir2)

ir1 <- IRanges(start = c(1,5,8), width = 5)
names(ir1) <- paste("A", 1:3, sep="")
ir2 <- IRanges(start = c(1,3,5), width = 3)
ir <- c(ir1,ir2)
plotRanges <- function(x, xlim=x, main = deparse(substitute(x)),
                       col ="black", sep =0.5, ...){
  height =1
  if(is(xlim, "Ranges"))
    xlim <-c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x)+1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins*(height + sep) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom+height, col=col,...)
  title(main)
  axis(1)
}

par(mfrows = c(2,1))


plotRanges(ir)
plotRanges(reduce(ir))
plotRanges(disjoin(ir))
plotRanges(ir1)
plotRanges(ir2)
plotRanges(intersect(ir1, ir2))


resize(ir, width = 1, fix = "center"  )

union(ir1,ir2)

reduce(c(ir1,ir2)) #Error in getListElement(x, i, ...) : IRanges objects don't support [[, as.list(), lapply(), or unlist() at the moment

ov <- findOverlaps(ir1,ir2)
ov
queryHits(ov)
unique(queryHits(ov))

# GRanges usage


gr = GRanges(seqnames = c("chr1"), strand = c("-", "+", "+"), ranges = IRanges(start = c(1,3,5), 
                                                                               width = 3))

flank(gr,4) #depends on the directionality of the strand
promoters(gr) #it uses the standart 2Kb length

seqinfo(gr) #not much we know about our sequence now, but..
seqlengths(gr) <- c("chr1"=10)
seqinfo(gr)
seqlevels(gr) # the answer [1] "chr1"

gaps(gr) #gives the all stuff on chromosome not covered, surprisingly 
        # * stand (means or we dont know the direction or it is for both directions) is entirely covered

seqnames(gr) <- c("chr1", "chr2") # we will get the error because oruginally we told that we have only one
seqlevels(gr) <- c("chr1", "chr2") # we tell that values can be assigned to those
# and now we name the ranges

seqnames(gr) <- c("chr1", "chr2", "chr1")
sort(gr)

seqlevels(gr) <- c("chr2", "chr1") #when we change order of levels, when sorting command works different way
sort(gr)

genome(gr) = "hg19"
seqinfo(gr)
gr2 = gr
genome(gr2)="hg18"
findOverlaps(gr,gr2) #will give an error because they are from different genomes


# practicing Genomic Ranges

df = DataFrame(ir = ir1, score = rnorm(3)) #DataFrame in Genomic Ranges is a specific and very useful formal
df[1,1]
df1 = data.frame(ir1 = ir1, score = rnorm(3))
df1[1,1]

values(gr) <-DataFrame(score = rnorm(3))

gr2 = GRanges(seqnames = c("chr1", "chr1", "chr2"), strand = c("-", "+", "+"), ranges = IRanges(start = c(1,3,5), width=2))
                                                                               
findOverlaps(gr,gr2) # with * strand we get overlap in "-" and "+" but if we want we can ignore that
findOverlaps(gr,gr2, ignore.strand = TRUE)

#we can choose a set by overlaps 

subsetByOverlaps(gr2, gr)

# we can convert classic data.frames into GRanges DataFrames

df2 = data.frame(chr = "chr1", start = 1:3, end = 4:6, score = rnorm(3))

makeGRangesFromDataFrame(df2) #"score column drops but we can keep it by...

makeGRangesFromDataFrame(df2, keep.extra.columns = TRUE)

sessionInfo()



### Biology use-cases

#Biology usecase I
#Suppose we want to identify transcription factor (TF) binding sites that overlaps known SNPs.

#Input objects are
#snps: a GRanges (of width 1)
#TF: a GRanges
  
# findOverlaps(snps, TF)

#!!! (watch out for strand)


# CASE II

#Suppose we have a set of differentially methylation regions (DMRs) (think genomic regions) and a set of CpG Islands and we want to find all DMRs within 10kb of a CpG Island.

#Input objects are
#dmrs: a GRanges
#islands: a GRanges


#big_islands <- resize(islands, width = 20000 + width(islands), fix = "center")
#findOverlaps(dmrs, big_islands)

#!!! (watch out for strand)

library(GenomeInfoDb)


gr <-GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))

seqlevels(gr, pruning.mode="coarse") <- "chr2"
gr

dropSeqlevels(gr, "chr2", pruning.mode="coarse")

# there is as well a way to dropchromosomes withw eird names (usually short contigs)

gr <-GRanges(seqnames = c("chr1", "chr23888f"), ranges = IRanges(start = 1:2, end = 4:5))

keepStandardChromosomes(gr, pruning.mode="coarse")
 
# different resources can use different names for chromosomes like chr01, chr1, I...
# we can convert them to one representation

gr <-GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))

newstyle = mapSeqlevels(seqlevels(gr), "NCBI")

newstyle

gr = renameSeqlevels(gr, newstyle)
gr


### AnnotationHub


ah = AnnotationHub()

ah[[1]] # it would retrieve the data from the net
ah = subset(ah, species == "Homo sapiens")

ah2 = display(ah)


#in GRanges for that we can look at signal value - the higher signal the more enrichment

### WE ARE INTERESTED IN METHYLATON AT PROMOTERS 

summary(width(gr1)) # we look at how big those peaks are

summary(width(gr2))
table(width(gr2)) #shows table of frequencies

peaks = gr2

qhs = query(ah, "RefSeq Genes") # there are 4 queries but we are interesred only in hg19 genome

genes = qhs[[1]] # <IRangesList> gives the info about exons
table(genes$name) #frequencies of names
table(table(genes$name)) # frequencies of frequencies

# above we can see how often we see single transctipt and so on

prom = promoters(genes)
table(width(prom)) # all promoters are exactly 2200 bp

args(promoters)

ov = findOverlaps(prom, peaks) #how many methylation marks overlap with promoters

length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
subsetByOverlaps(peaks, prom, ignore.strand=TRUE)
length(subsetByOverlaps(peaks, prom, ignore.strand=TRUE))

length(subsetByOverlaps(peaks, prom, ignore.strand=TRUE)) / length(peaks)
length(subsetByOverlaps(prom, peaks, ignore.strand=TRUE)) / length(prom)
# there are two different questions

sum(width(reduce(peaks, ignore.strand=TRUE))) / 10^6
sum(width(reduce(prom, ignore.strand=TRUE))) / 10^6

sum(width(intersect(peaks, prom,ignore.strand=TRUE ))) / 10^6


# how concluding significance of results we need to build a matrix

inOut = matrix(0, ncol=2, nrow = 2)
colnames(inOut) = c("In", "Out")
rownames(inOut) = c("In", "Out")


inOut[1,1] = sum(width(intersect(peaks, prom, ignore.strand=TRUE )))
inOut[1,2] = sum(width(setdiff(peaks, prom, ignore.strand=TRUE)))
inOut[2,1] = sum(width(setdiff(prom, peaks, ignore.strand=TRUE)))
inOut
colSums(inOut)
rowSums(inOut)
inOut[2,2] = 3*10^9 - sum(inOut)

# after we have the matrix we can calculate odds ratio

# normally fisher.test(inOut)$statistic  would work but in genome case the integers in R don't have such number
# we will calculate manually 

oddsRatio = inOut[1,1]*inOut[2,2] / (inOut[2,1]*inOut[1,2])
oddsRatio



### EXCERSICES

## Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
# Question 1: How many islands exists on the autosomes?

ah = AnnotationHub()

ah[[1]] # it would retrieve the data from the net
ah = subset(ah, species == "Homo sapiens")

cg = query(ah,"CpG Islands")
cg$dataprovider
cg$title
CG = cg[[1]]
sum(width(CG)) #21842742
CG <- keepStandardChromosomes(CG, pruning.mode="coarse")
sum(width(CG)) #21842742
auto <- extractSeqlevelsByGroup(species="Homo sapiens", style = "UCSC",group="auto")
CG_auto <- keepSeqlevels(CG,auto, pruning.mode="coarse")
length(CG_auto) #26641 right


# Question 2: How many CpG Islands exists on chromosome 4.

CG_on_4th <- CG[seqnames(CG)=="chr4"]
length(CG_on_4th) #1031


# Obtain the data for the H3K4me3 histone modification for the H1 cell line from Epigenomics Roadmap, 
# using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes 
# (chromosomes 1 to 22). Question 3: How many bases does these regions cover?

h3k4H1 = query(ah, c("H3K4me3", "H1")) #the way we specify the search , or we can put all "keywords" in a vector
h3k4H1$dataprovider
h3k4H1$title
H3K4H1 = h3k4H1[[7]]
H3K4H1 <- keepStandardChromosomes(H3K4H1, pruning.mode="coarse")

auto <- extractSeqlevelsByGroup(species="Homo sapiens", style = "UCSC",group="auto")
H3K4H1_auto <- keepSeqlevels(H3K4H1,auto, pruning.mode="coarse")

sum(width(H3K4H1_auto)) #41135164 



# Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, using the AnnotationHub package. 
# Subset these regions to only keep regions mapped to the autosomes. In the return data, each region has an associated "signalValue". 
# Question 4: What is the mean signalValue across all regions on the standard chromosomes?

h3k27H1 = query(ah, c("H3K27me3", "H1"))
h3k27H1$dataprovider
h3k27H1$title
H3K27H1 = h3k27H1[[7]]
H3K27H1 <- keepStandardChromosomes(H3K27H1, pruning.mode="coarse")
auto <- extractSeqlevelsByGroup(species="Homo sapiens", style = "UCSC",group="auto")
H3K27H1_auto <- keepSeqlevels(H3K27H1, auto, pruning.mode="coarse")

summary(H3K27H1_auto$signalValue) #4.771


# Bivalent regions are bound by both H3K4me3 and H3K27me3.
# Question 5: Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?

bivalent_auto <- intersect(H3K4H1_auto,H3K27H1_auto)
sum(width(bivalent_auto)) # 10289096

bivalent <- intersect(H3K4H1,H3K27H1)
sum(width(bivalent))

### Question 6: how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands? 

length(subsetByOverlaps(bivalent_auto, CG_auto)) / length(bivalent_auto) # 0.5383644


### Question 7: How big a fraction (expressed as a number between 0 and 1) of the bases which are part of CpG Islands, are also bivalent marked

sum(width(intersect(bivalent_auto, CG_auto))) / sum(width(CG_auto)) # 0.241688


### Question 8: How many bases are bivalently marked within 10kb of CpG Islands? Tip: consider using the "resize()"" function.

CG_auto_10k <-CG_auto + 10000
sum(width(subsetByOverlaps(bivalent_auto, CG_auto_10k, ignore.strand = TRUE))) #9785641 it's not right answer as well (closest answer is right)

### Question 9: How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?

sum(width(CG)) / (3*10^9) # 0.007050796


### Question 10: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.

inOut = matrix(0, ncol=2, nrow = 2)
colnames(inOut) = c("In", "Out")
rownames(inOut) = c("In", "Out")

inOut[1,1] = sum(width(intersect(bivalent, CG, ignore.strand=TRUE )))
inOut[1,2] = sum(width(setdiff(bivalent, CG, ignore.strand=TRUE)))
inOut[2,1] = sum(width(setdiff(CG,bivalent, ignore.strand=TRUE)))
inOut
colSums(inOut)
rowSums(inOut)
inOut[2,2] = 3*10^9 - sum(inOut)

oddsRatio = inOut[1,1]*inOut[2,2] / (inOut[2,1]*inOut[1,2])
oddsRatio # 168.0189


Oragnism's supported by GenomeInfoDb can be found by :

    names(genomeStyles())
  
    
    If one knows the organism one is interested in, then we can directly access 
  the information for the given organism along. Each function accepts an 
  argument called species which as "genus species", the default 
  is "Homo sapiens". In the following example we list out only the first five 
  entries returned by the code snippet.
  
  
  <<genomeStyles2>>=
    head(genomeStyles("Homo_sapiens"),5)
  @ 
    %% 
    
    We can also check if a given style is supported by GenomeInfoDb for a given
  species. For example, if we want to know if "UCSC" mapping is supported for 
  "Homo sapiens" we can ask :
    
    <<style-present>>=
    "UCSC" %in% names(genomeStyles("Homo_sapiens"))
  @
    
    
    \subsection{extractSeqlevels}
  We can also extract the desired seqlevelsStyle from a given organism using
  the \Rfunction{extractSeqlevels}
  <<extractSeqlevels>>=
    extractSeqlevels(species="Arabidopsis_thaliana", style="NCBI")
  @ 
    %%
    
    \subsection{extractSeqlevelsByGroup}
  We can also extract the desired seqlevelsStyle from a given organism based on 
  a group ( Group - 'auto' denotes autosomes, 'circular' denotes circular
            chromosomes and 'sex' denotes sex chromosomes; the default is all chromosomes
            are returned).
  <<extractSeqlevelsgroup>>=
    extractSeqlevelsByGroup(species="Arabidopsis_thaliana", style="NCBI",
                            group="auto")
  @ 
    %%
    
    \subsection{seqlevelsStyle}
  We can find the seqname Style for a given character vector by using the 
  \Rfunction{seqlevelsStyle}
  <<seqlevelsStyle>>=
    seqlevelsStyle(paste0("chr",c(1:30)))
  seqlevelsStyle(c("2L","2R","X","Xhet"))
  @
    %%
    
    \subsection{seqlevelsInGroup}
  We can also subset a given character vector containing seqnames  using the
  \Rfunction{seqlevelsInGroup}.
  We currently support 3 groups: 'auto' for autosomes, 'sex' for allosomes/sex 
  chromosomes and circular for 'circular' chromosomes.
  The user can also provide the style and species they are working with.
  In the following examples, we extract the sex, auto and circular chromosomes 
  for Homo sapiens :
    
    <<keepChr-txdb>>=
    newchr <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
  seqlevelsInGroup(newchr, group="sex")
  seqlevelsInGroup(newchr, group="auto")
  seqlevelsInGroup(newchr, group="circular")
  seqlevelsInGroup(newchr, group="sex","Homo_sapiens","UCSC")
  @
    %%
    
    if we have a vector containing seqnames and we want to verify the
  species and style for them , we can use:
    <<check2>>=
    seqnames <- c("chr1", "chr9", "chr2", "chr3", "chr10")
  all(seqnames %in% extractSeqlevels("Homo_sapiens", "UCSC"))
  @
    
    \subsection{orderSeqlevels}
  The \Rfunction{orderSeqlevels} can return the order of a given character vector 
  which contains seqnames.In the following example, we show how you can find the 
  order for a given seqnames character vector.
  <<orderSeqlevels>>=
    seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
  orderSeqlevels(seqnames)
  seqnames[orderSeqlevels(seqnames)]
  @
    %%
    
    \subsection{rankSeqlevels}
  The \Rfunction{rankSeqlevels} can return the rank of a given character vector 
  which contains seqnames.In the following example, we show how you can find the 
  rank for a given seqnames character vector.
  <<rankSeqlevels>>=
    seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
  rankSeqlevels(seqnames)
  @
    %%
    
    \subsection{mapSeqlevels}
  Returns a matrix with 1 column per supplied sequence name and 1 row per 
  sequence renaming map compatible with the specified style.  If \Rcode{best.only}
  is \Rcode{TRUE} (the default), only the "best" renaming maps (i.e. the rows with
                                                                less NAs) are returned.
  <<find>>=
    mapSeqlevels(c("chrII", "chrIII", "chrM"), "NCBI")
  @ 
    %%
    
    We also have several seqlevel utility functions.Let us construct a basic 
  GRanges and show how these functions can be used. .
  
  <<basic-gr>>=
    gr <- GRanges(paste0("ch",1:35), IRanges(1:35, width=5))
  gr
  @
    
    As you can see , we have "ch" instead of "chr" for chromosome names. 
  We can use \Rfunction{renameSeqlevels} to change the "ch" to "chr"
  
  \subsection{renameSeqlevels}
  As the first argument - it takes the object whose seqlevels we need to change, 
  and as the second argument it takes a named vector which has the changes.
  
  <<renameseqlevels>>=
    newnames <- paste0("chr",1:35)
  names(newnames) <- paste0("ch",1:35)
  head(newnames)
  gr <- renameSeqlevels(gr,newnames)
  gr
  @
    
    Humans have just 22 primary chromosomes - but here we have some extra seqlevels
  which we want to remove - there are several ways we can achieve this:
    
    \subsection{dropSeqlevels}  
  Here the second argument is the seqlevels that you want to drop. Because
  these seqlevels are in use (i.e. have ranges on them), the ranges on these
  sequences need to be removed before the seqlevels can be dropped. We call
  this {\it pruning}. The \Rcode{pruning.mode} argument controls how to prune
  \Rcode{gr}. Unlike for list-like objects (e.g. GRangesList) for which pruning
  can be done in various ways, pruning a GRanges object is straightforward and
  achieved by specifying \Rcode{pruning.mode="coarse"}.
  <<dropseqlevels>>=
    dropSeqlevels(gr, paste0("chr",23:35), pruning.mode="coarse")
  @
    
    \subsection{keepSeqlevels}  
  Here the second argument is the seqlevels that you want to keep.
  <<keepseqlevels>>=
    keepSeqlevels(gr, paste0("chr",1:22), pruning.mode="coarse")
  @
    
    \subsection{keepStandardChromosomes}  
  This function internally uses the pre-defined tables inside GenomeInfoDb to 
  find the correct seqlevels according to the sequence style of the object. 
  <<keepstdchr>>=
    keepStandardChromosomes(gr, pruning.mode="coarse")
  @
    One can also specify the optional species argument to bemore precise. 
  
  <<keepstdchr-2>>=
    plantgr <- GRanges(c(1:5,"MT","Pltd"), IRanges(1:7,width=5))
  keepStandardChromosomes(plantgr, species="Arabidopsis thaliana",
                          pruning.mode="coarse")
  @
    
    \section{Classes inside GenomeInfoDb package}
  
  \subsection{Genome-Description class}
  
  We also provide a Genome Description class which can be used in the following
  way:
    
    <<genome-description-class, message=FALSE>>=
    library(BSgenome.Celegans.UCSC.ce2)
  class(Celegans)
  is(Celegans, "GenomeDescription")
  provider(Celegans)
  seqinfo(Celegans)
  gendesc <- as(Celegans, "GenomeDescription")
  class(gendesc)
  gendesc
  provider(gendesc)
  seqinfo(gendesc)
  bsgenomeName(gendesc)
  @
    
    \subsection{Seqinfo class}
  
  <<Seqinfo-egs>>=
    ## Note that all the arguments (except 'genome') must have the
    ## same length. 'genome' can be of length 1, whatever the lengths
    ## of the other arguments are.
    x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
                 seqlengths=c(100, 200, NA, 15),
                 isCircular=c(NA, FALSE, FALSE, TRUE),
                 genome="toy")
  length(x)
  seqnames(x)
  names(x)
  seqlevels(x)
  seqlengths(x)
  isCircular(x)
  genome(x)
  
  x[c("chrY", "chr3", "chr1")]  # subset by names
  
  ## Rename, drop, add and/or reorder the sequence levels:
  xx <- x
  seqlevels(xx) <- sub("chr", "ch", seqlevels(xx))  # rename
  xx
  seqlevels(xx) <- rev(seqlevels(xx))  # reorder
  xx
  seqlevels(xx) <- c("ch1", "ch2", "chY")  # drop/add/reorder
  xx
  seqlevels(xx) <- c(chY="Y", ch1="1", "22")  # rename/reorder/drop/add
  xx
  
  y <- Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
               seqlengths=c(300, NA, 15))
  y
  merge(x, y)  # rows for chr3 and chrM are merged
  suppressWarnings(merge(x, y))
  
  ## Note that, strictly speaking, merging 2 Seqinfo objects is not
  ## a commutative operation, i.e., in general 'z1 <- merge(x, y)'
  ## is not identical to 'z2 <- merge(y, x)'. However 'z1' and 'z2'
  ## are guaranteed to contain the same information (i.e. the same
  ## rows, but typically not in the same order):
  suppressWarnings(merge(y, x))
  
  ## This contradicts what 'x' says about circularity of chr3 and chrM:
  isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
  y
  if (interactive()) {
    merge(x, y)  # raises an error
  }
  @
    
    \section{Examples}
  
  \subsection{converting seqlevel styles (eg:UCSC to NCBI)}
  
  A quick example using Drosophila Melanogaster. The txdb object contains 
  seqlevels in UCSC style, we want to convert them to NCBI 
  
  <<quick-style>>=
    txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
  seqlevels(txdb)
  genomeStyles("Drosophila melanogaster")
  mapSeqlevels(seqlevels(txdb), "NCBI")
  @
    
    \subsection{converting styles and removing unwanted seqlevels}
  
  Suppose we read in a Bam file or a BED file and the resulting GRanges have a lot
  of seqlevels which are not required by your analysis or you want to rename
  the seqlevels from the current style to your own style (eg:USCS to NCBI), we can 
  use the functionality provided by GenomeInfoDb to do that.
  
  Let us say that we have extracted the seqlevels of the Seqinfo object(say 
                                                                        GRanges from a BED file) in a variable called "sequence". 
  <<sequence, eval=FALSE>>=
    sequence <- seqlevels(x) 
  
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,"NCBI")
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  
  ## rename the seqlevels 
  x <- renameSeqlevels(x,newStyle)
  
  ## keep only the seqlevels you want (say autosomes)
  auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI", 
                                  group="auto")
  x <- keepSeqlevels(x,auto)
  @
    
    \section{Session Information}
  Here is the output of \Rfunction{sessionInfo} on the system on which
  this document was compiled:
    <<sessionInfo, results='asis', eval=TRUE>>=
    toLatex(sessionInfo())
  @
    
    \end{document}













