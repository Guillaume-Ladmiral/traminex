# TraMineR Manual — Collected Example Code

The following R code blocks were extracted from the **TraMineR** package manual's *Examples* sections.

---

## theTraMineR

```r
## load the mvad data
library(TraMineR)
data(mvad)
## create a state sequence object from columns 17 to 86
mvad.seq <- seqdef(mvad[,17:86])
## distribution plot by sex (male)
seqdplot(mvad.seq, group=mvad$male, border=NA)
## compute the LCS pairwise distance matrix
## among the first 10 sequences
mvad.lcs <- seqdist(mvad.seq[1:10,], method="LCS")
actcal Example data set: Activity calendar from the Swiss Household Panel
```

---

## Swiss

```r
data(actcal)
actcal.seq <- seqdef(actcal[,13:24])
## Defining the transition matrix
transition <- seqetm(actcal.seq, method="transition")
transition[1,1:4] <- c("FullTime" , "Decrease,PartTime",
"Decrease,LowPartTime", "Stop")
transition[2,1:4] <- c("Increase,FullTime", "PartTime" ,
"Decrease,LowPartTime", "Stop")
transition[3,1:4] <- c("Increase,FullTime", "Increase,PartTime",
"LowPartTime" , "Stop")
transition[4,1:4] <- c("Start,FullTime" , "Start,PartTime" ,
"Start,LowPartTime" , "NoActivity")
## Converting STS data to TSE
actcal.tse <- seqformat(actcal, 13:24, from = "STS",to = "TSE",
tevent = transition)
## Defining the event sequence object
actcal.eseq <- seqecreate(id=actcal.tse$id,
time=actcal.tse$time, event=actcal.tse$event)
```

---

## Alexis

```r
## Creating a sequence object with the columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Retrieving the alphabet
alphabet(actcal.seq)
## Setting the alphabet
alphabet(actcal.seq) <- c("FT", "PT", "LT", "NO")
## Event sequences
actcal.eseq <- seqecreate(actcal.seq)
alphabet(actcal.eseq)
bfspell Example data set: First 20 biofam sequences in SPELL form
```

---

## Alexis

```r
## Creating a sequence object with the columns 13 to 24
## in the 'actcal 'example data set
## The color palette is automatically set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Retrieving the color palette
cpal(actcal.seq)
seqiplot(actcal.seq)
## Setting a user defined color palette
cpal(actcal.seq) <- c("blue","red", "green", "yellow")
seqiplot(actcal.seq)
dissassoc Analysis of discrepancy from dissimilarity measures
```

---

## dissmfacw

```r
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])
## Building dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")
## R=1 implies no permutation test
da <- dissassoc(mvad.ham, group=mvad$gcse5eq, R=10)
print(da)
hist(da)
disscenter Compute distances to the center of a group
```

---

## dissmfacw

```r
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])
## Computing dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")
## Distance to center according to group gcse5eq
dc <- disscenter(mvad.ham, group=mvad$gcse5eq)
## Ploting distribution of dissimilarity to center
boxplot(dc~mvad$gcse5eq, col="cyan")
## Retrieving index of the first medoids, one per group
dc <- disscenter(mvad.ham, group=mvad$Grammar, medoids.index="first")
print(dc)
## Retrieving index of all medoids in each group
dc <- disscenter(mvad.ham, group=mvad$Grammar, medoids.index="all")
print(dc)
## Distances of all cases to each group center
dc <- disscenter(mvad.ham, group=mvad$Grammar, allcenter=TRUE)
print(cbind(dc,mvad$Grammar)[1:10,])
## Compare with the following
dc0 <- disscenter(mvad.ham, group=mvad$Grammar)
print(cbind(dc0,mvad$Grammar)[1:10,])
dissdomassoc Domain association measures
```

---

## Research

```r
data(biofam)
## Building one channel per type of event left, children or married
## Using only first 200 sequences
bf <- as.matrix(biofam[1:200, 10:25])
children <- bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6
## weights
weights <- biofam[1:200,"wp00tbgs"]
## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)
## distances by channel
dchild <- seqdist(child.seq, method="OM", sm="INDELSLOG")
dmarr <- seqdist(marr.seq, method="OM", sm="INDELSLOG")
dleft <- seqdist(left.seq, method="OM", sm="INDELSLOG")
dbiofam <- list(dchild,dmarr,dleft)
names(dbiofam) <- c("child","marr","left")
## Association between domains and R2 by domain
rass <- dissdomassoc(dbiofam, weights=weights)
## Joint distances using additive trick with domain INDELSLOG costs.
mcdist <- seqMD(channels=list(child.seq, marr.seq, left.seq), what="diss",
method="OM", sm =list("INDELSLOG", "INDELSLOG", "INDELSLOG"))
rassj <- dissdomassoc(dbiofam, jointdiss=mcdist, what=c("all"), weights=weights)
summary(rassj)
dissmergegroups Merging groups by minimizing loss of partition quality.
```

---

## social

```r
data(biofam)
## Building one channel per type of event (children, married, left home)
cases <- 1:40
bf <- as.matrix(biofam[cases, 10:25])
children <- bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6
## Creating sequence objects
child.seq <- seqdef(children, weights = biofam[cases, 'wp00tbgs '])
marr.seq <- seqdef(married, weights = biofam[cases, 'wp00tbgs '])
left.seq <- seqdef(left, weights = biofam[cases, 'wp00tbgs '])
## distances by domain
dchild <- seqdist(child.seq, method="OM", sm="INDELSLOG")
dmarr <- seqdist(marr.seq, method="OM", sm="INDELSLOG")
dleft <- seqdist(left.seq, method="OM", sm="INDELSLOG")
dnames <- c("child","marr","left")
## clustering each domain into 2 groups
child.cl2 <- cutree(hclust(as.dist(dchild)),k=2)
marr.cl2 <- cutree(hclust(as.dist(dmarr)),k=2)
left.cl2 <- cutree(hclust(as.dist(dleft)),k=2)
## Multidomain sequences
MD.seq <- seqMD(list(child.seq,marr.seq,left.seq))
d.expand <- seqdist(MD.seq, method="LCS")
clust.comb <- interaction(child.cl2,marr.cl2,left.cl2)
merged.grp <- dissmergegroups(d.expand, clust.comb,
weights=biofam[cases, 'wp00tbgs '])
## weighted size of merged groups
xtabs(biofam[cases, 'wp00tbgs '] ~ merged.grp)
dissmfacw Multi-factor ANOVA from a dissimilarity matrix
dissmfacw 21
```

---

## disscenter

```r
## Define the state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])
## Here, we use only first 100 sequences
mvad.seq <- mvad.seq[1:100,]
dissrep 23
## Compute dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")
## And now the multi-factor analysis
print(dissmfacw(mvad.ham ~ male + Grammar + funemp +
gcse5eq + fmpr + livboth, data=mvad[1:100,], R=10))
dissrep Extracting sets of representative objects using a dissimilarity matrix
```

---

## seqrep

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam[, 10:25], labels=biofam.lab)
## Computing the distance matrix
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", sm=costs)
## Representative set using the neighborhood density criterion
biofam.rep <- dissrep(biofam.om)
summary(biofam.rep)
## indexes of first occurrence of second representatives in original data
attr(biofam.rep,"idx.rep")[[2]][1]
dissrf Relative Frequency Groups.
```

---

## plot.seqrf

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
## Here, we use only 100 cases selected such that all elements
## of the alphabet be present.
## (More cases and a larger k would be necessary to get a meaningful example.)
biofam.seq <- seqdef(biofam[501:600, 10:25], labels=biofam.lab,
weights=biofam[501:600,"wp00tbgs"])
diss <- seqdist(biofam.seq, method="LCS")
## Using 12 groups, default MDS sorting,
## and original method by Fasang and Liao (2014)
dissrf(diss=diss, k=12, grp.meth="first")
## Using 12 groups, weights, default MDS sorting,
## and default "prop" method
w <- attr(biofam.seq, "weights")
dissrf(diss=diss, k=12, weights=w)
## With a user specified sorting variable
## Here time spent in parental home, which has ties
parentTime <- seqistatd(biofam.seq)[, 1]
b.srf <- seqrf(biofam.seq, diss=diss, k=12, sortv=parentTime)
## print, summary, and plot methods
summary(b.srf)
plot(b.srf)
plot(b.srf, which.plot="both")
```

---

## disscenter

```r
data(mvad)
## Defining a state sequence object
mvad.seq <- seqdef(mvad[, 17:86])
## Computing dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")
## Grow the tree using a low R value for illustration.
## For R=10, pval cannot be lower than 0.1
dt <- disstree(mvad.ham~ male + Grammar + funemp + gcse5eq + fmpr + livboth,
data=mvad, R = 10, pval = 0.1)
print(dt)
## Will only work if GraphViz is properly installed
## See seqtree for a simpler way to plot a sequence tree.
## Not run:
disstreedisplay(dt, image.fun = seqdplot, image.data = mvad.seq,
## Additional parameters passed to seqdplot
with.legend = FALSE, xaxis = FALSE, ylab = "", border=NA)
## End(Not run)
## Second method, using a specific function
myplotfunction <- function(individuals, seqs, ...) {
par(font.sub=2, mar=c(3,0,6,0), mgp=c(0,0,0))
## using mds to order sequence in seqIplot
mds <- suppressMessages(cmdscale(seqdist(seqs[individuals,], method="HAM"),k=1))
seqIplot(seqs[individuals,], sortv=mds,...)
}
## If image.data is not set, indexes of individuals are sent to image.fun
## Not run:
disstreedisplay(dt, image.fun = myplotfunction, cex.main = 3,
## additional parameters passed to myplotfunction
seqs = mvad.seq,
## additional parameters passed to seqIplot (through myplotfunction)
with.legend = FALSE, xaxis = FALSE, ylab = "")
## End(Not run)
## Retrieving terminal node membership
term.leaf <- disstreeleaf(dt)
table(term.leaf)
## Retrieving classification rules
rules <- disstree.get.rules(dt)
## Index of rule (terminal leaf) that applies to a specified profile
## covariates are: male, Grammar, funemp, gcse5eq, fmpr, livboth
profile <- data.frame(male="no", Grammar="yes", funemp="no", gcse5eq="yes", fmpr="no", livboth="no")
rules[disstree.assign(rules, profile=profile)]
```

---

## dissmfacw

```r
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])
## Building dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")
## Pseudo variance of the sequences
print(dissvar(mvad.ham))
ex1 Example data set with missing values and weights
```

---

## The

```r
data(ex2)
ex2w.seq <- seqdef(ex2.weighted, 1, weights=ex2.weighted$weight)
ex2u.seq <- seqdef(ex2.unweighted)
famform Example data set: sequences of family formation
```

---

## Gilbert

```r
## Creating a sequence object with the columns 13 to 24
## in the 'actcal 'example data set
data(biofam)
biofam <- biofam[sample(nrow(biofam),300),]
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam[,10:25], weights=biofam$wp00tbgs)
is.stslist(biofam.seq) #TRUE
attr(biofam.seq,"weights") <- NULL
is.stslist(biofam.seq) #TRUE
attr(biofam.seq,"weights") <- rep(1, nrow(biofam.seq))
is.stslist(biofam.seq) #FALSE
w <- rep(1, nrow(biofam.seq))
names(w) <- rownames(biofam.seq)
attr(biofam.seq,"weights") <- w
is.stslist(biofam.seq) #TRUE
mvad Example data set: Transition from school to work
```

---

## Gilbert

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam <- biofam[500:600,] ## using a subsample only
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)
## Plot of the 10 most frequent sequences
## with bar width proportional to the frequency
plot(biofam.seq)
## Plotting the whole data set
## with no borders
plot(biofam.seq, idxs=0, space=0, border=NA)
## =======
## Weights
## =======
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
plot.stslist.freq 49
plot(ex1.seq)
plot(ex1.seq, weighted=FALSE)
plot.stslist.freq Plot method for sequence frequency tables
```

---

## Alexis

```r
## Loading the 'actcal 'example data set
data(actcal)
## Defining a sequence object with data in columns 13 to 24
## (activity status from january to december 2000)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)
## 10 most frequent sequences in the data
actcal.freq <- seqtab(actcal.seq)
## Plotting the object
plot(actcal.freq, main="Sequence frequencies - actcal data set")
## Plotting all the distinct sequences without borders
## and space between sequences
actcal.freq2 <- seqtab(actcal.seq, idxs=0)
plot(actcal.freq2, main="Sequence frequencies - actcal data set",
border=NA, space=0)
plot.stslist.meant Plot method for objects produced by the seqmeant function
```

---

## automatic

```r
## Loading the mvad data set and creating a sequence object
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
"joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")
mvad.seq <- seqdef(mvad, 15:86, states=mvad.scodes, labels=mvad.labels)
## Computing the mean times
mvad.meant <- seqmeant(mvad.seq)
## Plotting
plot(mvad.meant, main="Mean durations in each state of the alphabet")
## Changing the y axis limits
plot(mvad.meant, main="Mean durations in each state of the alphabet",
ylim=c(0,40))
## Displaying error bars
52 plot.stslist.modst
mvad.meant.e <- seqmeant(mvad.seq, serr=TRUE)
plot(mvad.meant.e, main="Mean durations in each state of the alphabet",
ylim=c(0,40))
plot.stslist.modst Plot method for modal state sequences
```

---

## states

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)
## Modal state sequence
biofam.modst <- seqmodst(biofam.seq)
plot(biofam.modst)
plot.stslist.rep Plot method for representative sequence sets
```

---

## plot.stslist.statd

```r
## Loading the mvad data set and creating a sequence object
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
"joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")
## First 36 months trajectories
mvad.seq <- seqdef(mvad, 15:50, states=mvad.scodes, labels=mvad.labels)
## Computing Hamming distances
##
dist.ham <- seqdist(mvad.seq, method="HAM")
## Extracting a representative set using the sequence frequency
## as a representativeness criterion
mvad.rep <- seqrep(mvad.seq, diss=dist.ham)
## Plotting the representative set
plot(mvad.rep)
plot.stslist.statd Plot method for objects produced by the seqstatd function
```

---

## plot.subseqelist

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)
## State distribution
biofam.statd <- seqstatd(biofam.seq)
## State distribution plot (default type="d" option)
plot(biofam.statd)
## Entropy index plot
plot(biofam.statd, type="Ht")
## State distribution and entropy line
plot(biofam.statd, type="dH")
plot.subseqelist Plot frequencies of subsequences
```

---

## Matthias

```r
## loading data
data(actcal.tse)
## creating sequences
actcal.eseq <- seqecreate(actcal.tse)
## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.eseq,pmin.support=0.01)
## Frequence of first ten subsequences
plot(fsubseq[1:10], cex=2)
plot(fsubseq[1:10])
plot.subseqelistchisq Plot discriminant subsequences
```

---

## seqdef

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam <- biofam[500:600,] ## using a subsample only
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.shortlab <- c("P","L","M","LM","C","LC","LMC","D")
biofam.seq <- seqdef(biofam, 10:25, states=biofam.shortlab,
labels=biofam.lab)
## Print of first 5 sequences
print(biofam.seq[1:5,])
print(biofam.seq[1:5,], extended=TRUE)
print(biofam.seq[1:5,], format="SPS")
print(biofam.seq[1:5,], format="SPS", SPS.out = list(xfix = "", sdsep = "/"))
```

---

## help

```r
data(biofam)
biofam.seq <- seqdef(biofam, 10:25)
costs <- seqsubm(biofam.seq, method="TRATE")
sa <- seqalign(biofam.seq, 1:2, indel=1, sm=costs)
print(sa)
plot(sa)
sa <- seqalign(biofam.seq, c(1,5), indel=0.5, sm=costs)
print(sa)
plot(sa)
```

---

## seqconc

```r
data(mvad)
mvad.shortlab <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, states=mvad.shortlab, 15:86)
## Comparing sequences 1 and 2 in mvad.seq
seqcomp(mvad.seq[1,],mvad.seq[2,])
## Comparing sequences 176 and 211 in mvad.seq
seqcomp(mvad.seq[176,],mvad.seq[211,])
```

---

## seqdecomp

```r
data(actcal)
actcal.string <- seqconc(actcal,13:24)
head(actcal.string)
```

---

## seqtrate

```r
## Defining a sequence object with columns 10 to 25
## of a subset of the 'biofam 'example data set.
data(biofam)
biofam.seq <- seqdef(biofam[501:600,10:25])
## Indel and substitution costs based on log of inverse state frequencies
lifcost <- seqcost(biofam.seq, method="INDELSLOG")
## Here lifcost$indel is a vector
biofam.om <- seqdist(biofam.seq, method="OM", indel=lifcost$indel, sm=lifcost$sm)
## Optimal matching using transition rates based substitution-cost matrix
## and the associated indel cost
## Here trcost$indel is a scalar
trcost <- seqcost(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", indel=trcost$indel, sm=trcost$sm)
## Using costs based on FUTURE with a forward lag of 4
fucost <- seqcost(biofam.seq, method="FUTURE", lag=4)
biofam.om <- seqdist(biofam.seq, method="OM", indel=fucost$indel, sm=fucost$sm)
## Optimal matching using a unique substitution cost of 2
## and an insertion/deletion cost of 3
ccost <- seqsubm(biofam.seq, method="CONSTANT", cval=2)
biofam.om.c2 <- seqdist(biofam.seq, method="OM",indel=3, sm=ccost)
## Displaying the distance matrix for the first 10 sequences
## =================================
## Example with weights and missings
## =================================
data(ex1)
ex1.seq <- seqdef(ex1[,1:13], weights=ex1$weights)
## Unweighted
subm <- seqcost(ex1.seq, method="INDELSLOG", with.missing=TRUE, weighted=FALSE)
ex1.om <- seqdist(ex1.seq, method="OM", indel=subm$indel, sm=subm$sm, with.missing=TRUE)
## Weighted
subm.w <- seqcost(ex1.seq, method="INDELSLOG", with.missing=TRUE, weighted=TRUE)
ex1.omw <- seqdist(ex1.seq, method="OM", indel=subm.w$indel, sm=subm.w$sm, with.missing=TRUE)
ex1.om == ex1.omw
```

---

## seqconc

```r
# 1 sequence of length 4
seqdecomp("A-BB-C-DD")
# 2 sequences of length 6
seqdecomp(c("ABBCDD","BBCCAD"),sep="")
seqdef 69
seqdef Create a state sequence object
```

---

## seqformat

```r
## Creating a sequence object with the columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24,
labels=c("> 37 hours", "19-36 hours", "1-18 hours", "no work"))
## Displaying the first 10 rows of the sequence object
## Displaying the first 10 rows of the sequence object
## in SPS format
print(actcal.seq[1:10,], format="SPS")
## Plotting the first 10 sequences
plot(actcal.seq)
## Re-ordering the alphabet
actcal.seq <- seqdef(actcal,13:24,alphabet=c("B","A","D","C"))
alphabet(actcal.seq)
## Adding a state not appearing in the data to the
## alphabet
actcal.seq <- seqdef(actcal,13:24,alphabet=c("A","B","C","D","E"))
alphabet(actcal.seq)
## Adding a state not appearing in the data to the
## alphabet and changing the states labels
actcal.seq <- seqdef(actcal,13:24,
alphabet=c("A","B","C","D","E"),
states=c("FT","PT","LT","NO","TR"))
alphabet(actcal.seq)
## rbind and summary
seq1 <- actcal.seq[1:10,]
seq2 <- actcal.seq[20:25,]
seq <- rbind(seq1,seq2)
summary(seq)
## ============================
## Example with missing values
## ============================
data(ex1)
## With right="DEL" default value
seqdef(ex1,1:13)
## Eliminating 'left 'missing values
seqdef(ex1,1:13, left="DEL")
## Eliminating 'left 'missing values and gaps
seqdef(ex1,1:13, left="DEL", gaps="DEL")
## ====================
## Example with weights
## ====================
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
## weighted sequence frequencies
seqtab(ex1.seq)
```

---

## dissassoc

```r
## Define a state sequence object
data(mvad)
## First 12 months of first 100 trajectories
mvad.seq <- seqdef(mvad[1:100, 17:28])
## Position-wise discrepancy analysis using
## centered sliding windows of length 5.
mvad.diff <- seqdiff(mvad.seq, group=mvad$gcse5eq[1:100], cmprange=c(-2,2))
print(mvad.diff)
plot(mvad.diff, stat=c("Pseudo R2", "Levene"))
plot(mvad.diff, stat="discrepancy")
```

---

## additive

```r
## =========================
## Examples without missings
## =========================
## Defining a sequence object with columns 10 to 25
## of a subset of the 'biofam 'data set
data(biofam)
biofam.seq <- seqdef(biofam[501:600, 10:25])
## OM distances using the vector of indels and substitution
## costs derived from the estimated state frequencies
costs <- seqcost(biofam.seq, method = "INDELSLOG")
biofam.om <- seqdist(biofam.seq, method = "OM",
indel = costs$indel, sm = costs$sm)
## OM between sequences of transitions
biofam.omstran <- seqdist(biofam.seq, method = "OMstran",
indel = costs$indel, sm = costs$sm,
otto=.3, transindel="subcost")
seqdist 81
## Normalized LCP distances
biofam.lcp.n <- seqdist(biofam.seq, method = "LCP",
norm = "auto")
## Normalized LCS distances to the most frequent sequence
biofam.dref1 <- seqdist(biofam.seq, method = "LCS",
refseq = 0, norm = "auto")
## LCS distances to an external sequence
ref <- seqdef(as.matrix("(0,5)-(3,5)-(4,6)"), informat = "SPS",
alphabet = alphabet(biofam.seq))
biofam.dref2 <- seqdist(biofam.seq, method = "LCS",
refseq = ref)
## LCS distances between two subsets of sequences
set1 <- 1:10
set2 <- 31:36
biofam.dref2 <- seqdist(biofam.seq, method = "LCS",
refseq = list(set1,set2))
## Chi-squared distance over the full observed timeframe
biofam.chi.full <- seqdist(biofam.seq, method = "CHI2",
step = max(seqlength(biofam.seq)))
## Chi-squared distance over successive overlapping
## intervals of length 4
biofam.chi.ostep <- seqdist(biofam.seq, method = "CHI2",
step = 4, overlap = TRUE)
## ======================
## Examples with missings
## ======================
data(ex1)
## Ignore empty row 7
ex1.seq <- seqdef(ex1[1:6, 1:13])
## OM with indel and substitution costs based on
## log of inverse state frequencies
costs.ex1 <- seqcost(ex1.seq, method = "INDELSLOG",
with.missing = TRUE)
ex1.om <- seqdist(ex1.seq, method = "OM",
indel = costs.ex1$indel, sm = costs.ex1$sm,
with.missing = TRUE)
## Localized OM
ex1.omloc <- seqdist(ex1.seq, method = "OMloc",
sm = costs.ex1$sm, expcost=.1, context = .4,
with.missing = TRUE)
## OMspell with a scalar indel
82 seqdomassoc
indel <- max(costs.ex1$indel)
## OM of spells
ex1.omspell <- seqdist(ex1.seq, method = "OMspell",
indel = indel, sm = costs.ex1$sm,
with.missing = TRUE)
## Distance based on number of matching subsequences
ex1.nms <- seqdist(ex1.seq, method = "NMS",
with.missing = TRUE)
## Using the sequence vectorial representation metric
costs.fut <- seqcost(ex1.seq, method = "FUTURE", lag = 4,
proximities = TRUE, with.missing = TRUE)
ex1.svr <- seqdist(ex1.seq, method = "SVRspell",
prox = costs.fut$prox, with.missing = TRUE)
seqdomassoc Measures of association between domains of sequence data
```

---

## social

```r
data(biofam)
## Building one channel per type of event (left, children or married)
cases <- 1:50
bf <- as.matrix(biofam[cases, 10:25])
children <- bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6
## Building sequence objects
child.seq <- seqdef(children, weights = biofam[cases, 'wp00tbgs '])
marr.seq <- seqdef(married, weights = biofam[cases, 'wp00tbgs '])
left.seq <- seqdef(left, weights = biofam[cases, 'wp00tbgs '])
## distances by channel
dchild <- seqdist(child.seq, method="OM", sm="INDELSLOG")
dmarr <- seqdist(marr.seq, method="OM", sm="INDELSLOG")
dleft <- seqdist(left.seq, method="OM", sm="INDELSLOG")
dbiofam <- list(dchild,dmarr,dleft)
dnames <- names(dbiofam) <- c("child","marr","left")
seqdomassoc(list(child.seq,marr.seq,left.seq), dnames=dnames)
seqdomassoc(list(child.seq,marr.seq,left.seq), diss.dom=dbiofam,
rep.method="rep", assoc="V", dnames=dnames)
seqdomassoc(list(child.seq,marr.seq,left.seq), diss.dom=dbiofam,
rep.method="eq.group", assoc="V", dnames=dnames)
```

---

## seqdur

```r
## Creating a sequence object from columns 13 to 24
## in the 'actcal 'example data set
## Here we retain the first 10 sequences only.
data(actcal)
actcal.seq <- seqdef(actcal[1:10,13:24])
## Retrieving the DSS
actcal.dss <- seqdss(actcal.seq)
## Displaying the DSS for the first 10 sequences
## Example with with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1[, 1:13])
seqdss(ex1.seq)
seqdss(ex1.seq, with.missing=TRUE)
```

---

## seqdss

```r
## Creating a sequence object from columns 13 to 24
## in the 'actcal 'example data set
## Here we retain the first 10 sequences only.
data(actcal)
actcal.seq <- seqdef(actcal[1:10,13:24])
## Retrieving the spell durations
actcal.dur <- seqdur(actcal.seq)
## Displaying the durations for the first 10 sequences
```

---

## use

```r
## Loading data
data(actcal.tse)
## Creating the event sequence object
actcal.eseq <- seqecreate(actcal.tse)
## Printing sequences
## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.eseq,pmin.support=0.01)
## Counting the number of occurrences of each subsequence
msubcount <- seqeapplysub(fsubseq,method="count")
## First lines...
## Presence-absence of each subsequence
msubpres <- seqeapplysub(fsubseq,method="presence")
## First lines...
## Age at first appearance of each subsequence
msubage <- seqeapplysub(fsubseq,method="age")
## First lines...
```

---

## See

```r
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
##Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq <- seqefsub(actcal.eseq, pmin.support=0.01)
##searching for susbsequences discriminating the most men and women
data(actcal)
discr <- seqecmpgroup(fsubseq, group=actcal$sex, method="bonferroni")
##Printing the six most discriminating subsequences
print(discr[1:6])
##Plotting the six most discriminating subsequences
plot(discr[1:6])
```

---

## seqecreate

```r
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
##Searching for frequent subsequences, that is appearing at least 20 times
fsubseq <- seqefsub(actcal.eseq,min.support=20)
##looking for subsequence with FullTime
seqecontain(fsubseq,c("FullTime"))
seqecreate Create event sequence objects.
```

---

## length

```r
##Starting with states sequences
##Loading data
data(biofam)
## Creating state sequences
biofam.seq <- seqdef(biofam,10:25, informat= 'STS')
## Creating event sequences from biofam
biofam.eseq <- seqecreate(biofam.seq)
## Loading data
data(actcal.tse)
## Creating sequences
actcal.eseq <- seqecreate(id=actcal.tse$id, timestamp=actcal.tse$time,
event=actcal.tse$event)
##printing sequences
## Using the data argument
actcal.eseq <- seqecreate(data=actcal.tse)
## Example with missings
data(ex1) ## STS data with missing values
## Creating the state sequence object with by default
## the left missings and gaps coded as '*'and
## end missings coded as void ( '%')
sqex1 <- seqdef(ex1[,1:13])
## and without ignoring right missings (coded as '*')
sqex1b <- seqdef(ex1[,1:13], right=NA)
## Compare the outcome
seqecreate(sqex1)
seqecreate(sqex1, tevent= 'state ')
seqecreate(sqex1, tevent= 'state ', end.event=attr(sqex1, 'void '))
seqecreate(sqex1b, tevent= 'state ')
```

---

## is.seqelist

```r
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
## Searching for subsequences appearing at least 20 times
fsubseq <- seqefsub(actcal.eseq, min.support=20)
## The same using a percentage
fsubseq <- seqefsub(actcal.eseq, pmin.support=0.01)
## Getting a string representation of subsequences
## First ten most frequent subsequences
## Using time constraints
## Looking for subsequences starting in Summer (between June and September)
fsubseq <- seqefsub(actcal.eseq, min.support=10,
constraint=seqeconstraint(age.min=6, age.max=9))
##Looking for subsequences occurring in Summer (between June and September)
fsubseq <- seqefsub(actcal.eseq, min.support = 10,
constraint=seqeconstraint(age.min=6, age.max=9, age.max.end=9))
##Looking for subsequence enclosed in a 6 month period
## and with a maximum gap of 2 month
fsubseq <- seqefsub(actcal.eseq, min.support=10,
constraint=seqeconstraint(max.gap=2, window.size=6))
```

---

## Matthias

```r
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
sid <- seqeid(actcal.eseq)
length(sid)
head(sid)
```

---

## Matthias

```r
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
## Since end.event is not specified, contains no sequence lengths
## We set them manually as 12 for all sequences
sl <- numeric()
sl[1:2000] <- 12
seqelength(actcal.eseq) <- sl
## Retrieve lengths
slen <- seqelength(actcal.eseq)
summary(slen)
seqetm Create a transition-definition matrix
```

---

## for

```r
## Creating a state sequence object from columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24,
labels=c("FullTime", "PartTime", "LowPartTime", "NoWork"))
## Creating a transition matrix, one event per transition
seqetm(actcal.seq,method = "transition")
## Creating a transition matrix, single to-state events
seqetm(actcal.seq,method = "state")
## Creating a transition matrix, two events per transition
seqetm(actcal.seq,method = "period")
## changing the prefix of period start event.
seqetm(actcal.seq,method = "period", bp="begin")
```

---

## Matthias

```r
##Starting with states sequences
##Loading data
data(biofam)
## Creating state sequences
biofam.seq <- seqdef(biofam,10:25,informat= 'STS')
## Creating event sequences from biofam
biofam.eseq <- seqecreate(biofam.seq, weighted=FALSE)
## Using the weights
seqeweight(biofam.eseq) <- biofam$wp00tbgs
## Now seqefsub accounts for weights unless weighted is set to FALSE
fsubseq <- seqefsub(biofam.eseq, pmin.support=0.01)
## Searching for weighted susbsequences which best
## discriminate the birth cohort
discr <- seqecmpgroup(fsubseq, group=biofam$birthyr>=1940)
plot(discr[1:15])
```

---

## Alexis

```r
data(mvad)
mvad.shortlab <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, states=mvad.shortlab, 15:86)
## Finding occurrences of sequence 176 in mvad.seq
seqfind(mvad.seq[176,],mvad.seq)
## Finding occurrences of sequence 1 to 8 in mvad.seq
seqfind(mvad.seq[1:8,],mvad.seq)
seqformat 103
seqformat Conversion between sequence representation formats
```

---

## seqdef

```r
## ========================================
## Examples with raw STS sequences as input
## ========================================
## Loading a data frame with sequence data in the columns 13 to 24
data(actcal)
## Converting to SPS format
actcal.SPS.A <- seqformat(actcal, 13:24, from = "STS", to = "SPS")
head(actcal.SPS.A)
seqformat 107
## Converting to compressed SPS format with no
## prefix/suffix and with "/" as state/duration separator
actcal.SPS.B <- seqformat(actcal, 13:24, from = "STS", to = "SPS",
compress = TRUE, SPS.out = list(xfix = "", sdsep = "/"))
head(actcal.SPS.B)
## Converting to compressed DSS format
actcal.DSS <- seqformat(actcal, 13:24, from = "STS", to = "DSS",
compress = TRUE)
head(actcal.DSS)
## ==============================================
## Examples with a state sequence object as input
## ==============================================
## Loading a data frame with sequence data in the columns 10 to 25
data(biofam)
## Limiting the number of considered cases to the first 20
biofam <- biofam[1:20, ]
## Creating a state sequence object
biofam.labs <- c("Parent", "Left", "Married", "Left/Married",
"Child", "Left/Child", "Left/Married/Child", "Divorced")
biofam.short.labs <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
biofam.seq <- seqdef(biofam, 10:25, alphabet = 0:7,
states = biofam.short.labs, labels = biofam.labs)
## Converting to SPELL format
bf.spell <- seqformat(biofam.seq, from = "STS", to = "SPELL",
pdata = biofam, pvar = c("idhous", "birthyr"))
head(bf.spell)
## Converting to shifted replicated sequences (SRS)
bf.srs <- seqformat(biofam, var=10:25, from="STS", to="SRS",
covar=c("sex","plingu02"))
tail(bf.srs)
## ======================================
## Examples with SPELL sequences as input
## ======================================
## Loading two data frames: bfspell20 and bfpdata20
## bfspell20 contains the first 20 biofam sequences in SPELL format
## bfpdata20 contains the IDs and the years at which the
## considered individuals were aged 15
data(bfspell)
## Converting to STS format with alignement on calendar years
bf.sts.y <- seqformat(bfspell20, from = "SPELL", to = "STS",
id = "id", begin = "begin", end = "end", status = "states",
process = FALSE)
head(bf.sts.y)
## Converting to STS format with alignement on ages
bf.sts.a <- seqformat(bfspell20, from = "SPELL", to = "STS",
id = "id", begin = "begin", end = "end", status = "states",
process = TRUE, pdata = bfpdata20, pvar = c("id", "when15"),
limit = 16)
names(bf.sts.a) <- paste0("a", 15:30)
head(bf.sts.a)
## ==================================
## Examples for TSE and SPELL output
## in presence of missing values
## ==================================
data(ex1) ## STS data with missing values
## creating the state sequence object with by default
## the end missings coded as void ( '%')
sqex1 <- seqdef(ex1[,1:13])
as.matrix(sqex1)
## Creating state-event transition matrices
ttrans <- seqetm(sqex1, method= 'transition ')
tstate <- seqetm(sqex1, method= 'state ')
## Converting into time stamped events
seqformat(sqex1, from = "STS", to = "TSE", tevent = ttrans)
seqformat(sqex1, from = "STS", to = "TSE", tevent = tstate)
## Converting into vertical spell data
seqformat(sqex1, from = "STS", to = "SPELL", with.missing=TRUE)
seqformat(sqex1, from = "STS", to = "SPELL", with.missing=TRUE, right=NA)
seqformat(sqex1, from = "STS", to = "SPELL", with.missing=FALSE)
```

---

## Alexis

```r
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
## Searching for the first occurrence of state 1
## in each of the 5 first sequence of the biofam data set.
seqfpos(biofam.seq[1:5,],"1")
```

---

## Gilbert

```r
## End of spell in further education (FE) in first 10 mvad sequences
data(mvad)
m.seq <- seqdef(mvad[1:10,17:86])
seqfposend(m.seq, state="FE")
```

---

## seqhasmiss

```r
seq <- seqgen(1000,10,1:4,c(0.2,0.1,0.3,0.4))
seqstatd(seqdef(seq))
```

---

## Gilbert

```r
data(ex1)
s <- seqdef(ex1[,1:12])
res <- seqhasmiss(s)
```

---

## seqici

```r
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
seqibad(sx) ## using original alphabet order
seqibad(sx, stprec=c(1,2,3,6)) ## user defined undesirableness values
seqibad(sx, with.missing=TRUE, state.order=c( 'A','B','C','D'))
```

---

## For

```r
## Creating a sequence object from the mvad data set
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
"joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")
mvad.seq <- seqdef(mvad, 15:86, states=mvad.scodes, labels=mvad.labels)
##
mvad.ci <- seqici(mvad.seq)
summary(mvad.ci)
hist(mvad.ci)
## Example using with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13)
seqici(ex1.seq)
seqici(ex1.seq, with.missing=TRUE)
```

---

## by

```r
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Summarize and plots an histogram
## of the within sequence entropy
actcal.ient <- seqient(actcal.seq)
summary(actcal.ient)
hist(actcal.ient)
## Examples using with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
seqient(ex1.seq)
seqient(ex1.seq, with.missing=TRUE)
```

---

## seqintegr

```r
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
print(sx, format= 'SPS')
seqindic(sx, indic=c("lgth","nonm","visited","turbn","cplx"))
seqindic(sx, indic=c("lgth","nonm","visited","turbn","cplx"), with.missing=TRUE)
seqindic(sx, indic=c("lgth","dlgth","ppos","integr","prec"), with.missing=TRUE,
ipos.args=list(pos.states=c("A","B")),
prec.args=list(state.order=c("A","B","C"), state.equiv=list(c("C","D"))))
seqindic(sx, indic=c("volat","binary"), ipos.args=list(pos.states=c("A","B")))
seqindic(sx, indic=c("basic","integr"), ipos.args=list(pos.states="D"))
```

---

## seqipos

```r
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
seqintegr(sx)
seqintegr(sx, with.missing=TRUE)
seqintegr(sx, state="B")
seqintegr(sx, state="B", pow=1.5)
```

---

## Research

```r
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
nr <- attr(sx, 'nr') ## code for missing values
seqipos(sx, pos.states=c("A","B"))
seqipos(sx, pos.states=c("A","B"), with.missing=TRUE)
## ignoring state "D"
seqipos(sx, pos.states=c("A","B"), neg.states=c("C",nr), with.missing=TRUE)
seqipos(sx, pos.states=c("A","B"), dss=FALSE)
seqipos(sx, pos.states=c("A","B",nr), dss=FALSE, with.missing=TRUE)
seqipos(sx, pos.states=c("A","B"), index="volatility")
seqipos(sx, pos.states=c("A","B"), index="integr")
## retrieving the binary sequences
ip <- seqipos(sx, pos.states=c("A","B"))
attr(ip,"sbinary")
ip <- seqipos(sx, pos.states=c("A","B"), with.missing=TRUE)
attr(ip,"sbinary")
```

---

## Research

```r
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
seqistatd(actcal.seq[1:10,])
## Example using "with.missing" argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
seqistatd(ex1.seq)
seqistatd(ex1.seq, with.missing=TRUE)
```

---

## seqdur

```r
sq.dat <- c( 'AAAA ','AAAC ','ABC','ABAA ','AC')
sq <- seqdef(seqdecomp(sq.dat, sep= ''), right=NA)
seqivardur(sq, type=1)
seqivardur(sq, type=1, with.missing=TRUE)
vd2 <- seqivardur(sq, type=2)
print(vd2, stat= 'all')
vd2m <- seqivardur(sq, type=2, with.missing=TRUE)
print(vd2m, stat=c( 'var','vmax '))
```

---

## seqlegend

```r
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
seqivolatility(sx)
seqivolatility(sx, adjust=FALSE)
seqivolatility(sx, with.missing=TRUE)
seqivolatility(sx, w=.7, with.missing=TRUE)
seqlegend Plot a legend for the states in a sequence object
```

---

## Alexis

```r
## Loading the 'actcal 'example data set
## and defining a sequence object with
## activity statuses from Jan. to Dec. 2000 (columns 13 to 24)
## of first 100 cases.
data(actcal)
actcal.seq <- seqdef(actcal[1:100,13:24],
labels=c("> 37 hours", "19-36 hours", "1-18 hours", "no work"))
## Plotting the sequences frequency,
## the states distribution
## and the legend
par(mfrow=c(2,2))
seqIplot(actcal.seq, sortv="from.start", with.legend=FALSE)
seqfplot(actcal.seq, pbarw=TRUE, with.legend=FALSE)
seqdplot(actcal.seq, with.legend=FALSE)
seqlegend(actcal.seq)
```

---

## Alexis

```r
## Loading the 'famform 'example data set
data(famform)
## Defining a sequence object with the 'famform 'data set
ff.seq <- seqdef(famform)
## Retrieving the length of the sequences
## in the ff.seq sequence object
seqlength(ff.seq)
```

---

## Gilbert

```r
## Using the ex1 data set with sequences of different length
data(ex1)
s1 <- seqdef(ex1[,1:13])
seqlength(s1)
## sequence object s2 with a shorter 1st sequence
s2 <- s1
s2[1,8:13] <- attr(s2,"void")
seqlength(s2)
## aligning sequence lengths
seqlength.align(list(s1,s2))
seqLLCP Compute the length of the longest common prefix of two sequences
```

---

## Technical

```r
data(famform)
famform.seq <- seqdef(famform)
## The LCP 's length between sequences 1 and 2
## in the famform sequence object is 2
seqLLCP(famform.seq[1,],famform.seq[2,])
seqLLCS Compute the length of the longest common subsequence of two se-
```

---

## Technical

```r
LCS.ex <- c("S-U-S-M-S-U", "U-S-SC-MC", "S-U-M-S-SC-UC-MC")
LCS.ex <- seqdef(LCS.ex)
seqLLCS(LCS.ex[1,],LCS.ex[3,])
```

---

## Matthias

```r
## Creating the sequence objects using weigths
data(biofam)
biofam.seq <- seqdef(biofam, 10:25, weights=biofam$wp00tbgs)
## Computing sequence probabilities
biofam.prob <- seqlogp(biofam.seq)
## Comparing the probability of each cohort
cohort <- biofam$birthyr>1940
boxplot(biofam.prob~cohort)
```

---

## Gilbert

```r
data(biofam)
b.lab <- c("Parent",
"Left",
"Married",
"Left+Marr",
"Child",
"Left+Child",
"Left+Marr+Child",
"Divorced"
seqMD 135
)
b.short <- c("P","L","M","L+M","C","L+C","L+M+C","D")
set.seed(5)
cases <- sample(nrow(biofam),100)
b.seq <- seqdef(biofam[cases,10:25], labels=b.lab, states=b.short,
weights=biofam[cases,"wp00tbgs"])
## Tokens occurring at least once on average
alphabet(b.seq)[seqmaintokens(b.seq, mint=1)]
## Three more frequent tokens
main.tokens <- seqmaintokens(b.seq, k=3)
## Labels of main tokens
attr(b.seq,"labels")[main.tokens]
## Colors of main tokens
cpal(b.seq)[main.tokens]
seqMD Multidomain sequences
```

---

## seqcost

```r
data(biofam)
## Building one channel per type of event left home, married, and child
cases <- 200
bf <- as.matrix(biofam[1:cases, 10:25])
left <- bf==1 | bf==3 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
children <- bf==4 | bf==5 | bf==6
## Building sequence objects
left.seq <- seqdef(left)
marr.seq <- seqdef(married)
child.seq <- seqdef(children)
channels <- list(LeftHome=left.seq, Marr=marr.seq, Child=child.seq)
## CAT multidomain distances based on channel specific cost methods
MDdist <- seqMD(channels, method="OM",
sm =list("INDELSLOG", "INDELSLOG", "TRATE"), what="diss")
## Providing channel specific substitution costs
smatrix <- list()
smatrix[[1]] <- seqsubm(left.seq, method="TRATE")
smatrix[[2]] <- seqsubm(marr.seq, method="CONSTANT")
smatrix[[3]] <- seqsubm(child.seq, method="CONSTANT")
## Retrieving the MD sequences
MDseq <- seqMD(channels)
alphabet(MDseq)
## Retrieving the CAT multidomain substitution costs
## Using double weight for domain "child"
CATcost <- seqMD(channels,
sm=smatrix, cweight=c(1,1,2), what="cost")
## OMspell distances between MD sequences
MDdist2 <- seqdist(MDseq, method="OMspell",
sm = CATcost, indel=attr(CATcost,"indel"))
```

---

## argument

```r
## Defining a sequence object with columns 13 to 24
## in the actcal example data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal,13:24,labels=actcal.lab)
## Computing the mean time in the different states
seqmeant(actcal.seq)
## Mean times with their standard error
seqmeant(actcal.seq, serr=TRUE)
```

---

## plot.stslist.modst

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)
## Modal state sequence
seqmodst(biofam.seq)
## Examples using weights and with.missing arguments
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
seqmodst(ex1.seq)
seqmodst(ex1.seq, weighted=FALSE)
seqmodst(ex1.seq, weighted=FALSE, with.missing=TRUE)
```

---

## seqLLCP

```r
data(famform)
famform.seq <- seqdef(famform)
seqmpos(famform.seq[1,],famform.seq[2,])
seqmpos(famform.seq[2,],famform.seq[4,])
## Example with gaps in sequences
a <- c(NA,"A",NA,"B","C")
b <- c(NA,"C",NA,"B","C")
ex1.seq <- seqdef(rbind(a,b))
seqmpos(ex1.seq[1,], ex1.seq[2,])
seqmpos(ex1.seq[1,], ex1.seq[2,], with.missing=TRUE)
```

---

## seqdef

```r
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## The first 10 sequences in the actcal.seq
## sequence object
alphabet(actcal.seq)
## The first 10 sequences in the actcal.seq
## sequence object with numerical alphabet
seqnum(actcal.seq[1:10,])
## states A,B,C,D are now coded 0,1,2,3
alphabet(seqnum(actcal.seq))
seqpcplot Parallel coordinate plot for sequence data
```

---

## seqplot

```r
## ================
## plot biofam data
## ================
data(biofam)
lab <- c("Parent","Left","Married","Left+Marr","Child","Left+Child",
"Left+Marr+Child","Divorced")
## plot state sequences in STS representation
## ==========================================
## creating the weighted state sequence object.
biofam.seq <- seqdef(data = biofam[,10:25], labels = lab,
weights = biofam$wp00tbgs)
## select the first 20 weighted sequences (sum of weights = 18)
biofam.seq <- biofam.seq[1:20, ]
par(mar=c(4,8,2,2))
seqpcplot(seqdata = biofam.seq, order.align = "time")
## .. or
seqplot(seqdata = biofam.seq, type = "pc", order.align = "time")
## Distinct successive states (DSS)
## ==========================================
seqplot(seqdata = biofam.seq, type = "pc", order.align = "first")
## .. or (equivalently)
biofam.DSS <- seqdss(seqdata = biofam.seq) # prepare format
seqpcplot(seqdata = biofam.DSS)
## plot event sequences
## ====================
biofam.eseq <- seqecreate(biofam.seq, tevent = "state") # prepare data
## plot the time in the x-axis
seqpcplot(seqdata = biofam.eseq, order.align = "time", alphabet = lab)
## ordering of events
seqpcplot(seqdata = biofam.eseq, order.align = "first", alphabet = lab)
## ... or
plot(biofam.eseq, order.align = "first", alphabet = lab)
## additional arguments
## ====================
## non-embeddable sequences
seqpcplot(seqdata = biofam.eseq, ltype = "non-embeddable",
order.align = "first", alphabet = lab)
## align on last event
par(mar=c(4,8,2,2))
seqpcplot(seqdata = biofam.eseq, order.align = "last", alphabet = lab)
## use group variables
seqpcplot(seqdata = biofam.eseq, group = biofam$sex[1:20],
order.align = "first", alphabet = lab)
## color patterns (Parent)-(Married) and (Parent)-(Left+Marr+Child)
par(mfrow = c(1, 1))
seqpcplot(seqdata = biofam.eseq,
filter = list(type = "sequence",
value=c("(Parent)-(Married)",
"(Parent)-(Left+Marr+Child)")),
alphabet = lab, order.align = "first")
## color subsequence pattern (Parent)-(Left)
seqpcplot(seqdata = biofam.eseq,
filter = list(type = "subsequence",
value = "(Parent)-(Left)"),
alphabet = lab, order.align = "first")
## color sequences over 10% (within group) (function method)
seqpcplot(seqdata = biofam.eseq,
filter = list(type = "function",
value = "minfreq",
level = 0.1),
alphabet = lab, order.align = "first", seed = 1)
## .. same result using the convenience functions
seqpcplot(seqdata = biofam.eseq,
filter = 0.1,
alphabet = lab, order.align = "first", seed = 1)
seqpcplot(seqdata = biofam.eseq,
filter = seqpcfilter("minfreq", 0.1),
alphabet = lab, order.align = "first", seed = 1)
## highlight the 50% most frequent sequences
seqpcplot(seqdata = biofam.eseq,
filter = list(type = "function",
value = "cumfreq",
level = 0.5),
alphabet = lab, order.align = "first", seed = 2)
## .. same result using the convenience functions
seqpcplot(seqdata = biofam.eseq,
filter = seqpcfilter("cumfreq", 0.5),
alphabet = lab, order.align = "first", seed = 2)
148 seqplot
## linear gradient
seqpcplot(seqdata = biofam.eseq,
filter = list(type = "function",
value = "linear"),
alphabet = lab, order.align = "first", seed = 2)
seqpcplot(seqdata = biofam.eseq,
filter = seqpcfilter("linear"),
alphabet = lab, order.align = "first", seed = 1)
seqplot Plot state sequence objects
```

---

## plot.stslist.rep

```r
## ======================================================
## Creating state sequence objects from example data sets
## ======================================================
## biofam data set
data(biofam)
seqplot 153
## We use only a sample of 300 cases
set.seed(10)
biofam <- biofam[sample(nrow(biofam),300),]
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)
## actcal data set
data(actcal)
## We use only a sample of 300 cases
set.seed(1)
actcal <- actcal[sample(nrow(actcal),300),]
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal,13:24,labels=actcal.lab)
## ex1 using weights
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
## ====================
## Sequence index plots
## ====================
## First ten sequences
seqiplot(biofam.seq)
## All sequences sorted by age in 2000
## grouped by sex
seqIplot(actcal.seq, group=actcal$sex, sortv=actcal$age00)
## =======================
## State distribution plot
## =======================
## biofam grouped by sex
seqplot(biofam.seq, type="d", group=biofam$sex)
## actcal grouped by sex
seqplot(actcal.seq, type="d", group=actcal$sex)
## with overlayed entropy line
seqplot(actcal.seq, type="dH", group=actcal$sex)
## ===================
## Cross-sectional entropy plot
## ===================
seqplot(biofam.seq, type="Ht", group=biofam$sex)
## ========================
## Sequence frequency plots
## ========================
154 seqplot
## Plot of the 10 most frequent sequences
seqplot(biofam.seq, type="f")
## Grouped by sex
seqfplot(actcal.seq, group=actcal$sex)
## Unweighted vs weighted frequencies
seqfplot(ex1.seq, weighted=FALSE)
seqfplot(ex1.seq, weighted=TRUE)
## =====================
## Modal states sequence
## =====================
seqplot(biofam.seq, type="ms")
## same as
seqmsplot(biofam.seq)
## ====================
## Representative plots
## ====================
## Computing a distance matrix
## with OM metric
costs <- seqcost(actcal.seq, method="INDELSLOG")
actcal.om <- seqdist(actcal.seq, method="OM", sm=costs$sm, indel=costs$indel)
## Plot of the representative sets grouped by sex
## using the default density criterion
seqrplot(actcal.seq, group=actcal$sex, diss=actcal.om, coverage=.5)
## Plot of the representative sets grouped by sex
## using the "dist" (centrality) criterion
seqrplot(actcal.seq, group=actcal$sex, criterion="dist", diss=actcal.om, coverage=.33)
## ========================
## Relative frequency plots
## ========================
## Using default sorting by first MDS variable
seqrfplot(actcal.seq, diss=actcal.om, sortv=NULL, group=actcal$sex)
## ===============
## Mean time plot
## ===============
## actcal data set, grouped by sex
seqplot(actcal.seq, type="mt", group=actcal$sex)
## displaying mean times as bar labels
group <- factor(actcal$sex)
blab <- NULL
for (i in 1:length(levels(group))){
blab <- cbind(blab,seqmeant(actcal.seq[group==levels(group)[i],]))
}
seqmtplot(actcal.seq, group=group,
bar.labels = round(blab,digits=2), cex.barlab=1.2)
```

---

## seqpcplot

```r
data(biofam)
## Building one channel per type of event (left, children or married)
cases <- 200
bf <- as.matrix(biofam[1:cases, 10:25])
children <- bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6
## Building sequence objects
require(colorspace)
m.col <- sequential_hcl(2, palette = "PurpOr")
c.col <- sequential_hcl(2, palette = "Mint")
l.col <- sequential_hcl(2, palette = "OrYel")
child.seq <- seqdef(children, weights=biofam[1:cases,"wp00tbgs"], cpal=c.col)
marr.seq <- seqdef(married, weights=biofam[1:cases,"wp00tbgs"], cpal=m.col)
left.seq <- seqdef(left, weights=biofam[1:cases,"wp00tbgs"], cpal=l.col)
seqdom <- list(LeftHome=left.seq,Marr=marr.seq,Child=child.seq)
seqplotMD(seqdom, type="d", group=biofam[1:cases,"sex"], xaxis="bottom")
seqplotMD(seqdom, type="dH", group=biofam[1:cases,"sex"], dom.byrow=TRUE,
xaxis="bottom", yaxis="left")
seqplotMD(seqdom, type="mt", group=biofam[1:cases,"sex"])
seqplotMD(seqdom, type="I", group=biofam[1:cases,"sex"],
xaxis="bottom", sortv="from.end")
## sorting on first domain
seqplotMD(seqdom, type="I", group=biofam[1:cases,"sex"],
xaxis="bottom", sortv="from.start", dom.crit=1)
seqplotMD(seqdom, type="f", group=biofam[1:cases,"sex"],
xaxis="bottom", yaxis="left")
## distances betweem MD sequences
MDseq <- seqMD(seqdom, what="MDseq", ch.sep="+")
diss <- seqdist(MDseq, method="OM", sm="INDELSLOG")
seqplotMD(seqdom, type="rf", group=biofam[1:cases,"sex"],
xaxis="bottom", sortv="from.end", dom.crit=0, diss=diss, k=10)
seqplotMD(seqdom, type="r", group=biofam[1:cases,"sex"],
xaxis="bottom", dom.crit=0, diss=diss)
```

---

## Alexis

```r
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## search for pattern "DAAD"
## (no work-full time work-full time work-no work)
## results are stored in the 'daad 'object
daad <- seqpm(actcal.seq,"DAAD")
## Looking at the sequences
## containing the pattern
## search for pattern "AD"
## (full time work-no work)
seqpm(actcal.seq,"AD")
```

---

## seqici

```r
## Defining a sequence object with columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal <- actcal[1:20,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24], alphabet=c( 'A','B','C','D'))
## insecurity and precarity using the original state order
insec <- seqinsecurity(actcal.seq)
prec <- seqprecarity(actcal.seq)
ici <- seqici(actcal.seq) ## complexity
seqn <- seqformat(actcal.seq, to="SPS", compress=TRUE)
tab <- data.frame(seqn,ici,insec,prec)
names(tab) <- c("seq","ici","insec","prec")
## Assuming A and B as equivalent regarding insecurity
insec2 <- seqinsecurity(actcal.seq, state.equiv=list(c( 'A','B')))
tab <- cbind(tab,insec2)
names(tab)[ncol(tab)] <- "insec2"
## and letting C be non-comparable
insec3 <- seqinsecurity(actcal.seq, state.order=c("A","B","D"),
state.equiv=list(c( 'A','B')))
tab <- cbind(tab,insec3)
names(tab)[ncol(tab)] <- "insec3"
## bounding insecurity with undesirableness of best and worst state in sequence
insec4 <- seqinsecurity(actcal.seq, state.order=c("A","B","D"),
state.equiv=list(c( 'A','B')), bound=TRUE)
tab <- cbind(tab,insec4)
names(tab)[ncol(tab)] <- "insec4"
## degrading index
degr <- seqidegrad(actcal.seq, state.order=c("A","B","D"),
state.equiv=list(c( 'A','B')))
tab <- cbind(tab,degr)
names(tab)[ncol(tab)] <- "degr"
## Precarity with transition weights based on differences in state undesirableness
prec.rank <- seqprecarity(actcal.seq, method= 'RANK ')
###############
## Indexes in presence of missing values:
## missing state treated as an additional state
data(ex1)
## by default right missings are dropped from the sequences
s.ex1 <- seqdef(ex1[,1:13])
state.order=c("A","B","C","D") ## missing left as non-comparable
seqprecarity(s.ex1, state.order=state.order, with.missing=TRUE)
seqinsecurity(s.ex1, state.order=state.order, with.missing=TRUE)
## same result using the correction argument
dgp <- seqidegrad(s.ex1, with.missing=TRUE, state.order=state.order, method= 'TRATEDSS ')
seqprecarity(s.ex1, state.order=state.order, with.missing=TRUE, correction=dgp)
dg <- seqidegrad(s.ex1, with.missing=TRUE, state.order=state.order)
seqinsecurity(s.ex1, state.order=state.order, with.missing=TRUE, correction=dg)
## bounding insecurity with undesirableness of best and worst state in sequence
seqinsecurity(s.ex1, state.order=state.order, with.missing=TRUE, bound=TRUE)
```

---

## seqprecarity

```r
## Defining a sequence object with columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal <- actcal[1:200,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24])
## state precarity levels using the original state order
seqprecstart(actcal.seq)
## Assuming A and B as equivalent regarding precarity
seqprecstart(actcal.seq, state.equiv=list(c( 'A','B')))
## and letting C be non-comparable
seqprecstart(actcal.seq, state.order=c("A","B","D"), state.equiv=list(c( 'A','B')))
```

---

## seqdef

```r
## Recoding a state sequence object with seqrecode
data(actcal)
## Creating a state sequence object
actcal.seq <- seqdef(actcal,13:24, labels=c("> 37 hours", "19-36 hours",
"1-18 hours", "no work"))
## Regrouping states B and C and setting the whole alphabet to A BC D
actcal.new <-seqrecode(actcal.seq,
recodes = list("A"="A", "BC"=c("B", "C"), "D"="D"))
## Crosstabulate the first column of the recoded and
## original state sequence objects
table(actcal.new[,1], actcal.seq[,1])
## Same as before but using automatically original
## codes for unspecified states.
actcal.new2 <-seqrecode(actcal.seq,
recodes = list("BC"=c("B", "C")))
table(actcal.new2[,1], actcal.seq[,1])
## Same as before but using otherwise
actcal.new3 <-seqrecode(actcal.seq, recodes = list("A"="A", "D"="D"),
otherwise="BC")
table(actcal.new3[,1], actcal.seq[,1])
## Recoding factors
## Recoding the marital status to oppose married to all other case
maritalstatus <- recodef(actcal$civsta00,
recodes=list("Married"="married"), otherwise="Single")
summary(maritalstatus)
table(maritalstatus, actcal$civsta00)
## Recoding the number of kids in the household
## -2 is a missing value
nbkids <- recodef(actcal$nbkid00,
recodes=list("None"=0, "One"=1, "Two or more"=2:10), na=-2)
table(nbkids, actcal$nbkid00, useNA="always")
```

---

## seqplot

```r
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam[,10:25], labels=biofam.lab)
## Computing the distance matrix
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", sm=costs)
## Representative set using the neighborhood density criterion
biofam.rep <- seqrep(biofam.seq, diss=biofam.om, criterion="density")
summary(biofam.rep)
plot(biofam.rep)
## plot by groups represented by the representatives
seqdplot(biofam.seq, group=attr(biofam.rep,"Rep.group"), border=NA)
## indexes of sequences represented by 1st representative
r1.grp <- which(attr(biofam.rep,"Rep.group")==1)
## indexes of occurrences of the first representative sequence
attr(biofam.rep,"idx.rep")[[1]]
```

---

## seqdecomp

```r
seqsep("ABAAAAAAD")
```

---

## seqici

```r
## Loading the 'actcal 'example data set
data(actcal)
## Here we consider only the first 10 sequences
actcal <- actcal[1:10,]
## Defining a sequence object with data in columns 13 to 24
## (activity status from January to December 2000)
actcal.seq <- seqdef(actcal[,13:24], informat= 'STS')
## Computing the sequences turbulence
turb <- seqST(actcal.seq)
## Normalized turbulence
turb.norm <- seqST(actcal.seq, norm=TRUE)
## Normalized turbulence taking non-visited states into account.
turb2.norm <- seqST(actcal.seq, norm=TRUE, type=2)
```

---

## seqdHplot

```r
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
sd <- seqstatd(biofam.seq)
## Plotting the state distribution
plot(sd, type="d")
## Line of cross-sectional entropies
plot(sd, type="Ht")
## ====================
## example with weights
## ====================
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
## Unweighted
seqstatd(ex1.seq, weighted=FALSE)
seqstatd(ex1.seq, weighted=TRUE)
seqstatf State frequencies in the whole sequence data set
```

---

## within

```r
## Creating a sequence object from the actcal data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)
## States frequencies
seqstatf(actcal.seq)
## Example with weights
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
## Unweighted
seqstatf(ex1.seq, weighted=FALSE)
## Weighted
seqstatf(ex1.seq, weighted=TRUE)
seqstatl List of distinct states or events (alphabet) in a sequence data set.
```

---

## theTraMineR

```r
data(actcal)
seqstatl(actcal,13:24)
```

---

## seqdss

```r
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Number of subsequences with DSS=TRUE
seqsubsn(actcal.seq[1:10,])
## Number of subsequences with DSS=FALSE
seqsubsn(actcal.seq[1:10,],DSS=FALSE)
```

---

## Sequences

```r
## Creating a sequence object from the actcal data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)
## 10 most frequent sequences in the data
seqtab(actcal.seq)
## With idxs=0, we get all distinct sequences in the data set
## sorted in decreasing order of their frequency
stab <- seqtab(actcal.seq, idxs=0)
head(stab)
tail(stab)
## Example with weights
## from biofam data set using weigths
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
## Unweighted frequencies
seqtab(ex1.seq, weighted=FALSE)
## Weighted frequencies
seqtab(ex1.seq, weighted=TRUE)
```

---

## seqdss

```r
## Creating a sequence object from columns 13 to 24
## in the 'actcal 'example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Computing the number of transitions
actcal.trans <- seqtransn(actcal.seq)
## Displaying number of transitions in the first 10 sequences
## Example with with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13)
seqtransn(ex1.seq)
seqtransn(ex1.seq, with.missing=TRUE)
```

---

## seqdist

```r
## Loading the 'actcal 'example data set
data(actcal)
## Defining a sequence object with data in columns 13 to 24
## (activity status from January to December 2000)
actcal.seq <- seqdef(actcal[,13:24])
## Computing transition rates
seqtrate(actcal.seq)
## Computing transition rates between states "A" and "B" only
seqtrate(actcal.seq, c("A","B"))
## ====================
## Example with weights
## ====================
data(ex1)
ex1.seq <- seqdef(ex1[,1:13], weights=ex1$weights)
seqtrate(ex1.seq, weighted=FALSE)
seqtrate(ex1.seq, weighted=FALSE, count=TRUE)
## weights are accounted for by default
seqtrate(ex1.seq)
seqtrate(ex1.seq, count=TRUE)
```

---

## seqtreedisplay

```r
data(mvad)
## Defining a state sequence object
mvad.seq <- seqdef(mvad[, 17:86])
## Growing a seqtree from Hamming distances:
## Warning: The R=10 used here to save computation time is
## much too small and will generate strongly unstable results.
## We recommend to set R at least as R=1000.
## To comply with this small R value, we set pval = 0.1.
seqt <- seqtree(mvad.seq~ male + Grammar + funemp + gcse5eq + fmpr + livboth,
data=mvad, R=10, pval=0.1, seqdist.arg=list(method="HAM", norm="auto"))
print(seqt)
## Growing a seqtree from an existing distance matrix
mvad.dhd <- seqdist(mvad.seq, method="DHD")
seqt <- seqtree(mvad.seq~ male + Grammar + funemp + gcse5eq + fmpr + livboth,
data=mvad, R=10, pval=0.1, diss=mvad.dhd)
print(seqt)
### Following commands only work if GraphViz is properly installed
## Not run:
seqtreedisplay(seqt, type="d", border=NA)
seqtreedisplay(seqt, type="I", sortv=cmdscale(mvad.dhd, k=1))
## End(Not run)
```

---

## For

```r
## Creating a sequence object with the columns 13 to 24
## in the 'actcal 'example data set
## The color palette is automatically set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
## Retrieving the color palette
stlab(actcal.seq)
seqiplot(actcal.seq)
## Changing the state labels
stlab(actcal.seq) <- c("Full time","Part time (19-36 hours)",
"Part time (1-18 hours)", "No work")
seqiplot(actcal.seq)
```
