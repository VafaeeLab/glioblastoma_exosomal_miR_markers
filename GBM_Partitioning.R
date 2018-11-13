library(limma)
library(edgeR)
library(pROC)
library(caret)
library(randomForest)
library(DMwR)
library(gplots)
library(ggplot2)
library(e1071)
library(devtools)
library(ggbiplot)
library(robustbase)

# reading tabular data into R

pheno <- read.table("./Phenotype.txt",  header = TRUE, fill = TRUE, row.names = "Sample", sep = "\t")
raw   <- read.table("./Normalized_microRNA_counts.txt", header = TRUE, row.names = "miRNA",fill = TRUE,sep = "\t")

# Subsetting the input files:

set1 = subset(pheno, GBMvsCont == 'Control')
set2 = subset(pheno, GBMvsCont == 'GBM')

#set1 = subset(pheno, GliomavsMCont == 'Control')
#set2 = subset(pheno,  GliomavsMCont == 'Glioma')

set = data.frame(rep(0,dim(raw)[2]))
colnames(set) = "comparator"
rownames(set) = colnames(raw)
set[rownames(set1),1] = 1
set[rownames(set2),1] = 2
x.raw <- raw[ ,which(set > 0)]  
groups <- set[which(set > 0),]

# Filtering parameters

minLibraryCount = 50000
minNumberSamples = 0.5*ncol(x.raw)
minCount = 100

# Normalisation & Filtering: low abundant miRNAs with <200 read counts across 50% of samples were removed

keep <- colSums(x.raw) >= minLibraryCount  
x.filt <- x.raw[keep]
g <- groups[keep]
x.cpm <- cpm(x.filt) 
keep <- rowSums(x.cpm >=minCount) >= minNumberSamples
x.norm <- x.cpm[keep,] # 123 miRs, 29 samples


# Random Partitioning -----------------------------------------------------

source("./getSigMiRs.R")
N = 100
y    = as.factor(g)
data = x.norm

fit.rf = list()
DE.miRs = list()
DE.miRs.tets = list()
all.miRs = c()

acc = matrix(data = NA,N,1) 
sen = matrix(data = NA,N,1)
spe = matrix(data = NA,N,1)

for (i in 1:N) {
  print(i) 
  
# Seperate data to test and train
  
Train <- createDataPartition(y, p=0.5, list=FALSE)
  
#Identify features (DE miRs)

data.train = data[,Train]
y.train = y[Train]
DE.miRs[[i]] <- rownames(getSigMiRs(data.train, y.train, use.FDR = FALSE, pval_Cutoff = 0.05, logFC_cutoff = log2(2)))
  
  # fold change and p-value of DE.miRs on Test samples

  data.test = data[,-Train]
  y.test = y[-Train]
  DE.miRs.tets[[i]] <- getSigMiRs(data.test, y.test, use.FDR = FALSE, pval_Cutoff = 0.05, logFC_cutoff = log2(2))
  
  # Set up expressions for classifiers
  x = t(data[DE.miRs[[i]],])
  all.miRs = c(all.miRs,DE.miRs[[i]])
  
  # Run Random Forest
  
  x.train <- x[ Train, ]
  x.test <-  x[ -Train, ]
  y.train <- y[Train]
  y.test <-  y[-Train]
  
  fit.rf[[i]] <- randomForest(x.train,as.factor(y.train), xtest =x.test, ytest = as.factor(y.test), importance=TRUE, proximity = TRUE)
  t.rf = fit.rf[[i]]$test$confusion
   
  acc[i,1] = sum(diag(t.rf))/sum(t.rf)
  sen[i,1] = t.rf[1,1]/rowSums(t.rf[,1:2])[1]
  spe[i,1] = t.rf[2,2]/rowSums(t.rf[,1:2])[2]
  
}


t = table(all.miRs)
outFile = paste("results",".GBMvsCont-partitioning.csv", sep="")
write.csv(t, file =outFile, row.names=TRUE,quote = FALSE)
hist(t,20)
keep = t >= 75
t = t[keep]
res = matrix(data = NA, nrow = N,ncol = length(t))
for (i in 1:N) {
  imp.rf = importance(fit.rf[[i]], type = 1)
  for (j in 1:length(t)) {
    ind = which(rownames(imp.rf) == rownames(t)[j])
    if (length(ind) > 0) {
      res[i,j] = imp.rf[ind] 
    }
  }
}

colnames(res) = rownames(t)

heatmap.2(t(res),dendrogram = "none", Rowv = FALSE, Colv = "Rowv", density.info='none',
          main = 'Feature Importance',trace = "none",scale = "none",
          col = colorpanel(100, "yellow", "white", "purple"), 
          cexRow = 1.25, na.rm = FALSE, na.color = "black")

res.fc = matrix(data = NA, nrow = length(t),ncol = N)
rownames(res.fc) = rownames(t)

for (i in 1:length(t)) {
  miR = rownames(t)[i]
  for (j in 1:N) {
    tmp = DE.miRs.tets[[j]]
    res.fc[i,j] = tmp[miR,"logFC"] # change if using voom version
  }
}

fc = rowMedians(res.fc, na.rm = TRUE)
fc = as.matrix(fc)
rownames(fc) = rownames(t)
heatmap.2(cbind(fc,fc),dendrogram = "none", Rowv = FALSE, Colv = "Rowv", density.info='none',main = 'fold change', trace = "none",scale = "none",col = colorpanel(100, "blue", "white", "red"),cexRow = 1.25, na.rm = FALSE, na.color = "black")

#,colsep = 1:ncol(res), rowsep = 1:nrow(res), sepcolor="white",sepwidth=c(0.05,0.05))


res.pval = matrix(data = NA, nrow = length(t),ncol = N)
rownames(res.pval) = rownames(t)

for (i in 1:length(t)) {
  miR = rownames(t)[i]
  for (j in 1:N) {
    tmp = DE.miRs.tets[[j]]
    m = min(tmp[miR,c("ttest.pval","wilcox.pval","edgeR.pval")])
    res.pval[i,j] = m
  }
}
pval = rowMedians(res.pval, na.rm = TRUE)
pval = as.data.frame(pval)
rownames(pval) = rownames(t)
barplot(-log10(rowMedians(res.pval, na.rm = TRUE)), main="-Log10Pvalues", horiz=FALSE, beside = TRUE, names.arg=rownames(pval))

#Randomforest on selected miRs using all data
z = data.frame(cbind(t(data[colnames(res),]),g))

#select optimum mtry 
tuneRF(z[,1:ncol(z) -1], as.factor(z$g), mtryStart = 3, trace=TRUE, plot=TRUE, doBest=FALSE)

rf <- randomForest(as.factor(g) ~ ., data = z, mtry = 6, importance=TRUE, proximity = TRUE)
print(rf)

# Importance of each predictor.
varImpPlot(rf, type = 1)

# Try all combinations, Violin plot ---------------------------------------
z = data.frame(cbind(t(data[rownames(t),]),g))
err <- c()
s <- c()
miRs <- c()
n = 2^length(t)-1
for (i in 1:n) {
  print(i)
  v = as.integer(intToBits(i))[1:length(t)]
  zz = subset(z,select=c(as.logical(v),TRUE))
  rf <- randomForest(as.factor(g) ~ ., data = zz, importance=TRUE, proximity = TRUE)
  err = c(err,rf$err.rate[50,1])
  s = c(s,length(which(v==1)))
  tmp = zz
  st = paste(colnames(tmp),collapse=", ")
  miRs = c(miRs,st)
  cat(i,"\t",rf$err.rate[50,1],"\t",colnames(zz),"\n")
}

es = as.data.frame(cbind(err,s))

ggplot(es, aes(factor(s), err)) + geom_violin(aes(fill = factor(s)))+geom_jitter(size = 1.5,width = 0.2)#+coord_flip()


