library(limma)
library(edgeR)
library(calibrate)
library(data.table)
library(gplots)
library(ggplot2)
library(e1071)
library(randomForest)
library(gridExtra)
library(boot)
library(pROC)
library(RColorBrewer)

# Read files: 

pheno <- read.table("./Phenotype.txt",  header = TRUE, fill = TRUE, row.names = "Sample", sep = "\t")
raw   <- read.table("./Normalized_microRNA_counts.txt", header = TRUE, row.names = "miRNA",fill = TRUE,sep = "\t")


# Subsetting the input files:

set1 = subset(pheno, GBMvsCont == 'Control'    )
set2 = subset(pheno, GBMvsCont == 'GBM'  )

#set1 = subset(pheno, GliomavsMCont == 'Control'    )
#set2 = subset(pheno, GliomavsMCont == 'Glioma'  )

#set1 = subset(pheno, GBMvsGlioma == 'Glioma'    )
#set2 = subset(pheno, GBMvsGlioma == 'GBM'  )

set = data.frame(rep(0,dim(raw)[2]))
colnames(set) = "comparator"
rownames(set) = colnames(raw)
set[rownames(set1),1] = 1
set[rownames(set2),1] = 2
x.raw <- raw[ ,which(set > 0)]  
groups <- set[which(set > 0),]


# Normalisation & Filtering -----------------------------------------------

minLibraryCount = 50000
minNumberSamples = 0.5*ncol(x.raw)
minCount = 100

keep <- colSums(x.raw) >= minLibraryCount  
x.filt <- x.raw[keep]
g <- groups[keep]

x.cpm <- cpm(x.filt) 

keep <- rowSums(x.cpm >= minCount) >= minNumberSamples
x <- x.cpm[keep,] 

# Testing for DE miRNAs ---------------------------------------------------

# Fold Change calculation
fc1 = data.frame(round(rowMeans(x[,g==2])/rowMeans(x[,g==1]),digits = 1))
fc = log2(fc1)
colnames(fc) <-"logFC"  

# 1) EdgR:
y <- DGEList(counts=x,group=factor(g)) 
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
de <- exactTest(y)
tmp <- topTags(de,n=Inf)$table
de.EdgR <- tmp[ order(row.names(tmp)), ]

# 2) TTEST:
de.ttest <- c()
for (i in 1:nrow(x)) {
  de.ttest[i] <- t.test(x = x[i,g==1], y = x[i,g==2], var.equal = FALSE)$p.value
}
de.ttest <- data.frame(de.ttest)
row.names(de.ttest) <-rownames(x)
colnames(de.ttest) <-"PValue"

# 3) Wilcoxon:
de.wilcox <- c()
for (i in 1:nrow(x)) {
  de.wilcox[i] <- wilcox.test(x = x[i,g==1], y = x[i,g==2])$p.value
}
de.wilcox <- data.frame(de.wilcox)
row.names(de.wilcox) <-rownames(x)
colnames(de.wilcox) <-"PValue"

# Output 
res <- cbind(de.EdgR[,"PValue"],de.ttest,de.wilcox,fc,rowMeans(x[,g==2]),rowMeans(x[,g==1]))
colnames(res) <- c("edgeR.PValue", "ttest.PValue","wilcox.PValue", "FC","mean.CPM.GBM","mean.CPM.Cont")
outFile = paste("results",".GBMvsControl.csv", sep="")
write.csv(res, file =outFile, row.names=TRUE,quote = FALSE)


# VennDiagram -------------------------------------------------------------
pval_Cutoff = 0.05
logFC_cutoff = log2(2)
A = row.names(subset(res, edgeR.PValue < pval_Cutoff))
B = row.names(subset(res, ttest.PValue < pval_Cutoff))
C = row.names(subset(res, wilcox.PValue < pval_Cutoff))
D = row.names(subset(res, abs(fc[,1])>= logFC_cutoff))
venn( list(EdgeR=A,TTest=B,Wilcox=C) )

de <- unique(union(union(intersect(A,B),intersect(B,C)),intersect(A,C)))
de.fc <- intersect(de,D)

# Heatmap -----------------------------------------------------------------
x.de <- x[c(de.fc),]
x.de.reorder <-cbind(x.de[,g==2],x.de[,g==1])
category.color<-c(rep("red",length(g[g==2])),rep("blue",length(g[g==1])))


ColSideColors = category.color
heatmap.2(x.de.reorder,dendrogram = "row", Colv = FALSE,density.info='none', colCol = category.color, trace = "none",scale = "row",col = colorpanel(100, "dodgerblue", "white", "red"), cexRow = 1.25, lmat = rbind(c(0,3,4),c(2,1,0),c(0,0,0)), lwid = c(0.5,4,1.5), lhei = c(1,1,1), sepwidth=c(0.01,0.01),sepcolor="white",colsep=1:ncol(x.de.reorder),rowsep=1:nrow(x.de.reorder))

# Boxplots ----------------------------------------------------------------

p = list()
for (i in 1:nrow(x.de)) {
  z = cbind(g,x.de[i,])
  colnames(z) = c("Groups","CMP")
  z = data.frame(z)
#  p[[i]]<-ggplot(z, aes(factor(g), expr, fill = factor(g)))+ geom_boxplot()+scale_fill_manual(values=alpha(c("yellow", "green"), 0.7))+ geom_jitter(width = 0.1)+ggtitle(rownames(x.de)[i]) + theme(panel.background = element_rect(fill = "grey95"), plot.title = element_text(size=15),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")#+ labs(x="group1 vs group2",y="miR Expression") 
  p[[i]]<-ggplot(z, aes(factor(Groups), CMP))+ geom_boxplot(fatten = 0.8, lwd = 0.25,width = 0.5)+ geom_jitter(width = 0.1)+ggtitle(rownames(x.de)[i]) +
    scale_x_discrete(labels = c("Glio", "GBM"))+
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          plot.title = element_text(size=15),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position="none")
  
  #+ labs(x="group1 vs group2",y="miR Expression") 
  
  }
do.call(grid.arrange,p)


# GLM and ROC -------------------------------------------------------------

par(mfrow=c(5, 7))

glmFit = list()
p.glm = list()
cv = list()

ci = list() 

err = c()
for (c in 1:nrow(x.de)) {
  z = cbind(g==1,x.de[c,])
  colnames(z) = c("g","expr")
  z = data.frame(z)
  glmFit[[c]] <- glm(g~expr,data=z,family=binomial())
  
  prob=predict(glmFit[[c]],type=c("response"))
  z$prob = prob
  r = roc(g ~ prob, data = z)
  
  tmp = ci.auc(roc(g ~ prob, data = z, smooth = FALSE), conf.level = 0.95)
  ci[[c]] = tmp[1:3]
  
  plot(r,main = paste(rownames(x.de)[c]," AUC: ", round(r$auc,3)),xlim =c(1:0),ylim = c(0:1))
  cv[[c]] <- cv.glm(z, glmFit[[c]]) # leav-one-out CV
  err = c(err, cv[[c]]$delta[1])
}

ci.m <- as.data.frame(matrix(unlist(ci), ncol = 3, byrow = TRUE))
ci.m = cbind(cbind(rownames(x.de)),ci.m)
colnames(ci.m) <- c("miR", "low", "auc", "high")

View(ci.m)
#View(ci.m) to see AUC and lower/higher 95% confidence interval

z = cbind(as.data.frame(rownames(x.de)),err)
colnames(z) = c("miR","err")
ggplot(z, aes(x = factor(miR), y = err)) + geom_bar(aes(fill = factor(miR)), colour = "grey50", width = 0.5, stat = "identity")+  
 scale_fill_manual(values = colorRampPalette(brewer.pal(6, "YlOrRd"))(nrow(z))) +
 coord_cartesian(ylim = c(0.10, 0.25))+
  theme(panel.background = element_rect(fill = "white", color = "grey"), 
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none")

outFile = paste("LOOCV",".ControlvsLGG.csv", sep="")
write.csv(z, file =outFile, row.names=TRUE,quote = FALSE)

# RandomForest ------------------------------------------------------------
z = data.frame(cbind(t(x.de),g))

#select optimum mtry 
tuneRF(z[,1:ncol(z) -1], as.factor(z$g), mtryStart = 3, trace=TRUE, plot=TRUE, doBest=FALSE)

rf <- randomForest(as.factor(g) ~ ., data = z, mtry = 6, importance=TRUE, proximity = TRUE)
print(rf)

# Importance of each predictor.
varImpPlot(rf, type = 1)

zz = subset(z,select=c("hsa.miR.182.5p", "hsa.miR.328.3p", "hsa.miR.485.3p", "hsa.miR.486.5p", "g"))
##zz = subset(z,select=c("hsa.miR.182.5p", "hsa.miR.339.5p", "hsa.miR.340.5p", "hsa.miR.485.3p", "hsa.miR.486.5p", "g"))
##zz = subset(z,select=c("hsa.miR.328.3p", "hsa.miR.485.3p" , "g"))
##zz = subset(z,select=c("hsa.miR.339.5p", "hsa.miR.485.3p" , "g"))
rf <- randomForest(as.factor(g) ~ ., data = zz, importance=TRUE, proximity = TRUE)
tmp = ci.auc(roc(g,as.numeric(rf$predicted), data = z, smooth = FALSE), conf.level = 0.95)
ci[[c]] = tmp[1:3] # this will give you AUC and 95% CI
View(ci[[c]])

