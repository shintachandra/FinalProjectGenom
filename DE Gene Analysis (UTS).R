if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("Biobase", "GEOquery"))
update.packages("GEOquery")  # Update paket GEOquery ke versi terbaru

library(Biobase)
library(GEOquery)
library(hgu133a.db)
library(limma)
library("annotate")
library("hgu133a.db") 


dtgeo <- getGEO('GSE10072', destdir=".", GSEMatrix = TRUE)[[1]]
print(dtgeo)

eset <- dtgeo
print(eset)

phdtgeo <- pData(eset)
head(phdtgeo)
names(phdtgeo)

expdtgeo <- exprs(eset)
dim(expdtgeo)
head(expdtgeo)
annotation(dtgeo)

#Meta(dtgeo)$platform

annotation(eset)  <- "hgu133a"

#anotasi untuk microarray yang digunkana
annotation(eset) <- "hgu133a.db"

#install microarray
BiocManager::install("hgu133a.db", ask=FALSE)
library(hgu133a.db)

#filterisasi data 
BiocManager::install ("genefilter")
require(genefilter)

esetFilt = nsFilter(eset) 
print(esetFilt)


### Extract the Expression of the Filtered Dataset ###
expdtgeoFilt <- exprs(esetFilt$eset)

###  plot the  original and filtered data ###
par(mfrow=c(1,2))
hist(expdtgeo, main ='original')
hist(expdtgeoFilt, main="filtered")

dim(expdtgeoFilt)
View(phdtgeo)

fit
### Analysis ###

vargrp <- phdtgeo$source_name_ch1
table(vargrp)

###  Assign Group Code ###
group <- ifelse(vargrp=="Adenocarcinoma of the Lung", 0 , 1)
group
# 0 =sakit cancer 1= normal

BiocManager::install ("limma")  ## Instal limma ##
###  Limma Analysis###
library(limma)
design <- model.matrix(~group)
fit <- eBayes(lmFit(expdtgeoFilt,design))
fit

## Show only top 50 genes ##
topResult <- topTable(fit, coef=2, number=50)

# Install dan load enhanced volcano 
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)
EnhancedVolcano(topResult,
                lab = rownames(topResult),   # Label dari nama baris (gene names)
                x = 'logFC',                 # Sumbu X: Log2 Fold Change
                y = 'adj.P.Val',             # Sumbu Y: Adjusted p-value
                xlab = bquote(Log[2] ~ "Fold Change"),  # Memperbaiki ekspresi untuk xlab
                ylab = "Adjusted p-value")  # Label untuk sumbu Y

###  Selected Genes###

rownames(topResult)

## Extract selected gene names ##
selected  <- rownames(expdtgeoFilt) %in% rownames(topResult)

## Extract the expression of the selected genes
expdtgeosel <- expdtgeoFilt [selected, ]

###  heatmap of the top genes ###

heatmap(expdtgeosel)
par(mfrow=c(2,2))
for ( i in 1:4) plot(vargrp, expdtgeosel[i,], 
                     main= rownames(expdtgeosel)[i])
str(vargrp)
table(vargrp, useNA = "ifany")
vargrp <- as.factor(vargrp)
summary(expdtgeosel[i, ])
par(mfrow=c(2,2))
for (i in 1:4) {
  y <- as.numeric(expdtgeosel[i, ])  # pastikan y numeric
  if (all(is.finite(y))) {           # hanya plot kalau semua nilai valid
    plot(vargrp, y, main = rownames(expdtgeosel)[i])
  } else {
    message("Skipping ", rownames(expdtgeosel)[i], " because of NA or Inf values")
  }
}


#### See gene name and description  ###

library("annotate")
library("hgu133a.db")    ## Note: It depends on your chip 

GeneSelected <- select(hgu133a.db, rownames(topResult), 
                       c("SYMBOL","ENTREZID", "GENENAME"))
GeneSelected

ids <- rownames(topResult)
ids

GeneSelected <- select(hgu133a.db, ids, c("SYMBOL","ENTREZID", "GENENAME", "GO"))

## Lakukan jika belum meng install GO.db
BiocManager::install ("GO.db")

###  Gene Ontology for the top genes###
library(GO.db)
GOselected <- select(GO.db, GeneSelected$GO, c("TERM", "GOID"))
head(GOselected)

#simpan grup label dan sampel dalam format csv
data_filter_trans <- t(expdtgeoFilt)
write.csv(data_filter_trans, "data_filter1.csv", row.names = TRUE)

#file nyimpan label grup yang benar
group_df <- data.frame(SampleID = colnames(expdtgeoFilt),
                       group = group)
write.csv(group_df, "label_grup_sampel1.csv", row.names = FALSE)