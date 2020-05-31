#################################################################
## Direct inhibition of the NOTCH transcription factor complex ##
#################################################################
## Autor: Sara Dorado Alfaro                                   ##
#################################################################

### Working directory
mypath="~/Documents/Master_bioinfo/TRRPP/Trabajo_GSEA"
setwd(mypath)
data_path = "GSE18198_RAW/"

##1. Load libraries
library(affy) #para affymetrix
library(limma) #obtaining DEGs
library(genefilter) #filtering
library(annotate) #annotation
library(hgu133plus2.db) #annotation human genome
library(ggplot2) #plotting

## carga de datos
#2. Import targets.txt file
targets <- readTargets("GSE18198_RAW/targets.txt", row.names="FileName")
#3. Import .CEL files
# linea celular KOPT-K1
targets.KOPT_K1 <- targets[targets$LineaCelular=="KOPT-K1", ]
files.KOPT_K1 <- paste(data_path, paste(targets.KOPT_K1$Fichero, ".CEL", sep=""), sep="")
data.KOPT_K1 <- ReadAffy(filenames=files.KOPT_K1)
# linea velular HPB-ALL
targets.HPB_ALL <- targets[targets$LineaCelular=="HPB-ALL", ]
files.HPB_ALL<- paste(data_path, paste(targets.HPB_ALL$Fichero, ".CEL", sep=""), sep="")
data.HPB_ALL <- ReadAffy(filenames=files.HPB_ALL)
#phenodata for plots
ph.KOPT_K1 = data.KOPT_K1@phenoData
ph.HPB_ALL = data.HPB_ALL@phenoData

## Preprocesamiento
#4. Normalize with RMA 
#generates object eset (class ExprSet), 
#expresso function provides intensities in log scale
# norm_data.KOPT_K1 <- expresso(data.KOPT_K1,
#                         bg.correct=TRUE,
#                         bgcorrect.method="rma",
#                         normalize=TRUE,
#                         normalize.method="quantiles",
#                         pmcorrect.method="pmonly",
#                         summary.method="medianpolish",
#                         verbose=TRUE)
# norm_data.HPB_ALL <- expresso(data.HPB_ALL,
#                         bg.correct=TRUE,
#                         bgcorrect.method="rma",
#                         normalize=TRUE,
#                         normalize.method="quantiles",
#                         pmcorrect.method="pmonly",
#                         summary.method="medianpolish",
#                         verbose=TRUE)
# save results
#save(norm_data.KOPT_K1, file="norm_KOPT-K1.RData")
#save(norm_data.HPB_ALL, file="norm_HPB-ALL.RData")
# load data
load("norm_KOPT-K1.RData")
load("norm_HPB-ALL.RData")

#boxplots
boxplot(data.KOPT_K1, main="Datos crudos (KOPT-K1)",
        ylab="Intensidad", col="lightgreen")
boxplot(data.HPB_ALL, main="Datos crudos (HPB-ALL)",
        ylab="Intensidad", col="lightgreen")
boxplot(exprs(norm_data.KOPT_K1), main="Datos normalizados (KOPT-K1)",
        ylab="Intensidad (log scale)", col="lightgreen")
boxplot(exprs(norm_data.HPB_ALL), main="Datos normalizados (HPB-ALL)",
        ylab="Intensidad (log scale)", col="lightgreen")

#histograms
color=c('green','green','green','red','red','red')
hist(data.KOPT_K1[,1:6],lwd=2,which='pm',col=color,
     ylab='Densidad',xlab='Log2 intensidad',main='Datos crudos KOPT-K1')
hist(data.HPB_ALL[,1:6],lwd=2,which='pm',col=color,
     ylab='Densidad',xlab='Log2 intensidad',main='Datos crudos HPB-ALL')

#6.Data filtering using IQR.
esetIQR.KOPT_K1 <- varFilter(norm_data.KOPT_K1, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
esetIQR.HPB_ALL <- varFilter(norm_data.HPB_ALL, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)

#######Análisis de expresión diferencial#######
#7. Design matrix.
design.KOPT_K1<-cbind(DMSO=c(1,1,1,0,0,0), 
              SAHM1=c(0,0,0,1,1,1)) #generamos el diseño a mano
rownames(design.KOPT_K1)<-targets.KOPT_K1$Fichero
design.HPB_ALL<-cbind(DMSO=c(1,1,1,0,0,0), 
                      SAHM1=c(0,0,0,1,1,1)) #generamos el diseño a mano
rownames(design.HPB_ALL)<-targets.HPB_ALL$Fichero

#8. Contrasts matrix.
#prepara contrastes estadísticos asignando 
#que muestras van en cada contraste
cont_matrix.KOPT_K1 <- makeContrasts(SAHM1vsDMSO=SAHM1-DMSO, levels=design.KOPT_K1)
cont_matrix.HPB_ALL <- makeContrasts(SAHM1vsDMSO=SAHM1-DMSO, levels=design.HPB_ALL)

#9. Obtaining differentially expressed genes (DEGs)
#Linear model and eBayes 
#KOPT-K1
fit.KOPT_K1 <- lmFit(esetIQR.KOPT_K1, design.KOPT_K1)
fit2.KOPT_K1 <- contrasts.fit(fit.KOPT_K1, cont_matrix.KOPT_K1)
fit2.KOPT_K1 <- eBayes(fit2.KOPT_K1)
# Results (KOPT-K1)
toptableIQR.KOPT_K1 <- topTable(fit2.KOPT_K1, number=dim(exprs(esetIQR.KOPT_K1))[1], adjust.method="BH", sort.by="p")
#HPB-ALL
fit.HPB_ALL <- lmFit(esetIQR.HPB_ALL, design.HPB_ALL) 
fit2.HPB_ALL <- contrasts.fit(fit.HPB_ALL, cont_matrix.HPB_ALL)
fit2.HPB_ALL <- eBayes(fit2.HPB_ALL)
# Results (HPB-ALL)
toptableIQR.HPB_ALL <- topTable(fit2.HPB_ALL, number=dim(exprs(esetIQR.HPB_ALL))[1], adjust.method="BH", sort.by="p")

##10. Filter DEGs (FDR < 0.005)
filtered.KOPT_K1 <- subset(toptableIQR.KOPT_K1, adj.P.Val<=0.001)
filtered.HPB_ALL <- subset(toptableIQR.HPB_ALL, adj.P.Val<=0.001)


####### Annotation #######
##11. Annotate
list_GeneSymbol.KOPT_K1 <- getSYMBOL(rownames(filtered.KOPT_K1), "hgu133plus2.db")
toptable_annotated.KOPT_K1 <- data.frame(cbind(filtered.KOPT_K1, GeneSymbol=list_GeneSymbol.KOPT_K1))
toptable_annotated.KOPT_K1$Probe <- rownames(toptable_annotated.KOPT_K1)

data_plot = toptable_annotated.KOPT_K1[order(toptable_annotated.KOPT_K1$logFC),][1:18,]

ggplot(data_plot, aes(y=logFC, x=GeneSymbol)) + 
  geom_col(position="dodge") +
  #scale_x_discrete(limit = head(KOPT_down$plot_label,21)) +
  coord_flip() +
  ggtitle("Genes diferencialmente subregulados en KOPT-K1 con tratamiento SAHM1")

list_GeneSymbol.HPB_ALL <- getSYMBOL(rownames(filtered.HPB_ALL), "hgu133plus2.db")
toptable_annotated.HPB_ALL <- cbind(filtered.HPB_ALL, GeneSymbol=list_GeneSymbol.HPB_ALL)
toptable_annotated.HPB_ALL$Probe <- rownames(toptable_annotated.HPB_ALL)

data_plot = toptable_annotated.HPB_ALL[order(toptable_annotated.HPB_ALL$logFC),][1:18,]

ggplot(data_plot, aes(y=logFC, x=GeneSymbol)) + 
  geom_col(position="dodge") +
  #scale_x_discrete(limit = head(KOPT_down$plot_label,21)) +
  coord_flip() +
  ggtitle("Genes diferencialmente subregulados en HPB-ALL con tratamiento SAHM1")

## 13. Common downregulated genes between 'KOPT-K1' and 'HPB-ALL'
common_genes <- merge(toptable_annotated.KOPT_K1[toptable_annotated.KOPT_K1$logFC < 0, c("Probe", "GeneSymbol", "logFC")],
                      toptable_annotated.HPB_ALL[toptable_annotated.HPB_ALL$logFC < 0, c("Probe", "GeneSymbol", "logFC")],
                      by=c('Probe', 'GeneSymbol'), suffixes=c('(KOPT-K1)', '(HPB_ALL)'))
common_genes

#Plot results
KOPT_subset <- common_genes[,1:3]
HPB_subset <- common_genes[,c(1,2,4)]
colnames(KOPT_subset) <- c("Probe", "GeneSymbol", "logFC")
colnames(HPB_subset) <-  c("Probe", "GeneSymbol", "logFC")
KOPT_subset$Cell_line <- 'KOPT-K1'
HPB_subset$Cell_line <- 'HPB-ALL'
HPB_subset
KOPT_subset
common_rbind <- rbind(KOPT_subset, HPB_subset)
ggplot() + geom_col(data=common_rbind, aes(x=GeneSymbol, y=logFC, fill=Cell_line),
                    position='dodge') +
  ggtitle("DEGs for KOPT-K1 and HPB-ALL after SAHM1 treatment") +
  theme(plot.title=element_text(hjust=0.5))

##13. Check against GSI-NOTCH Geneset
GSI_geneset <- read.csv(paste(data_path,"GSI-NOTCH_gene_set.csv", sep=""),
                        header=TRUE,
                        sep=";")
GSI_geneset
KOPT_subset <- subset(toptable_annotated.KOPT_K1, Probe %in% GSI_geneset$Probe) #nothing, which is good
HPB_subset <- subset(toptable_annotated.HPB_ALL, Probe %in% GSI_geneset$Probe) #nothing, which is good
# pintamos los resultados (GeneSymbol vs logFC)
KOPT_subset$Cell_line <- 'KOPT-K1'
HPB_subset$Cell_line <- 'HPB-ALL'
down_rbind <- rbind(KOPT_subset, HPB_subset)
ggplot() + geom_col(data=down_rbind, aes(x=GeneSymbol, y=logFC, fill=Cell_line),
                    position='dodge') +
  ggtitle("DEGs for KOPT-K1 and HPB-ALL after SAHM1 treatment (Intersection with GSI genes)") +
  theme(plot.title=element_text(hjust=0.5))



