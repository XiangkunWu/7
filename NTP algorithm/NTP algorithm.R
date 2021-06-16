# set working path
workdir <- "C:\\Users\\18718\\Desktop\\新建文件夹"; setwd(workdir)

#-----------#
# using ntp #
library(limma)
library(CMScaller)

# set label
dirct <- "up" # find up-regulated markers, change to 'dn' if you want to use down-regulated biomarkers (not suggested)
p.cutoff <- 0.05 # cutoff to find marker
p.adj.cutoff <- 0.05 # cutoff to find marker
n.marker <- 100 # find 10 markers for each subtype

# 1. load data (if you provide subtype on your own, please make sure you have the same format of 'subtype' file)
train.exp <- read.table("TCGA_matrix.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
test.exp <- read.table("E-MTAB-1980_matrix.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
subtype <-  read.table("CNMF_subtype.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)

# 2. extract common samples
comsam <- intersect(colnames(train.exp), rownames(subtype))
train.exp <- train.exp[,comsam]
subtype <- subtype[comsam, , drop = F]

# 3. perform differential expression analysis using limma in training dataset
createList  <- function(res = NULL) {
  
  tumorsam <- res$samID
  sampleList <- list()
  treatsamList <- list()
  treatnameList <- c()
  ctrlnameList <- c()
  
  subt <- unique(res$Subtype)
  
  for (i in subt) {
    sampleList[[i]] <- tumorsam
    treatsamList[[i]] <- intersect(tumorsam, res[which(res$Subtype == i), "samID"])
    treatnameList[i] <- i
    ctrlnameList[i] <- "Others"
  }
  
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}

complist <- createList(res = subtype)
sampleList <- complist[[1]]
treatsamList <- complist[[2]]
treatnameList <- complist[[3]]
ctrlnameList <- complist[[4]]
allsamples <- colnames(train.exp)

# log transformation (use range() to see if your own data needs log transformation)
if(range(train.exp)[2] > 100) {
  gset <- log2(train.exp + 1)
} else {gset <- train.exp}

for (k in 1:length(sampleList)) {
  samples <- sampleList[[k]]
  treatsam <- treatsamList[[k]]
  treatname <- treatnameList[k]
  ctrlname <- ctrlnameList[k]
  
  compname <- paste(treatname, "_vs_", ctrlname, sep="")
  tmp <- rep("others", times = length(allsamples))
  names(tmp) <- allsamples
  tmp[samples] <- "control"
  tmp[treatsam] <- "treatment"
  
  outfile <- file.path(workdir, paste("limma_test_result.", compname, ".txt", sep = ""))
  
  pd <- data.frame(Samples = names(tmp),
                   Group = as.character(tmp),
                   stringsAsFactors = FALSE)
  
  design <-model.matrix(~ -1 + factor(pd$Group, levels = c("treatment","control")))
  colnames(design) <- c("treatment","control")
  
  fit <- limma::lmFit(gset, design = design);
  contrastsMatrix <- limma::makeContrasts(treatment - control, levels = c("treatment", "control"))
  fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
  fit2 <- limma::eBayes(fit2, 0.01)
  resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
  resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
  resData$id <- rownames(resData)
  colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
  resData$fc <- 2^resData$log2fc
  resData <- resData[order(resData$padj),c("id","fc","log2fc","pvalue","padj")]
  write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
  cat(paste0("limma of ",compname, " done...\n"))
}

# 4. find subtype-specific up-regulated biomarkers
DEpattern <- ".*._vs_Others.txt$"
DEfiles <- dir(workdir,pattern = DEpattern)

if(dirct == "up") { outlabel <- "unique_upexpr_marker.txt" }
if(dirct == "down") { outlabel <- "unique_downexpr_marker.txt" }

genelist <- c()
for (filek in DEfiles) {
  DEres <- read.table(file.path(workdir, filek), header=TRUE, row.names=NULL, sep="\t", quote="", stringsAsFactors=FALSE)
  DEres <- DEres[!duplicated(DEres[, 1]),]
  DEres <- DEres[!is.na(DEres[, 1]), ]
  rownames(DEres) <- DEres[, 1]
  DEres <- DEres[, -1]
  
  #rownames(DEres) <- toupper(rownames(DEres))
  if (dirct == "up") {
    genelist <- c( genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]) )
  }
  if (dirct == "down") {
    genelist <- c( genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]) )
  }
}
unqlist <- setdiff(genelist,genelist[duplicated(genelist)])

marker <- list()
for (filek in DEfiles) {
  DEres <- read.table(file.path(workdir, filek), header=TRUE, row.names=NULL, sep="\t", quote="", stringsAsFactors=FALSE)
  DEres <- DEres[!duplicated(DEres[, 1]),]
  DEres <- DEres[!is.na(DEres[, 1]), ]
  rownames(DEres) <- DEres[, 1]
  DEres <- DEres[, -1]
  
  #rownames(DEres) <- toupper(rownames(DEres))
  if(dirct == "up") {
    outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]) )
    outk <- DEres[outk,]
    outk <- outk[order(outk$log2fc, decreasing = TRUE),]
    
    if(nrow(outk) > n.marker) {
      marker[[filek]] <- outk[1:n.marker,]
    } else {
      marker[[filek]] <- outk
    }
    marker$dirct <- "up"
  }
  if(dirct == "down") {
    outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]) )
    outk <- DEres[outk,]
    outk <- outk[order(outk$log2fc, decreasing = FALSE),]
    
    if(nrow(outk) > n.marker) {
      marker[[filek]] <- outk[1:n.marker,]
    } else {
      marker[[filek]] <- outk
    }
    marker$dirct <- "down"
  }
  # write file
  write.table(outk, file=file.path(workdir, paste(gsub("_vs_Others.txt","", filek, fixed = TRUE), outlabel, sep = "_")), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
}

# generate templates for nearest template prediction
templates <- NULL
for (filek in DEfiles) {
  tmp <- data.frame(probe = rownames(marker[[filek]]),
                    class = sub("_vs_Others.txt","",sub(".*.result.","",filek)),
                    dirct = marker$dirct,
                    stringsAsFactors = FALSE)
  templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
}
write.table(templates, file=file.path(workdir,paste0(dirct,"regulated_marker_templates.txt")), row.names = FALSE, sep = "\t", quote = FALSE)

# 5. apply template to training set to see agreement
emat <- t(scale(t(gset), scale = T, center = T))
ntp.train <- ntp(emat      = emat,
                 templates = templates,
                 doPlot    = T,
                 nPerm     = 1000,
                 distance  = "cosine",
                 nCores    = 1,
                 seed      = 2020104,
                 verbose   = T)
dev.copy2pdf(file = "ntp heatmap in train.pdf", width = 6, height = 6)
subtype$ntp <- ntp.train[rownames(subtype), "prediction"]
table(subtype$Subtype, subtype$ntp)
#   C1  C2
C1 166  10
C2  54 279
write.table(subtype, "train subtype with ntp prediction.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(ntp.train, "ntp prediction in training.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 6. apply template to testing set
emat <- t(scale(t(test.exp), scale = T, center = T)) # test.exp has been already normalized so no log transformation will do
ntp.test <- ntp(emat      = emat,
                templates = templates,
                doPlot    = T,
                nPerm     = 1000,
                distance  = "cosine",
                nCores    = 1,
                seed      = 202095,
                verbose   = T)
dev.copy2pdf(file = "ntp heatmap in test.pdf", width = 6, height = 6)
write.table(ntp.test, "ntp prediction in testing.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
