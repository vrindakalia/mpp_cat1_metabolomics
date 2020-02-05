####
# Revisiting cat-1 and MPP+ metabolomics data for IUTOX and for worm paper
####

###########
# N2 worms exposed to MPP+ for 4 hours
# There are 5 biological samples for exposure and control
# M9 threshold = 1.5
###########

############
# HILICPOS #
############

## Ignoring bacterial blank since exposures were made in M9, and worms were bacteria free for nearly 5 hours
setwd("Raw_files/Hilic_pos")

map_all <- read.table("Hilicpos_sample_id_mapfile.txt", header = T, sep = "\t")

feat <- read.table("Merged_pall_averaged_featuretable_repmax2of3.txt", header = T, sep = "\t")

## Subset MPP+ samples
mpp_map_all <- map_all[157:234,]


### Keeping only those ID's that match the averaged sample ID
mpp_map_all$num <- c(rep(1:3))
mpp_map <- mpp_map_all[which(mpp_map_all$num == 1),]

#get rid of .mzXML in column names
#names(feat) <- gsub(pattern = ".mzXML", replacement = "", x = names(feat), fixed = TRUE) 

# Using the map file to extract index of samples from feature table
a <- as.vector(NA)
for(i  in 1:26){
    a[i] <- which(colnames(feat) == mpp_map$File.Name[i])
}

# Extracting samples from feature table
feat_mpp <- feat[,c(1:2,a)]

# CVs of all samples from this run
hist(feat$median_CV, main = "MPP+ CVs", xlab = "Median CV")

# Extarcting MPP+ rows from map file
mpp_samples <- mpp_map[-grep("Q", mpp_map$Sample.ID),]

b <- as.vector(NA)
for(i  in 1:20){
    b[i] <- which(colnames(feat_mpp) == mpp_samples$File.Name[i])
}
feat_mpp_samples <- feat_mpp[,c(1:2,b)]


feat_mpp_samples$mztime <- paste(feat_mpp_samples$mz, feat_mpp_samples$time, sep = "/") ## For metaboanalyst

mpp_metab <- feat_mpp_samples[, c(23, 3:22)] ## For metaboanalyst

mpp_label <- mpp_samples[,c(1:2)]

mpp_label$sample <- unlist(lapply(strsplit(as.character(mpp_label$Sample.ID),split = "_"), function(i) i[1]))

mpp_label$exp <- unlist(lapply(strsplit(as.character(mpp_label$Sample.ID),split = "_"), function(i) i[2]))

## Adding appropriate labels to blank samples
for(i in 1:nrow(mpp_label)){
    if(mpp_label$sample[i]=="M9"){
        mpp_label$exp[i] <- "M9"
    } else if(mpp_label$sample[i] == "Bact"){
        mpp_label$exp[i] <- "Bact"
    }
}

### Adding labels to feature table

### Converting abundance values to numeric
for(i in 1:22){
    feat_mpp_samples[[i]] <- as.numeric(as.character(feat_mpp_samples[[i]]))
}

#label <- c("mz", "time", mpp_map$exp, "mztime")
#label <- c("mz", "time", mpp_map$exp)
#lab <- matrix(NA, nrow = 1, ncol = ncol(mpp_feat))
#colnames(lab) <- names(mpp_feat)
#mpp_featlab <- rbind(lab, mpp_feat)
#mpp_featlab[1,] <- label

### Adding labels to metaboanalyst raw feature file
#metab_lab <- mpp_featlab[1,c(23, 3:22)]
#mpp_metab_lab <- rbind(metab_lab, mpp_metab)

#write.table(mpp_metab_lab, "mpp_metaboanalyst_raw_06232018.txt", col.names = T, row.names = F, sep = "\t")

############
# M9 BLANK #
############

### Using M9 as threshold
# Creating a subset of M9 features
m9_filenames <- mpp_samples[grep("M9", mpp_samples$Sample.ID),]
m <- as.vector(NA)
for(i  in 1:5){
    m[i] <- which(colnames(feat_mpp_samples) == m9_filenames$File.Name[i])
}
m9_feat <- feat_mpp_samples[,c(m)]

for(i in 1:1){
    m9_feat[[i]] <- as.numeric(as.character(m9_feat[[i]]))
}

m9_feat$mean <- rowMeans(m9_feat)

m9_mean <- as.data.frame(t(m9_feat$mean))

mpp_worms <- mpp_label[which(mpp_label$exp == "MPP" | mpp_label$exp == "C"),]
w <- as.vector(NA)
for(i  in 1:10){
    w[i] <- which(colnames(feat_mpp_samples) == mpp_worms$File.Name[i])
}
mpp_worm_feat <- feat_mpp_samples[,c(w)]

tmpp_worm_feat <- as.data.frame(t(mpp_worm_feat))
tmpp_worm_feat[1:5,1:5]
## M9 threshold: 1.5 
## This step takes a long time! 
mat_m9 <- matrix(NA, ncol = ncol(tmpp_worm_feat), nrow=nrow(tmpp_worm_feat))
for(j in 1:ncol(tmpp_worm_feat)){
    for(i in 1: nrow(tmpp_worm_feat)){
        if(tmpp_worm_feat[i,j] > 1.5*m9_mean[1,j]){
            mat_m9[i,j] <- 1
        }
        else {
            mat_m9[i,j] <- 0
        }
    }
}

m9_sum <- colSums(mat_m9)

index.threshold <- which(m9_sum >=1)

mpp_mz <- feat_mpp_samples[,1:2]
mpp_mz_feat <- rbind(t(mpp_mz), tmpp_worm_feat)
mpp_mz_feat[1:5,1:5]

mpp_mz_filt <- mpp_mz_feat[,index.threshold]

mpp_mz_filt[1:5, 1:5]


### Replacing missing/ non-detects with min/2
### Generate minimum values for each mzcat
mz_filt_m9 <- mpp_mz_filt[1:2,]
mpp_feat <- mpp_mz_filt[-c(1:2),]

for(j in 1:ncol(mpp_feat)){
    for(i in 1:nrow(mpp_feat)){
        if(mpp_feat[i,j] == 0 | is.null(mpp_feat[i,j] | is.na(mpp_feat[i,j] | is.nan(mpp_feat[i,j] | is.infinite(mpp_feat[i,j] ))))){
            mpp_feat[i,j] <- (0.5*min(mpp_feat[which(mpp_feat[,j] > 0),j]))
        }
    }
}

mpp_feat[1:5,1:5]
# Log transform abundances to approach normality 
mpp_feat.log <- log10(mpp_feat)

# Checking log-transformation for normality
ggplot(data = mpp_feat, aes(x = V10)) + geom_histogram()
ggplot(data = mpp_feat.log, aes(x = V10)) + geom_histogram()

# Binding mzs back with imputed features
mpp.filt.imp.mz <- rbind(mz_filt_m9, mpp_feat)
mpp.filt.imp.mz.log <- rbind(mz_filt_m9, mpp_feat.log)
# Changing columns and removing mz and time rows
feat.mzs <- paste(mpp.filt.imp.mz[1,], mpp.filt.imp.mz[2,], sep = "_")
names(mpp.filt.imp.mz) <- feat.mzs
names(mpp.filt.imp.mz.log) <- feat.mzs
mpp.filt.imp <- mpp.filt.imp.mz[-c(1:2),]
mpp.filt.imp.log <- mpp.filt.imp.mz.log[-c(1:2),]

# Adding in labels
sample.order <- rownames(mpp_feat)
label <- mpp_worms$exp[order(match(mpp_worms$File.Name, sample.order))]
mpp.filt.imp$label <- c(label)
mpp.filt.imp.log$label <- c(label)

mpp.metaboanalyst <- as.data.frame(t(mpp.filt.imp)) %>%
    .[c(14735, 1:14734),] 

write.table(mpp.metaboanalyst, "/Users/vk2316/Documents/Miller_lab/Metabolomics_data/MPP/20181116/mpp_metaboanalyst.txt", col.names = T, row.names = T, sep = "\t")

# Running multiple regression models, regressing on all the features

# Running regressions with multiple dependent variables
# OUTCOME
out_start=1 # position of first mz
out_end= 14734 # position of last mz
out_nvar=out_end-out_start+1 # total number of dependent variables

out_variable=rep(NA, out_nvar) # matrix to hold all outcome data

# data from label variable
beta=rep(NA, out_nvar) # coeeficient for status
se = rep(NA, out_nvar)
pvalue=rep(NA, out_nvar)


# Use this code for multiple dependent variables
# EXPOSURE 
#exp_start=4
#exp_end=5
#exp_nvar=exp_end-exp_start+1
#exp_nvar = 2
#exp_variable=rep(NA, exp_nvar)
#exp_beta=rep(NA, exp_nvar)
#exp_se = rep(NA, out_nvar)
#exp_pvalue=rep(NA, exp_nvar)

number=1 # counter

for (i in out_start:out_end){
    outcome = colnames(mpp.filt.imp.log)[i]
    model <- lm(get(outcome) ~ label,
                na.action = na.exclude,
                data=mpp.filt.imp.log)
    
    Vcov <- vcov(model, useScale = FALSE)
    beta <- model$coefficients
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    beta[number] = as.numeric(beta[2])
    se[number] = as.numeric(se[2])
    pvalue[number] = as.numeric(pval[2])
    out_variable[number] = outcome
    
    #exp_beta[number] = as.numeric(beta[2])
    #exp_se[number] = as.numeric(se[2])
    #exp_pvalue[number] = as.numeric(pval[2])
    # exp_variable[number] = exposure
    number = number + 1
}

# Creating data frame with relevant regression output
outcome = data.frame(out_variable, beta, se, pvalue)


# Subsetting data based on p-value cutoff
outcome.sub <- outcome %>%
    filter(status_pvalue < 0.1) %>%
    mutate(thresh = as.numeric(status_pvalue < 0.05))

# Plotting subset of mzs to show features significant at p < 0.05
ggplot(data = outcome.sub, aes(x = out_variable, y = pvalue, color = factor(thresh))) + 
    geom_point() + scale_y_reverse() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# total number of features that met threshold of p<0.05
sum(outcome.sub$thresh) # 363 variables with p < 0.05

# Creating separate file for mummichog
outcome.mcg <- outcome.sub %>%
    select(out_variable, status_pvalue, status_beta)



##################################################################
# EDIT THIS CODE NEXT
##################################################################
## PANDA
mpp_m9_panda <- cbind(mz_filt_m9, mpp_filt_m9_feat)
lab <- c("mz", "time", "MPP_1", "Control_1","MPP_2", "Control_2", "MPP_3", "Control_3",  "MPP_4", "Control_4", "MPP_5", "Control_5")
label <-matrix(NA, nrow = 1, ncol = 12)
colnames(label) <- names(mpp_m9_panda)
mpp_m9_pan <- rbind(label, mpp_m9_panda)
mpp_m9_pan[1,] <- lab


lab <- c("mz", "time", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "mppsum", "controlsum")
label <-matrix(NA, nrow = 1, ncol = 14)
colnames(label) <- names(mpp_filt_m9)
mpp_m9_fc <- rbind(label, mpp_filt_m9)
mpp_m9_fc[1,] <- lab

setwd("C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP")

write.table(mpp_m9_pan,"mpp_m9_15_panda.txt", col.names = F, row.names = F, sep = "\t")

## Metaboanalyst
mpp_m9_pan$mztime <- paste(mpp_m9_pan$mz, mpp_m9_pan$time, sep = "/")
mpp_m9_metab <- mpp_m9_pan[,c(13, 3:12)]
lab <- c("mztime", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control")
mpp_m9_metab[1,] <- lab
write.table(mpp_m9_metab,"mpp_m9_15_metab.txt", col.names = T, row.names = F, sep = "\t")

library(xmsPANDA)

feature_table_file<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\mpp_m9_15_panda.txt"
class_labels_file<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\classlabels_m9.txt"
outloc<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\PANDA"

demetabs_res<-diffexp(feature_table_file=feature_table_file,
                      parentoutput_dir=outloc,
                      class_labels_file=class_labels_file,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="zeros",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=NA, group.missing.thresh=NA, input.intensity.scale="raw", 
                      log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE, 
                      quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
                      rsd.filt.list = c(0), pairedanalysis = FALSE, featselmethod="limma",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="pearson", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("royalblue1","darkgrey"), 
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2, 
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="BER",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=FALSE,max_comp_sel=1,saveRda=FALSE,pca.cex.val=4,pls.permut.count=NA,
                      pca.ellipse=TRUE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=5)


sink(file=NULL)
#####################################################

###########
# C18 NEG #
###########
## Ignoring bacterial blank since exposures were made in M9, and worms were bacteria free for nearly 5 hours
setwd("C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\c18_neg")

map_all <- read.table("C18neg_sample_id_mapfile.txt", header = T, sep = "\t")

feat <- read.table("Merged_pall_averaged_featuretable_repmax2of3.txt", header = T, sep = "\t")

### Keeping only those ID's that match the averaged sample ID
map_all$num <- c(rep(1:3))
map <- map_all[which(map_all$num == 1),]

### Keeping only MPP+ labels
map_mpp <- map[c(53:78),] 

#get rid of .mzXML in column names
names(feat) <- gsub(pattern = ".mzXML", replacement = "", x = names(feat), fixed = TRUE) 

# From feature table, extracting the mpp samples
mpp_feat <- feat[,c(1:10,63:88)]

hist(mpp_feat$median_CV, main = "MPP+ CVs", xlab = "Median CV")

mpp_feat <- mpp_feat[,-c(3:10)]

mpp_feat$mztime <- paste(mpp_feat$mz, mpp_feat$time, sep = "/") ## For metaboanalyst

mpp_metab <- mpp_feat[, c(23, 3:22)] ## For metaboanalyst

mpp_map <- map_mpp[,c(1:2)]

mpp_map$sample <- unlist(lapply(strsplit(as.character(mpp_map$Sample.ID),split = "_"), function(i) i[1]))

mpp_map$exp <- unlist(lapply(strsplit(as.character(mpp_map$Sample.ID),split = "_"), function(i) i[2]))

## Adding appropriate labels to blank samples
for(i in 1:nrow(mpp_map)){
    if(mpp_map$sample[i]=="M9"){
        mpp_map$exp[i] <- "M9"
    } else if(mpp_map$sample[i] == "Bact"){
        mpp_map$exp[i] <- "Bact"
    }
}

### Adding labels to feature table

### Converting abundance values to numeric
for(i in 1:22){
    mpp_feat[[i]] <- as.numeric(as.character(mpp_feat[[i]]))
}

#label <- c("mz", "time", mpp_map$exp, "mztime")
label <- c("mz", "time", mpp_map$exp)
lab <- matrix(NA, nrow = 1, ncol = ncol(mpp_feat))
colnames(lab) <- names(mpp_feat)
mpp_featlab <- rbind(lab, mpp_feat)
mpp_featlab[1,] <- label

### Adding labels to metaboanalyst raw feature file
metab_lab <- mpp_featlab[1,c(23, 3:22)]
mpp_metab_lab <- rbind(metab_lab, mpp_metab)

write.table(mpp_metab_lab, "mpp_metaboanalyst_raw_06232018.txt", col.names = T, row.names = F, sep = "\t")

############
# M9 BLANK #
############

### Using M9 as threshold
tmpp_featlab <- as.data.frame(t(mpp_featlab))
tm9 <- tmpp_featlab[grep("M9", tmpp_featlab$V1),]
m9 <- as.data.frame(t(tm9))
m9 <- m9[-1,]

for(i in 1:ncol(m9)){
    m9[[i]] <- as.numeric(as.character(m9[[i]]))
}

m9$mean <- rowMeans(m9[,1:5])
mztime <- mpp_feat[,1:2]
## Extracting both MPP and Control M9 features
mpp_m9 <- tmpp_featlab[which(tmpp_featlab$V1 == "MPP" | tmpp_featlab$V1== "C"),]
## Transposing the dataframe back
mpp_m9 <- as.data.frame(t(mpp_m9))
mpp_m9 <- mpp_m9[-1,]
for(i in 1:10){
    mpp_m9[[i]] <- as.numeric(as.character(mpp_m9[[i]]))
}

## M9 threshold: 1.5
for(j in 1:10){
    for(i in 1: nrow(mpp_m9)){
        if(mpp_m9[i,j] > 1.5*m9$mean[i]){
            mpp_m9[i,j+10] <- 1
        }
        else {
            mpp_m9[i,j+10] <- 0
        }
    }
}

for(i in 1:nrow(mpp_m9)){
    mpp_m9[i,21] <- sum(mpp_m9[i, 11:15])
}

for(i in 1:nrow(mpp_m9)){
    mpp_m9[i,22] <- sum(mpp_m9[i, 16:20])
}

#mpp_m9 <- t(tmpp_m9)
mpp_m9_mz <- cbind(mztime, mpp_m9)

### Keeping features that meet sum >=1 criteria
mpp_filt_m9 <- mpp_m9_mz[which(mpp_m9_mz[,23] >=1 | mpp_m9_mz[,24] >=1),]

mpp_filt_m9 <- mpp_filt_m9[,c(1:12,23,24)]

lab <- c("mz", "time", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "mppsum", "controlsum")
label <-matrix(NA, nrow = 1, ncol = 14)
colnames(label) <- names(mpp_filt_m9)
mpp_m9_fc <- rbind(label, mpp_filt_m9)
mpp_m9_fc[1,] <- lab

### Replacing missing/ non-detects with min/2
## Replace zeros in cat_filt1 with min/2
## removing the sum columns from cat_filt1
mpp_filt_m9_feat <- mpp_filt_m9[,-c(1,2,13,14)]

for(i in 1:ncol(mpp_filt_m9_feat)){
    mpp_filt_m9_feat[[i]] <- as.numeric(as.character(mpp_filt_m9_feat[[i]]))
}

### Generate minimum values for each mzcat
mz_filt_m9 <- mpp_filt_m9[,c(1:2)]

for(j in 1:ncol(mpp_filt_m9_feat)){
    for(i in 1:nrow(mpp_filt_m9_feat)){
        if(mpp_filt_m9_feat[i,j] == 0 | is.null(mpp_filt_m9_feat[i,j])){
            mpp_filt_m9_feat[i,j] <- (0.5*min(mpp_filt_m9_feat[1, which(mpp_filt_m9_feat[1,] > 0)]))
        }
    }
}

## PANDA
mpp_m9_panda <- cbind(mz_filt_m9, mpp_filt_m9_feat)
lab <- c("mz", "time", "MPP_1", "Control_1","MPP_2", "Control_2", "MPP_3", "Control_3",  "MPP_4", "Control_4", "MPP_5", "Control_5")
label <-matrix(NA, nrow = 1, ncol = 12)
colnames(label) <- names(mpp_m9_panda)
mpp_m9_pan <- rbind(label, mpp_m9_panda)
mpp_m9_pan[1,] <- lab

setwd("C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\c18")

write.table(mpp_m9_pan,"mpp_m9_15_c18_panda.txt", col.names = F, row.names = F, sep = "\t")

## Metaboanalyst
mpp_m9_pan$mztime <- paste(mpp_m9_pan$mz, mpp_m9_pan$time, sep = "/")
mpp_m9_metab <- mpp_m9_pan[,c(13, 3:12)]
lab <- c("mztime", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control", "MPP", "Control")
mpp_m9_metab[1,] <- lab
write.table(mpp_m9_metab,"mpp_m9_15_c18_metab.txt", col.names = T, row.names = F, sep = "\t")

library(xmsPANDA)

feature_table_file<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\c18\\mpp_m9_15_c18_panda.txt"
class_labels_file<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\c18\\classlabels_m9.txt"
outloc<-"C:\\Users\\vkalia2\\Desktop\\Miller_lab\\Metabolomics\\Metabolomics\\Data_analysis\\cat1_murphy_mpp\\MPP\\20180815\\c18\\PANDA"

demetabs_res<-diffexp(feature_table_file=feature_table_file,
                      parentoutput_dir=outloc,
                      class_labels_file=class_labels_file,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="zeros",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=NA, group.missing.thresh=NA, input.intensity.scale="raw", 
                      log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE, 
                      quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
                      rsd.filt.list = c(0), pairedanalysis = FALSE, featselmethod="limma",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="pearson", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("royalblue1","darkgrey"), 
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2, 
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="BER",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=FALSE,max_comp_sel=1,saveRda=FALSE,pca.cex.val=4,pls.permut.count=NA,
                      pca.ellipse=TRUE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=5)


sink(file=NULL)
#####################################################
