###############################################
##       MPP+ treatment in worms             ##
###############################################

library(tidyverse)
library(janitor)
library(mixOmics)
library(ggbiplot)
library(gplots)
library(preprocessCore)
library(RFmarkerDetector)

# Call in raw files
map <- read.table("Raw_files/Hilic_pos/Hilicpos_sample_id_mapfile.txt", header = T, sep = "\t")

feat <- read.table("Raw_files/Hilic_pos/Merged_pall_averaged_featuretable_repmax2of3.txt", header = T, sep = "\t")

## Extracting the MPP+ features
## From sequence list, mpp samples were run together, with internal randomization
## Changing the map file
map.mpp <- map[157:234,] %>%
    filter(str_detect(Sample.ID, "MPP|C|M9_B")) %>% #select eexposure, control, and M9 blank  samples
    mutate(rep = rep(1:3, 15)) %>% # add number to indicate injection number
    filter(rep == 1) %>% #extract only the first injection to  match the names  in feature table
    mutate(Sample = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[2]])) %>% # create sample type variable
    mutate(Rep = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[3]])) # create variable to decide rep

## Extract mpp samples from the feature table
#names(feat) <- gsub(pattern = ".mzXML", replacement = "", x = names(feat), fixed = TRUE) 
names.mpp <- names(feat)[(names(feat) %in% map.mpp$File.Name)]
mz.mpp <- feat[,1:2]
feat.mpp <- feat %>%
    subset(select = names.mpp) %>%
    bind_cols(mz.mpp, .)

## Extracting data on M9 blank to use for filtering features
mpp.m9 <- map.mpp[grep("B", map.mpp$Sample),]
names.mpp.m9 <- names(feat.mpp)[(names(feat.mpp) %in% mpp.m9$File.Name)]
feat.mpp.m9 <- feat.mpp %>%
    subset(select = names.mpp.m9)
feat.mpp.m9$m9.mean <- rowMeans(feat.mpp.m9)

## Extracting the worm samples
mpp.worms <- map.mpp[!grepl("B", map.mpp$Sample),]
names.mpp.worms <- names(feat.mpp)[(names(feat.mpp) %in% mpp.worms$File.Name)]
feat.mpp.worms <- feat.mpp %>%
    subset(select = names.mpp.worms)

mpp.m9 <- cbind(feat.mpp.worms, feat.mpp.m9$m9.mean)

tmpp.m9 <- as.data.frame(t(mpp.m9))

tcat.m9[1:5,1:5]

sum.mpp <- data.frame()

for(j in 1:10){
    for(i in 1:ncol(tcat.m9)){
        if(tmpp.m9[j,i] > 1.5 * tmpp.m9[11,i]){
            sum.mpp[j,i] <- 1
        }
        else{
            sum.mpp[j,i] <- 0
        }
    }
}

sum.mpp[1:5, 1:5]
dim(sum.mpp)
col.sum.mpp <- colSums(sum.mpp)
hist(col.sum.mpp)

names(col.sum.mpp) <- names(tmpp.m9)

mz.mpp <- as.data.frame(t(mz.mpp))
names(mz.mpp) <- names(tmpp.m9)
mpp.sum <- rbind(mz.mpp, tmpp.m9[1:10,], col.sum.mpp)

dim(mpp.sum) #13 x 21479

mpp.sum[1:5,1:5]
mpp.filt <- mpp.sum[-13, which(mpp.sum[13,] >= 1)] 

dim(mpp.filt) # 12 x 14734

## Labels
mpp.filt <- as.data.frame(t(mpp.filt))
mz.filt <- mpp.filt[,1:2]
lab <- c("mz", "time", mpp.worms$Sample)
names(lab) <- names(mpp.filt)
mpp.metab <- rbind(lab, mpp.filt)

mpp.metab[1:5,1:5]

mz.time <- paste(mpp.metab$mz[-1], mpp.metab$time[-1],sep = "_")
#mz.time[1:5]

mpp.metab.comp <- mpp.metab  %>%
    dplyr::select(-mz, -time) %>%
    t(.)  %>%
    as.data.frame() %>%
    rownames_to_column("sample")

names(mpp.metab.comp) <- c("sample", "strain", mz.time)
mpp.metab.comp[1:5,1:5]

# Convert factor to numeric
mpp.comp <- as.data.frame(apply(mpp.metab.comp[,3:ncol(mpp.metab.comp)], 2, function(x) as.numeric(as.character(x))))
mpp.sample.strain <- mpp.metab.comp[,1:2]

### Data imputation
# Replacing zeros and null vales with min/2
# Function to impute with min / 2: https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r
replacezero <- function(x) "[<-"(x, !x|x==0, min(x[x>0], na.rm = T)/2)
mpp.imp <- as.data.frame(apply(mpp.comp[,3:ncol(mpp.comp)], 2, replacezero)) #down a column
mpp.imp[1:5,1:5]
mpp.feat.imp <- as.data.frame(cbind(mpp.sample.strain, mpp.imp))
mpp.feat.imp[1:10,1:5]

write.table(mpp.feat.imp, "results/feat.mpp.box.txt", col.names = T, row.names = F, sep = "\t")

# PCA
feat.pca.mpp <- prcomp(mpp.feat.imp[,3:ncol(mpp.feat.imp)], center = T, scale = T)
plot(feat.pca.mpp)
summary(feat.pca.mpp)


ggbiplot(feat.pca.mpp,ellipse=TRUE,  groups=factor(mpp.feat.imp$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Treatment", values= c("red3", "grey40"))+
    theme_minimal()

# Quantile normalize and log transform, auto scale
# Convert to  matrix form: samples in columns, features in rows
feat.norm.mpp <- as.matrix(t(mpp.feat.imp[,3:ncol(mpp.feat.imp)]))
feat.norm.mpp[1:5,1:5]
#feat.norm <- apply(feat.norm, 2, function(x) as.numeric(x))
feat.normalized.mpp <- feat.norm.mpp %>%
    normalize.quantiles(.,copy=TRUE) %>%
    as.data.frame(.) 
feat.normalized.mpp[1:5,1:5]
names(feat.normalized.mpp) <- mpp.feat.imp$sample
row.names(feat.normalized.mpp) <- row.names(feat.norm.mpp)

#feat.normalized.mpp[1:5,1:5]

# Maybe plot how total ion intensities changed across samples?
# Change in top 20 features?
#feat.normalized[1:5,1:5]
#feat.norm[1:5,1:5]

sum.int.mpp <- as.data.frame(apply(feat.norm.mpp,1,sum))
sum.int.mpp$feat <- row.names(feat.norm.mpp)
names(sum.int.mpp)[1] <- "sum.int"
sum.mpp <- sum.int.mpp %>%
    arrange(desc(sum.int)) 

sum.int.sample <- as.data.frame(apply(feat.norm.mpp,2,sum))
names(sum.int.sample)[1] <- "sum.int"
sum.sample.mpp <- sum.int.sample %>%
    mutate(sample = c(1:10))

plot(sum.int~sample, data = sum.sample.mpp)

feat.norm.sub <- as.data.frame(t(feat.norm.mpp[rownames(feat.norm.mpp) %in% sum.mpp$feat[1:20],]))
feat.normalized.sub <- as.data.frame(t(feat.normalized.mpp[rownames(feat.normalized.mpp) %in% sum.mpp$feat[1:20],]))
feat.norm.sub$feat <- row.names(feat.norm.sub)
feat.norm.plot <- feat.norm.sub %>%
    gather(key = feat, value = intensity)%>%
    mutate(norm = "pre")

### log transformation
feat.normalized.log10 <- feat.normalized.mpp %>%
    log10(.) %>%
    t(.) %>%
    as.data.frame()      

feat.normalized.log10[1:5,1:5]

### Auto-scaling
#install.packages("RFmarkerDetector")
feat.norm.log.scale <- feat.normalized.log10 %>%
    autoscale(., exclude = F) %>%
    as.data.frame()

feat.norm.log.scale[1:10,1:5]

sum.int.normalized <- as.data.frame(apply(feat.norm.log.scale,2,sum))
sum.int.normalized$feat <- names(feat.norm.log.scale)
names(sum.int.normalized)[1] <- "sum.int"
sum.post <- sum.int.normalized %>%
    arrange(desc(sum.int)) 

feat.norm.log.scale.sub <- as.data.frame(feat.norm.log.scale[,names(feat.norm.log.scale) %in% sum$feat[1:20]])

feat.normalized.plot <- feat.norm.log.scale.sub %>%
    gather(key = feat, value = intensity)%>%
    mutate(norm = "post")

feat.both <- rbind(feat.norm.plot, feat.normalized.plot)
feat.both$norm2 <- factor(feat.both$norm,levels = c("pre", "post"))
ggplot(data = feat.both, aes(x=feat, y  = intensity)) +
    geom_boxplot() +
    facet_wrap(~norm2, scale = "free_x") +
    coord_flip() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

# Adding strain information
feat.norm.log.scale$strain <- mpp.worms[order(mpp.worms$File.Name),]$Sample
feat.norm.log.scale[1:5,1:5]

# PCA with own centering and scaling
feat.pca.own <- prcomp(feat.norm.log.scale[,1:ncol(feat.norm.log.scale)-1], center = F, scale. = F)
library(ggbiplot)
ggbiplot(feat.pca.own,ellipse=TRUE,  groups=factor(feat.norm.log.scale$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Strain", values= c("red3", "grey40"))+
    theme_minimal()

### PLS-DA
library(mixOmics)

tiff("figures/pls-da.mpp.tiff", width = 4.5, height = 4.5, units = 'in', res = 300)

## PLS-DA across individual samples
metab <- feat.norm.log.scale[,c(1:ncol(feat.norm.log.scale)-1)]
strain <- as.factor(feat.norm.log.scale$strain)           

## PLS-DA function
plsda.strain <- plsda(metab, strain, ncomp = 5) 


## Plot PLS-DA

plotIndiv(plsda.strain, comp = c(1,2), ind.names = FALSE, legend =F, ellipse = TRUE,
          title = "",style = "graphics",
          size.title = rel(1.2), size.xlabel = rel(1),
          size.ylabel = rel(1), size.axis = rel(0.7), size.legend = rel(0.8),
          legend.title = "", legend.title.pch = "",
          legend.position = "right",
          #xlim = c(-42, 25), ylim = c(-70, 50),
          col.per.group =  c("darkgray", "orange")
)

text(x  = 75, y = -60, labels = c("Control"), adj = 1, col = "darkgray")
text(x  = -60, y = 100, labels = c("MPP+"), adj = 1, col = "orange")

dev.off()

# Look at mzs and rts that feed into components 1, 2, and 10

# Extract VIP scores from PLS-DA to get features of interest
status.vip <- as.data.frame(vip(plsda.strain))
status.vip.1210 <- status.vip[,c(1,2)] #rownames contain mz and rt information


# T-test
# Break up feature file by strain
feat.mpp <- feat.norm.log.scale %>%
    filter(strain == "MPP")
#feat.norm.log.scale[1:5,1:5]
feat.control <- feat.norm.log.scale %>%
    filter(strain == "C")
#feat.norm.log.scale[1:5,1:5]

# OUTCOME
out_start= 1 # position of first mz
out_end= ncol(feat.norm.log.scale)-1 # position of last mz
out_nvar=out_end-out_start+1 # total number of dependent variables

out_variable=rep(NA, out_nvar) # matrix to hold all outcome data

# data from t-test
mean_estimate = rep(NA, out_nvar) # coeeficient for status
mean_se = rep(NA, out_nvar)
mean_pvalue=rep(NA, out_nvar)

number=1 # counter

for (i in out_start:out_end){
    outcome = colnames(feat.mpp)[i]
    
    x <- feat.mpp %>%
        pull(i)
    
    y = feat.control %>%
        pull(i) 
    
    test <- t.test(x, y, paired = F, alternative = "two.sided")
    
    beta <- test$statistic
    se <- test$stderr
    zval <- beta / se
    pval <- test$p.value
    
    mean_estimate[number] = as.numeric(beta)
    mean_se[number] = as.numeric(se)
    mean_pvalue[number] = as.numeric(pval)
    
    out_variable[number] = outcome
    
    #exp_beta[number] = as.numeric(beta[2])
    #exp_se[number] = as.numeric(se[2])
    #exp_pvalue[number] = as.numeric(pval[2])
    # exp_variable[number] = exposure
    number = number + 1
}


# Creating data frame with relevant regression output
outcome = data.frame(out_variable, mean_estimate, mean_se, mean_pvalue)

outcome$time <- unlist(lapply(strsplit(as.character(outcome$out_variable), "_"), function(x) as.numeric(x[2])))

outcome$out_variable <- gsub("X", "",outcome$out_variable) #remove X from name of columns to retain mz_rt ids
outcome$mz <- unlist(lapply(strsplit(as.character(outcome$out_variable), "_"), function(x) as.numeric(x[1])))

outcome$logp <- -log10(outcome$mean_pvalue)
outcome$col <- ifelse(outcome$mean_pvalue < 0.05, 1, 2)
sum(outcome$col == 1) #199
redblack  <- c("red3", "black")


# Plotting t-test results
tiff("figures/mwas.mpp.tiff", width =4.5, height = 4, units = 'in', res = 300)

ggplot(data = outcome, aes(x = time, y = logp, color = factor(col))) +
    geom_point() +
    scale_color_manual(values=c("red3", "darkgrey")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) +
    geom_hline(yintercept = 1.30103, color = "blue") + #p = 0.05
    geom_hline(yintercept = 1, color = "blue", linetype = "dashed") + #p = 0.1
    theme_classic() +
    ylab("-log10(p)") +
    xlab("retention time") +
    theme(legend.position = "none")

dev.off()

attach(outcome)
outcome.order <- outcome[order(-logp),]
detach(outcome)
class(outcome$mean_estimate)

outcome %>%
    dplyr::select(mz, mean_pvalue, mean_estimate) %>%
    mutate(m.z = mz, p.value = mean_pvalue, t.score = mean_estimate) %>%
    dplyr::select(m.z, p.value, t.score) %>%
    write.table("mpp.t.txt", sep = "\t", col.names = T, row.names = F)

outcome %>%
    dplyr::select(mz, time, mean_pvalue, mean_estimate) %>%
    mutate(m.z = mz, p.value = mean_pvalue, t.score = mean_estimate) %>%
    dplyr::select(m.z, time, p.value, t.score) %>%
    write.table("mpp.t.mcg.txt", sep = "\t", col.names = T, row.names = F)

# Pathway analysis
## Bubble plot
mem <- read.table("results/mcg/mummichog_pathway_enrichment_mpp.csv",
                  header = T, sep = ",")

mem$logp <- -log10(mem$FET)
mem$enrich <- mem$Hits.sig/mem$Expected
#eth.hispcauc <- merge(hisp.sub, cauc.sub, by = "pathway", all = T)    
#eth.all <- merge(eth.hispcauc, afr.sub, by = "pathway", all = T)
mem.sub <- mem %>%
    filter(logp > 0.5) %>%
    filter(enrich > 0.3)

#tiff("figures/pathways.mpp.tiff", width = 5.7, height = 4, units = 'in', res = 300)
ggplot(mem.sub, aes(x=logp, y = reorder(X, logp), size = Hits.total,  col = enrich)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low="blue", high="red")+
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways altered in worms exposed to MPP+", 
            subtitle = "Size of bubble represents number of significant hits \nEnrichment is calculated as (Total Hits/Pathway size)") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7),
          axis.text=element_text(size=7), 
          axis.title=element_text(size=9,face="bold"),
          strip.text = element_text(size=7),
          legend.text=element_text(size=7),
          legend.title=element_text(size=8),
          legend.position="bottom") +
    guides(size=guide_legend("Overlap size")) +
    labs(col = "Enrichment")
#dev.off()

tiff("figures/pathways.mpp.black.tiff", width = 6.7, height = 4, units = 'in', res = 300)

### No overlap size on graph, enrichment as size of bubble
mem.sub$label <- paste0(mem.sub$X," (",mem.sub$Hits.sig,"/", mem.sub$Pathway.total,")")
ggplot(mem.sub, aes(x=logp, y = reorder(label, logp), size = enrich)) +
    geom_point() +
    #scale_color_gradient(low="blue", high="red") +
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways altered in worms exposed to MPP+", 
            subtitle = "Enrichment is calculated as (Total Hits/Expected number of hits)") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7),
          axis.text=element_text(size=9), 
          axis.title=element_text(size=9,face="bold"),
          strip.text = element_text(size=7),
          legend.text=element_text(size=7),
          legend.title=element_text(size=8),
          legend.position="bottom") +
    labs(size = "Enrichment")

dev.off()

# Heatmap
# Select p-value < 0.05 from regression
outcome.mwas <- outcome %>%
    filter(mean_pvalue < 0.05)

outcome.mwas.mz <- as.data.frame(outcome.mwas[,1])

names(outcome.mwas.mz) <- "mz_time"

metab.mwas <- metab[,which(names(metab) %in% outcome.mwas.mz$mz_time)]
feat.norm.log.scale$strain
heatmap.rows <- paste0(feat.norm.log.scale$strain, ".", c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5"))
heatmap_rows <- heatmap.rows

row.names(metab.mwas) <- heatmap_rows

tiff("figures/heatmap_mwas_p005.mpp.tiff", width = 9, height = 7, units = 'in', res = 300)
library(RColorBrewer)
my_group <- as.numeric(as.factor(substr(as.character(row.names(metab.mwas)),1, 2)))
colSide <- c("darkgray", "orange")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "RdYlBu"))(7)
col.hard <- c("#4575B4","#4575B4" ,"#4575B4","black","#D73027","#D73027","#D73027")
names(metab.mwas)<-c(NA)
#distCor <- function(x) as.dist(1-cor(t(x)))
#hclustAvg <- function(x) hclust(x, method="average")
library(gplots)

heatmap.2(as.matrix(metab.mwas),
          Rowv=T,
          Colv=TRUE,
          trace='none',        # turns off trace lines inside the heat map
          density.info="none",
          #distfun=distCor, 
          #hclustfun=hclustAvg,
          RowSideColors = colSide, margins = c(5,10),
          xlab = "", ylab =  "",
          col = col.hard)   

dev.off()



