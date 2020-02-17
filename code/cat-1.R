# cat-1 metabolomics

# Figures for cat-1 paper
# Heatmap
# MWAS
# Pathway analysis (metaboanalyst)
# box plots of top hits

library(tidyverse)

# Call in raw files
map <- read.table("Raw_files/Hilic_pos/Hilicpos_sample_id_mapfile.txt", header = T, sep = "\t")

feat <- read.table("Raw_files/Hilic_pos/Merged_pall_averaged_featuretable_repmax2of3.txt", header = T, sep = "\t")

## Extracting the cat1 features
## From sequence list, cat-1 samples were run together, with internal randomization
## Changing the map file
map.cat <- map[1:114,] %>%
    filter(str_detect(Sample.ID, "cat|N2")) %>%
    mutate(rep = rep(1:3, 32)) %>%
    filter(rep == 1) %>%
    mutate(Sample = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[1]])) %>%
    mutate(Rep = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[2]])) 

# To split string can also use the following. But will lose original column 
#separate(Sample.ID, c("Sample", "Replicate", "Injection"))

## Changing the feature table
#names(feat) <- gsub(pattern = ".mzXML", replacement = "", x = names(feat), fixed = TRUE) 
names.cat <- names(feat)[(names(feat) %in% map.cat$File.Name)]
mz.cat <- feat[,1:2]
cv.cat <- feat$median_CV
feat.cat <- feat %>%
    subset(select = names.cat) %>%
    bind_cols(mz.cat, .)

## Extracting data on M9 blank (ignoring bacterial blank)
map.m9 <- map.cat[grep("M9", map.cat$Sample),]
names.m9 <- names(feat.cat)[(names(feat.cat) %in% map.m9$File.Name)]
feat.m9 <- feat.cat %>%
    subset(select = names.m9)
feat.m9$m9.mean <- rowMeans(feat.m9)

## Extracting the worm samples
map.worms <- map.cat[!grepl("M9|Bact", map.cat$Sample),]
names.worms <- names(feat.cat)[(names(feat.cat) %in% map.worms$File.Name)]
feat.worms <- feat.cat %>%
    subset(select = names.worms)

cat.m9 <- cbind(feat.worms, feat.m9$m9.mean)

tcat.m9 <- as.data.frame(t(cat.m9))
#tcat.m9[1:5,1:5]

sum <- data.frame()

for(j in 1:12){
    for(i in 1:ncol(tcat.m9)){
        if(tcat.m9[j,i] > 1.5 * tcat.m9[13,i]){
            sum[j,i] <- 1
        }
        else{
            sum[j,i] <- 0
        }
    }
}

sum[1:5, 1:5]
dim(sum)
col.sum <- colSums(sum)
hist(col.sum)

names(col.sum) <- names(tcat.m9)
mz.cat <- as.data.frame(t(mz.cat))
names(mz.cat) <- names(tcat.m9)
cat.sum <- rbind(mz.cat,tcat.m9[1:12,], col.sum)
dim(cat.sum) #15 x 21479
cat.sum[1:5,1:5]
cat.filt <- cat.sum[-15, which(cat.sum[15,] >= 1)] 
dim(cat.filt) # 14 x 17782

## Labels
cat.filt <- as.data.frame(t(cat.filt))
mz.filt <- cat.filt[,1:2]
lab <- c("mz", "time", map.worms$Sample)
names(lab) <- names(cat.filt)
cat.metab <- rbind(lab, cat.filt)


cat.metab[1:5,1:5]

mz.time <- paste(cat.metab$mz[-1], cat.metab$time[-1],sep = "_")
#mz.time[1:5]

library(janitor)
cat.metab.comp <- cat.metab  %>%
    dplyr::select(-mz, -time) %>%
    t(.)  %>%
    as.data.frame() %>%
    rownames_to_column("sample")

names(cat.metab.comp) <- c("sample", "strain", mz.time)
cat.metab.comp[1:5,1:5]

# Convert factor to numeric
cat.comp <- as.data.frame(apply(cat.metab.comp[,3:ncol(cat.metab.comp)], 2, function(x) as.numeric(as.character(x))))
cat.sample.strain <- cat.metab.comp[,1:2]

### Data imputation
# Replacing zeros and null vales with min/2
# Function to impute with min / 2: https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r
replacezero <- function(x) "[<-"(x, !x|x==0, min(x[x>0], na.rm = T)/2)
cat.imp <- as.data.frame(apply(cat.comp[,3:ncol(cat.comp)], 2, replacezero)) #down a column
cat.imp[1:5,1:5]
cat.feat.imp <- as.data.frame(cbind(cat.sample.strain, cat.imp))
cat.feat.imp[1:5,1:5]

write.table(cat.feat.imp, "results/feat.cat1.box.txt", col.names = T, row.names = F, sep = "\t")
# PCA
feat.pca <- prcomp(cat.feat.imp[,3:ncol(cat.feat.imp)], center = T, scale = T)
plot(feat.pca)
#biplot(feat.pca)

library(ggbiplot)
ggbiplot(feat.pca,ellipse=TRUE,  groups=factor(cat.feat.imp$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Strain", values= c("red3", "grey40"))+
    theme_minimal()

# Quantile normalize and log transform, auto scale
library(preprocessCore)
#row.names(feat.70.imp) <- samples

# Convert to  matrix form: samples in columns, features in rows
feat.norm <- as.matrix(t(cat.feat.imp[,3:ncol(cat.feat.imp)]))
feat.norm[1:5,1:5]
#feat.norm <- apply(feat.norm, 2, function(x) as.numeric(x))
feat.normalized <- feat.norm %>%
    normalize.quantiles(.,copy=TRUE) %>%
    as.data.frame(.) 
#feat.normalized[1:5,1:5]
names(feat.normalized) <- cat.feat.imp$sample
row.names(feat.normalized) <- row.names(feat.norm)


#feat.normalized[1:5,1:5]

# Maybe plot how total ion intensities changed across samples?
# Change in top 20 features?
#feat.normalized[1:5,1:5]
#feat.norm[1:5,1:5]

sum.int <- as.data.frame(apply(feat.norm,1,sum))
sum.int$feat <- row.names(feat.norm)
names(sum.int)[1] <- "sum.int"
sum <- sum.int %>%
    arrange(desc(sum.int)) 

sum.int.sample <- as.data.frame(apply(feat.norm,2,sum))
names(sum.int.sample)[1] <- "sum.int"
sum.sample <- sum.int.sample %>%
    mutate(sample = c(1:12))

plot(sum.int~sample, data = sum.sample)

feat.norm.sub <- as.data.frame(t(feat.norm[rownames(feat.norm) %in% sum$feat[1:20],]))
feat.normalized.sub <- as.data.frame(t(feat.normalized[rownames(feat.normalized) %in% sum$feat[1:20],]))
feat.norm.sub$feat <- row.names(feat.norm.sub)
feat.norm.plot <- feat.norm.sub %>%
    gather(key = feat, value = intensity)%>%
    mutate(norm = "pre")

### log transformation
feat.normalized.log10 <- feat.normalized %>%
    log10(.) %>%
    t(.) %>%
    as.data.frame()      

feat.normalized.log10[1:5,1:5]

### Auto-scaling
#install.packages("RFmarkerDetector")
library(RFmarkerDetector)
?autoscale
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
feat.norm.log.scale$strain <- map.worms[order(map.worms$File.Name),]$Sample
feat.norm.log.scale[1:5,1:5]

# PCA with own centering and scaling
feat.pca <- prcomp(feat.norm.log.scale[,1:ncol(feat.norm.log.scale)-1], center = F, scale. = F)
library(ggbiplot)
ggbiplot(feat.pca,ellipse=TRUE,  groups=factor(feat.norm.log.scale$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Strain", values= c("red3", "grey40"))+
    theme_minimal()

### PLS-DA
library(mixOmics)

tiff("figures/pls-da.tiff", width = 4.5, height = 4.5, units = 'in', res = 300)

## PLS-DA across individual samples
metab <- feat.norm.log.scale[,c(1:ncol(feat.norm.log.scale)-1)]
strain <- as.factor(feat.norm.log.scale$strain)           

## PLS-DA function
plsda.strain <- plsda(metab, strain, ncomp = 5) 


## Plot PLS-DA

plotIndiv(plsda.strain, comp = c(1,2), ind.names = FALSE, legend = F, ellipse = TRUE,
          title = "",style = "graphics",
          size.title = rel(1.2), size.xlabel = rel(1),
          size.ylabel = rel(1), size.axis = rel(0.7), size.legend = rel(0.8),
          legend.title = "", legend.title.pch = "",
          legend.position = "right",
          #xlim = c(-42, 25), ylim = c(-70, 50),
          col.per.group =  c("navy", "lightblue4")
          )

text(x  = -15, y =40, labels = c("cat1"), adj = 1, col = "navy")
text(x  =20, y = 30, labels = c("N2"), adj = 1, col = "lightblue4")

dev.off()

# Look at mzs and rts that feed into components 1, 2, and 10

# Extract VIP scores from PLS-DA to get features of interest
status.vip <- as.data.frame(vip(plsda.strain))
status.vip.1210 <- status.vip[,c(1,2)] #rownames contain mz and rt information


# T-test
# Break up feature file by strain
feat.cat1 <- feat.norm.log.scale %>%
    filter(strain == "cat1")
#feat.norm.log.scale[1:5,1:5]
feat.n2 <- feat.norm.log.scale %>%
    filter(strain == "N2")
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
    outcome = colnames(feat.cat1)[i]
    
    x <- feat.cat1 %>%
        pull(i)
    
    y = feat.n2 %>%
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

#outcome$out_variable <- gsub("X", "",outcome$out_variable) #remove X from name of columns to retain mz_rt ids
outcome$mz <- unlist(lapply(strsplit(as.character(outcome$out_variable), "_"), function(x) as.numeric(x[1])))
outcome$logp <- -log10(outcome$mean_pvalue)
outcome$col <- ifelse(outcome$mean_pvalue < 0.05, 1, 2)
sum(outcome$col == 1) #114
redblack  <- c("red3", "black")


# Plotting t-test results
tiff("figures/mwas.tiff", width =4.5, height = 4, units = 'in', res = 300)

ggplot(data = outcome, aes(x = time, y = logp, color = factor(col))) +
    geom_point() +
    scale_color_manual(values=c("red3", "darkgrey")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) +
    geom_hline(yintercept = 1.30103, color = "blue") + #p = 0.05
    geom_hline(yintercept = 0.6020, color = "blue", linetype = "dashed") + #p = 0.25
    theme_classic() +
    ylab("-log10(p)") +
    xlab("retention time") +
    theme(legend.position = "none")

dev.off()

sum(col)
attach(outcome)
outcome.order <- outcome[order(-logp),]
detach(outcome)
class(outcome$mean_estimate)

outcome %>%
    dplyr::select(mz, mean_pvalue, mean_estimate) %>%
    mutate(m.z = mz, p.value = mean_pvalue, t.score = mean_estimate) %>%
    dplyr::select(m.z, p.value, t.score) %>%
    write.table("cat1.t.txt", sep = "\t", col.names = T, row.names = F)

outcome %>%
    dplyr::select(mz, time, mean_pvalue, mean_estimate) %>%
    mutate(m.z = mz, p.value = mean_pvalue, t.score = mean_estimate) %>%
    dplyr::select(m.z, time, p.value, t.score) %>%
    write.table("cat1.t.mcg.txt", sep = "\t", col.names = T, row.names = F)

# Pathway analysis
## Bubble plot
mem <- read.table("results/mcg/mummichog_pathway_enrichment.csv",
                  header = T, sep = ",")

mem$logp <- -log10(mem$FET)
mem$enrich <- mem$Hits.sig/mem$Expected
#eth.hispcauc <- merge(hisp.sub, cauc.sub, by = "pathway", all = T)    
#eth.all <- merge(eth.hispcauc, afr.sub, by = "pathway", all = T)
mem.sub <- mem %>%
    filter(logp > 0.5) %>%
    filter(enrich > 0.5)

#tiff("figures/pathways.tiff", width = 5.5, height = 4, units = 'in', res = 300)
ggplot(mem.sub, aes(x=logp, y = reorder(X, logp), size = Hits.total,  col = enrich)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low="blue", high="red")+
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways altered in cat-1 worms", 
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

### No overlap size on graph, enrichment as size of bubble
mem.sub$label <- paste0(mem.sub$X," (",mem.sub$Hits.sig,"/", mem.sub$Pathway.total,")")
ggplot(mem.sub, aes(x=logp, y = reorder(label, logp), size = enrich)) +
    geom_point() +
    #scale_color_gradient(low="blue", high="red") +
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways altered in cat-1 worms", 
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



# Heatmap
# Select p-value < 0.05 from regression
outcome.mwas <- outcome %>%
    filter(mean_pvalue<0.05)

outcome.mwas.mz <- as.data.frame(outcome.mwas[,1])

names(outcome.mwas.mz) <- "mz_time"

metab.mwas <- metab[,which(names(metab) %in% outcome.mwas.mz$mz_time)]
feat.norm.log.scale$strain
heatmap.rows <- paste0(feat.norm.log.scale$strain, ".", c("1", "1", "2", "2", "3", "3", "4", "5", "4", "5", "6", "6"))
heatmap_rows <- heatmap.rows

row.names(metab.mwas) <- heatmap_rows

tiff("figures/heatmap_mwas_p005.tiff", width = 9, height = 7, units = 'in', res = 300)
library(RColorBrewer)
my_group <- as.numeric(as.factor(substr(as.character(row.names(metab.mwas)),1, 2)))
colSide <- c("navy", "lightblue4")[my_group]
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

# For metaboanalyst
cat.metab.1.5 <- cat.metab %>%
    mutate(mztime = paste(mz, time, sep = "/")) %>%
    select(mztime, everything()) %>%
    select(-mz, -time)

write.table(cat.metab.1.5, "/Users/vk2316/Documents/Miller_lab/Metabolomics_data/MPP_Murphy_CAT1/Exploratory_analysis/Metaboanalyst/cat_1/cat1_metaboanalyst.txt", col.names = T, row.names = F, sep = "\t")  

# PCA with all sample types
map.all <- map[1:114,] %>%
    mutate(rep = rep(1:3, 38)) %>%
    filter(rep == 1) %>%
    mutate(Sample = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[1]])) %>%
    mutate(Rep = sapply(strsplit(as.character(Sample.ID), "_"), function(i) i[[2]])) 


# PCA with own centering and scaling
feat.pca <- prcomp(feat.norm.log.scale[,1:ncol(feat.norm.log.scale)-1], center = F, scale. = F)
library(ggbiplot)
ggbiplot(feat.pca,ellipse=TRUE,  groups=factor(feat.norm.log.scale$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Strain", values= c("red3", "grey40"))+
    theme_minimal()

names.all <- names(feat)[(names(feat) %in% map.all$File.Name)]
mz.all <- feat[,1:2]
cv.all <- feat$median_CV
feat.all <- feat %>%
    subset(select = names.all) %>%
    bind_cols(mz.all, .)

# Extract features with intensity 1.5 times that of M9
feat.all.filt <- feat.all[feat.all$mz %in% cat.filt$mz,]
feat.all.filt[1:5,1:5]
feat.all.filt.mz.time <- feat.all.filt[,1:2]
feat.all.filt.mztime <- paste0(feat.all.filt$mz, "_", feat.all.filt$time)

# Transpose and make amenable for computation
library(janitor)
feat.comp <- feat.all.filt  %>%
    dplyr::select(-mz, -time) %>%
    t(.)  %>%
    as.data.frame() %>%
    rownames_to_column("sample")
feat.comp[1:5,1:5]

names(feat.comp) <- c("sample", feat.all.filt.mztime)
feat.comp[1:5,1:5]
dim(feat.comp) #8 x 4258
sample <- feat.comp[,1]

# Convert factor to numeric
feat.all.comp <- as.data.frame(apply(feat.comp[,2:ncol(feat.comp)], 2, function(x) as.numeric(as.character(x))))

### Data imputation
# Replacing zeros and null vales with min/2
# Function to impute with min / 2: https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r
replacezero <- function(x) "[<-"(x, !x|x==0, min(x[x>0], na.rm = T)/2)
feat.imp <- as.data.frame(apply(feat.all.comp, 2, replacezero)) #down a column
feat.imp[1:5,1:5]
feat.all.imp <- as.data.frame(cbind(sample, feat.imp))
feat.all.imp[1:5,1:5]

# Adding strain information
feat.all.imp$strain <- map.all[order(map.all$File.Name),]$Sample
feat.all.imp[1:5,1:5]
class(feat.all.imp$`85.0476951464037_79.4874455703363`)

# PCA
#feat.all <- as.data.frame(apply(feat.all.imp[,2:ncol(feat.all.imp)], 2, function(x) as.numeric(as.character(x))))
feat.pca <- prcomp(feat.all.imp[,2:4257], center = T, scale = T)
plot(feat.pca)

library(ggbiplot)
ggbiplot(feat.pca, ellipse=F,  groups=factor(feat.all.imp$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Type", values= c("red3", "grey40", "green", "blue", "brown", "yellow", "purple"))+
    theme_minimal()


### Box plots
# Pathway related? trytophan? tyrosine?
# mcg annotation? 
# Top 5 features from MWAS?



#################################################
##               DEVELOPMENT                   ##
#################################################
#### Change thresholding criteria??
cat.filt <- cat.sum[-15, which(cat.sum[15,] >= 0.3 *12)] 
dim(cat.filt) # 14 x 4257
