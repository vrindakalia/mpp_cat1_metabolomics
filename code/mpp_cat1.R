###############################################
##         MPP+ compared to cat-1            ##
###############################################
# Two ways to do this: metabnet and mummichog on overlapping features
# xMSanalyzer to get overlapping features

library(MetabNet)
library(xMSanalyzer)
library(tidyverse)
library(janitor)

cat1_mztime <- read.table("results/mcg/cat1.t.mcg.txt", header =  T, sep = "\t") %>%
    arrange(m.z)

mpp_mztime <- read.table("results/mcg/mpp.t.mcg.txt", header =  T, sep = "\t") %>%
    arrange(m.z)

#getVenn(cat1_mztime, "CAT-1", mpp_mztime , "MPP+", mz.thresh = 10, time.thresh=30,
#        alignment.tool = NA, 
#        xMSanalyzer.outloc = "/Users/vk2316/Documents/Miller_lab/Metabolomics_data/mpp_cat1_metabolomics/comparisons",
#        use.unique.mz=F, plotvenn=TRUE)

#overlap.mz <- find.Overlapping.mzs(cat1_mztime, mpp_mztime, mz.thresh = 10, time.thresh = NA,
#                                   alignment.tool=NA)

# Using features that meet p < 0.05 criteria to compare overlapping features
cat1_top <- cat1_mztime %>%
    filter(p.value < 0.25)

mpp_top <- mpp_mztime %>%
    filter(p.value < 0.2)

overlap_005 <- getVenn(cat1_top, "cat-1 \n(p < 0.25)", mpp_top, "MPP+ (p < 0.2)", mz.thresh = 10, time.thresh=30,
                       alignment.tool = NA, 
                       xMSanalyzer.outloc = "/Users/vk2316/Documents/Miller_lab/Metabolomics_data/mpp_cat1_metabolomics/results/comparisons/overlap",
                       use.unique.mz=F, plotvenn=TRUE)

overlap_mztime <- overlap_005$common
overlap_mztime$mz_time <- paste0(overlap_mztime$mz.data.A,"_", overlap_mztime$time.data.A)

##### Creating a file that lets us see the pathways that these overlapping features belong to: 
# give these features a cut-off p-value, change p-value on others

# features not in the overlap
cat1_mztime$mz_time <- paste0(cat1_mztime$m.z,"_", cat1_mztime$time)
overlap_mztime <- cat1_mztime[which(cat1_mztime$mz_time %in% overlap_mztime$mz_time),]
nonoverlap_mztime <- cat1_mztime[-which(cat1_mztime$mz_time %in% overlap_mztime$mz_time),]

overlap_mztime$p.dum <- rep(0.01, times = dim(overlap_mztime)[1])

nonoverlap_mztime$p.dum <- rep(0.5, times = dim(nonoverlap_mztime)[1])

dum.mcg <- rbind(overlap_mztime, nonoverlap_mztime) %>%
    dplyr::select(m.z, p.dum, t.score) %>%
    mutate(p.value = p.dum) %>%
    dplyr::select(m.z, p.value, t.score) %>%
    write.table("results/comparisons/mcg_hilicpos_cat1mppcompare.txt", col.names = T, sep = "\t", row.names = F)

overlap.feat <- overlap_mztime[,1:2] %>%
    write.table("results/comparisons/metabnet/sig.metab.file.txt", col.names = T, row.names = F, sep = "\t")
mz_time  <- (overlap_mztime[,1:2])

#######################################################
# Look for features of interest in MPP and CAT-1 FT  ##
#######################################################

compounds <- read.table("results/comparisons/mummichog_matched_compound_all_compare.csv", header = T, sep = ",")
pathways <- read.table("results/comparisons/mummichog_pathway_enrichment_compare.csv", header = T, sep = ",")

pathways$logp <- -log10(pathways$FET)
pathways$enrich <- pathways$Hits.sig/pathways$Expected
#eth.hispcauc <- merge(hisp.sub, cauc.sub, by = "pathway", all = T)    
#eth.all <- merge(eth.hispcauc, afr.sub, by = "pathway", all = T)
mem.sub <- pathways %>%
    filter(logp > 0.35)

#tiff("figures/pathways.comparison.tiff", width = 5.7, height = 4, units = 'in', res = 300)
ggplot(mem.sub, aes(x=logp, y = reorder(X, logp), size = Hits.total,  col = enrich)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low="blue", high="red")+
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways enriched by overlapping features", 
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

tiff("figures/pathways.comparison.black.tiff", width = 6.5, height = 4, units = 'in', res = 300)

### No overlap size on graph, enrichment as size of bubble
mem.sub$label <- paste0(mem.sub$X," (",mem.sub$Hits.sig,"/", mem.sub$Pathway.total,")")
ggplot(mem.sub, aes(x=logp, y = reorder(label, logp), size = enrich)) +
    geom_point() +
    #scale_color_gradient(low="blue", high="red") +
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    ggtitle("Pathways enriched by overlapping features", 
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

# From  mummichog metaboanalyst interface, compounds that are part of pathways and in overlap:
#1. Glycero-3-phosphocholine C00670
#2. Choline C00114
#3. N-Acetyl-D-glucosamine-6-phosphate  C00357
#4. D-glucosamine C00329
#5. Malate C00149
#6. 4-hydroxyphenylacetate C00642
#7. Carnitine C00487

comp.interest <- c("C00670", "C00114", "C00357", "C00329", "C00149", "C00642", "C00487")

info.interest <- compounds[compounds$Matched.Compound %in% comp.interest,] %>%
    mutate(Name = ifelse(Matched.Compound == comp.interest[1], "Glycero-3-phosphocholine",
                         ifelse(Matched.Compound == comp.interest[2], "Choline",
                                ifelse(Matched.Compound == comp.interest[3], "N-Acetyl-D-glucosamine-6-phosphate",
                                       ifelse(Matched.Compound == comp.interest[4], "D-glucosamine",
                                              ifelse(Matched.Compound == comp.interest[5], "Malate",
                                                     ifelse(Matched.Compound == comp.interest[6], "4-hydroxyphenylacetate", "Carnitine")))))))

info.interest #contains: mz, kegg id, matched form, mass.diff, Name
info.box <- info.interest %>%
    group_by(Name, Matched.Form) %>%
    mutate(min.adduct = min(Mass.Diff)) %>%
    ungroup() %>%
    filter(Mass.Diff == min.adduct)

write.table(info.box, "results/info.box.adduct.table.txt",  col.names = T, row.names = F, sep =  "\t")
# Next: subset the feature table to find these mzs and create box plots to see changes in levels between 
# cat1 and N2 and 
# MPP exposure v control
#Call in feature table with intensities



# Find mzs also  present in overlap file
ind <- vector(mode = "logical", length = nrow(info.box))
for(i in 1: nrow(info.box)) {
    ind[i] <- which(mz_time$m.z %in% info.box.1$Query.Mass[i])
    print(i)
}

ind


feat.mpp <- read.table("results/feat.mpp.box.txt", header = T, sep = "\t")
feat.cat1 <- read.table("results/feat.cat1.box.txt", header = T, sep = "\t")

mz.mpp <- unlist(lapply(strsplit(as.character(names(feat.mpp)), "_"), function(x) x[1]))
time.mpp <- unlist(lapply(strsplit(as.character(names(feat.mpp)), "_"), function(x) x[2]))
mz.mpp <- gsub("X", "", mz.mpp)

mz.cat1 <- unlist(lapply(strsplit(as.character(names(feat.cat1)), "_"), function(x) x[1]))
time.cat1 <- unlist(lapply(strsplit(as.character(names(feat.cat1)), "_"), function(x) x[2]))
mz.cat1 <- gsub("X", "", mz.cat1)

#cbind(info.box[1,],feat.cat1[,c(2,which(mz.cat1 %in% info.box$Query.Mass[1]))])

#feat.mpp[,c(2,which(mz.mpp %in% info.box$Query.Mass[1]))]

#box.mpp <- cbind(info.box[1,], feat.mpp[,c(2,which(mz.mpp %in% info.box$Query.Mass[1]))])

data.plots <- data.frame()

for(i in 1:7){
    # CAT-1
    box <- feat.cat1[,c(2,which(mz.cat1 %in% info.box$Query.Mass[i]))]
    box.plot <- box[,1:2]
    names(box.plot) <- c("Strain", "Intensity")
    box.plot$Strain <- as.factor(box.plot$Strain)
    #box.plot$Intensity <- as.numeric(levels(box.plot$Intensity))[box.plot$Intensity]
    box.plot$log.intensity <- log2(box.plot$Intensity)
    box.plot$neuro <- rep("cat1", 12)
    # MPP
    box.mpp <- feat.mpp[,c(2,which(mz.mpp %in% info.box$Query.Mass[i]))]
    box.plot.mpp <- box.mpp[,1:2]
    names(box.plot.mpp) <- c("Strain", "Intensity")
    box.plot.mpp$Strain <- as.factor(box.plot.mpp$Strain)
    #box.plot$Intensity <- as.numeric(levels(box.plot$Intensity))[box.plot$Intensity]
    box.plot.mpp$log.intensity <- log2(box.plot.mpp$Intensity)
    box.plot.mpp$neuro <- rep("MPP", 10)
    box.plot.all <- rbind(box.plot, box.plot.mpp)
    box.plot.all$Strain.label <- case_when(
        box.plot.all$Strain == "C" ~ "Control",
        box.plot.all$Strain == "MPP" ~ "MPP+",
        box.plot.all$Strain == "cat1" ~ "cat-1",
        box.plot.all$Strain == "N2" ~ "N2"
    )
    box.plot.all$Strain.label <- factor(box.plot.all$Strain.label, levels = c("cat-1", "N2", "MPP+", "Control"))
    name <- info.box$Name[i]
    mz <- info.box$Query.Mass[i]
    adduct <- info.box$Matched.Form[i]
    boxes <-  ggplot(data = box.plot.all, aes(x = Strain.label,  y = log.intensity, fill = Strain.label)) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_manual(values = c("royalblue3", "lightblue4", "orange", "darkgray")) +
        ggtitle(name, subtitle=paste0("mz = ",  mz, "", "  Adduct = ", adduct)) +
        theme_minimal() +
        xlab("") +
        ylab("Log Intensity")
    tiff(paste0("results/comparisons/figures/comparison", name,".tiff"), width =4.5, height = 4, units = 'in', res = 300)
    print(boxes)
    dev.off()
    dat.meta <- cbind(box.plot.all, name, mz, adduct)
    data.plots <- rbind(data.plots, dat.meta)
    #boxplot(log.intensity~Strain, data = box.plot.all, ylab = "log2(Intensity)", xlab="Strain", 
     #       cex.lab = 1, cex.axis = 1, col = c('gray', "white"))
    #title(main = paste(strwrap(substring(box$Name[1], 1, 100), width = 60), collapse = "\n"), cex.main = 1)
    #text(x = 2.25, y = max(box.plot$log.intensity), labels = "mz=", cex =0.75)
    #text(x = 2.4, y = max(box.plot$log.intensity), labels = round(as.numeric(box$Query.Mass),4), cex =0.75)
    #text(x = 2.4, y = max(box.plot$log.intensity), labels = box$Matched.Form, pos = 1,cex =0.75)
    print(i)
    }

write.table(data.plots, "results/comparisons/data_comparison_plots.txt", col.names = T,  row.names = F, sep = "\t")

##############
# Remove line 9 from info.box since not present in MPP+ feature table
#info.box <- info.box[-9,]
# Remove line 10 from info.box since not present in MPP+ feature table
#info.box <- info.box[-10,]

box.plot.all$Strain.label <- case_when(
    box.plot.all$Strain == "C" ~ "Control",
    box.plot.all$Strain == "MPP" ~ "MPP+",
    box.plot.all$Strain == "cat1" ~ "cat-1",
    box.plot.all$Strain == "N2" ~ "N2"
)
