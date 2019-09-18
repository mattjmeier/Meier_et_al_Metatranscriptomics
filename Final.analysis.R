#### Script to reproduce analysis and metatranscriptomics figures from Meier et al.

library(biomformat)
library(edgeR)
library(DESeq2)
library(phyloseq)
library(scales)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# BIOM format from MG-RAST needs to have "type" changed to a compatible term such as "OTU table"

function.biom <- read_biom("./Subsystems function.biom")

function.matrix <- as.matrix(biom_data(function.biom))

function.matrix <- t(function.matrix)

function.matrix   <- function.matrix[c("0_ppm_r1",
                                       "0_ppm_r2",
                                       "0_ppm_r3",
                                       "60_ppm_r1",
                                       "60_ppm_r2",
                                       "60_ppm_r3",
                                       "833_ppm_r1",
                                       "833_ppm_r2",
                                       "833_ppm_r3",
                                       "2000_ppm_r1",
                                       "2000_ppm_r2",
                                       "2000_ppm_r3"), ]

colnames(function.matrix) <- paste(observation_metadata(function.biom)[[1]], ":",
                                   observation_metadata(function.biom)[[2]] , ":",
                                   observation_metadata(function.biom)[[3]] , ":",
                                   observation_metadata(function.biom)[[4]])

namesOfAllSubsystems <- paste(observation_metadata(function.biom)[[1]], ":",
                              observation_metadata(function.biom)[[2]] , ":",
                              observation_metadata(function.biom)[[3]] , ":",
                              observation_metadata(function.biom)[[4]])

map <- data.frame(row.names = c("0_ppm_r1",
                                "0_ppm_r2",
                                "0_ppm_r3",
                                "60_ppm_r1",
                                "60_ppm_r2",
                                "60_ppm_r3",
                                "833_ppm_r1",
                                "833_ppm_r2",
                                "833_ppm_r3",
                                "2000_ppm_r1",
                                "2000_ppm_r2",
                                "2000_ppm_r3"),
                  Silver_concentration=factor(c("0","0","0",
                                                "60","60","60",
                                                "833","833","833",
                                                "2000","2000","2000"),
                                              levels = c("0","60","833","2000")),
                  replicate=rep(c("1","2","3"),4))

# Get into Phyloseq
otumat <- as(biom_data(function.biom), "matrix")
otumat <- otumat[,c("0_ppm_r1",
                    "0_ppm_r2",
                    "0_ppm_r3",
                    "60_ppm_r1",
                    "60_ppm_r2",
                    "60_ppm_r3",
                    "833_ppm_r1",
                    "833_ppm_r2",
                    "833_ppm_r3",
                    "2000_ppm_r1",
                    "2000_ppm_r2",
                    "2000_ppm_r3")]
OTU    <- otu_table(otumat, taxa_are_rows=TRUE)
taxmat <- as.matrix(observation_metadata(function.biom), rownames.force=TRUE)
TAX    <- tax_table(taxmat)
physeq <- phyloseq(OTU, TAX, sample_data(map))

# Run some transformations and make dataframe of tax table
taxdf <- as.data.frame(tax_table(physeq))
taxdf$taxid <- row.names(taxdf)
ps_transformed <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

metalresistance <-subset_taxa(physeq, functionalHierarchy3=="Cobalt-zinc-cadmium_resistance")
resistance_anti_and_tox <-subset_taxa(physeq, functionalHierarchy2=="Resistance to antibiotics and toxic compounds")
oxstress <-subset_taxa(ps_transformed, functionalHierarchy2=="Oxidative stress")

searchterm=c("nitro")
searched_names <- subset(taxdf, grepl(searchterm, taxdf$functionalHierarchy1, ignore.case=T) | grepl(searchterm, taxdf$functionalHierarchy2, ignore.case=T) | grepl(searchterm, taxdf$functionalHierarchy3, ignore.case=T)| grepl(searchterm, taxdf$functionalHierarchy4, ignore.case=T))
searched_tax_ids <- row.names(searched_names)
physeq_searched = prune_taxa(searched_tax_ids, ps_transformed)
physeq_searched_glom <- tax_glom(physeq_searched, taxrank="functionalHierarchy4")

plot_heatmap(resistance_anti_and_tox, sample.order=row.names(map),   taxa.label="functionalHierarchy4", sample.label="Silver_concentration", trans=log_trans(2))
plot_heatmap(resistance_anti_and_tox, sample.order=row.names(map),   taxa.label="functionalHierarchy3", sample.label="Silver_concentration", trans=log_trans(2))
plot_bar(physeq_searched_glom, x="Silver_concentration", fill="functionalHierarchy4")


### dput(as.character(yourVector)) ## Useful way to get things formatted after pasting in...
### These are genes within Level 3 pathways of BMDExpress passing the filters
genesOfInterest <- c("3696","2039","1623","1072","166","163",
                     "3530","3529","3528","1071","108","2570",
                     "723","722","639","3301","2796","2794",
                     "1045","1070","542","105","2018","764",
                     "763","1237","1227","1216","1972","1821",
                     "1571","3638","3379","1632")

### These are all the genes that passed strict filters for BMD analysis
genesOfInterest_BMD_passed_filters  <- c("1263", "2781", "3653", "904", "2820", "443", "1742", "1401", 
                                         "3441", "1116", "1369", "2076", "2454", "1641", "3339", "1608", 
                                         "2854", "1738", "3210", "1066", "1952", "1630", "3276", "1946", 
                                         "1209", "2189", "2027", "1635", "1264", "3569", "199", "3668", 
                                         "918", "3225", "2940", "2433", "2239", "1848", "3656", "3135", 
                                         "1892", "725", "3157", "3344", "3209", "219", "3204", "3249", 
                                         "14", "840", "2572", "476", "3306", "3649", "3434", "2469", "1389", 
                                         "2435", "1967", "1622", "36", "185", "874", "2481", "2213", "2445", 
                                         "677", "938")

### These are all the genes with a modeled BMD <2000 ppm
genesOfInterest_broad <- c("3441", "3209", "677", "1263", "1021", "904", "2820", "1116", 
                           "1264", "1641", "918", "14", "3135", "1814", "3668", "36", "1635", 
                           "2940", "235", "3276", "2076", "2379", "323", "3653", "1066", 
                           "1401", "874", "3344", "1738", "1946", "1045", "1630", "3157", 
                           "2027", "1369", "1952", "3126", "1667", "854", "2239", "1571", 
                           "1257", "2785", "3129", "2796", "1608", "3693", "3204", "3312", 
                           "2454", "443", "2377", "2572", "1892", "840", "938", "1216", 
                           "848", "672", "863", "1926", "655", "3370", "906", "1129", "1838", 
                           "2433", "1010", "1130", "3573", "1848", "3473", "3225", "1209", 
                           "2189", "542", "3649", "246", "2459", "1140", "725", "1389", 
                           "3515", "1185", "1643", "1133", "1851", "1821", "2493", "2435", 
                           "2322", "1967", "1562", "105", "3580", "1742", "185", "2399", 
                           "2488", "1832", "1766", "1816", "808", "1501", "512", "1237", 
                           "2039", "2807", "2060", "2402", "3210", "1972", "1230", "33", 
                           "1794", "2068", "1003", "773", "1662", "2101", "3346", "1633", 
                           "824", "1726", "1831", "2213", "3569", "63", "199", "2481", "228", 
                           "1227", "2570", "428", "3434", "1451", "3696", "2306", "2445", 
                           "1259", "371", "2226", "2244", "2238", "166", "2898", "81", "1368", 
                           "277", "332", "445", "476", "163", "836", "3306", "3301", "2059", 
                           "516", "2018", "32", "3339", "2142", "2844", "219", "600", "3638", 
                           "3671", "2023", "1485", "2802", "2466", "763", "3404", "2298", 
                           "331", "1120", "1590", "3699", "1826", "2206", "2794", "2363", 
                           "2686", "454", "678", "1623", "2804", "2431", "2432", "1489", 
                           "2284", "2564", "1174", "2364", "1703", "2222", "2758", "2854", 
                           "2468", "1166", "316", "1923", "2781", "1347", "2576", "2592", 
                           "135", "1690", "1622", "2700", "2720", "2225", "1473", "3289", 
                           "1070", "2469", "2359", "3471", "1919", "3104", "843", "764", 
                           "1985", "1072", "1936", "2889", "531", "587", "3544", "893", 
                           "2614", "3018", "761", "108", "1071", "3528", "3154", "3529", 
                           "3530", "766", "2320", "1519", "575", "2688", "1153", "1605", 
                           "2025", "2218", "3340", "3249", "2973", "1941", "3170", "2504", 
                           "3379", "3656", "202", "21", "2818", "2393", "3430", "87", "1858", 
                           "3428", "1275", "2869", "2739", "2922", "1576", "2419", "3659", 
                           "1685", "3582", "3432", "710", "733", "2321", "3199", "1735", 
                           "2928", "2024", "765", "804")

ps_genes_dose_response <- prune_taxa(genesOfInterest, ps_transformed)
ps_genes_dose_response2 <- prune_taxa(genesOfInterest_BMD_passed_filters, ps_transformed)
ps_genes_dose_response3 <- prune_taxa(genesOfInterest_broad, ps_transformed)
# Order by expression @ 0 ppm AgNP
genesOfInterest_ordered <- names(sort(rowMeans(ps_genes_dose_response@otu_table@.Data[,1:3])))
genesOfInterest_BMD_passed_filters_ordered <- names(sort(rowMeans(ps_genes_dose_response2@otu_table@.Data[,1:3])))
genesOfInterest_BMD_passed_filters_ordered2 <- names(sort(rowMeans(ps_genes_dose_response3@otu_table@.Data[,1:3])))

dr_heatmap1 <- plot_heatmap(ps_genes_dose_response,
             sample.order=row.names(map),
             taxa.order=rev(genesOfInterest_ordered),
             taxa.label="functionalHierarchy4",
             sample.label="Silver_concentration")
            
dr_heatmap1$scales$scales[[1]]$name <- "Silver concentration"
dr_heatmap1$scales$scales[[2]]$name <- "Gene"
print(dr_heatmap1)
 
dr_heatmap2 <- plot_heatmap(ps_genes_dose_response2,
             sample.order=row.names(map),
             taxa.order=rev(genesOfInterest_BMD_passed_filters_ordered),
             taxa.label="functionalHierarchy4",
             sample.label="Silver_concentration")
dr_heatmap2$scales$scales[[1]]$name <- "Silver concentration"
dr_heatmap2$scales$scales[[2]]$name <- "Gene"
print(dr_heatmap2)

dr_heatmap3 <- plot_heatmap(ps_genes_dose_response3,
             sample.order=row.names(map),
             taxa.order=rev(genesOfInterest_BMD_passed_filters_ordered2),
             taxa.label="functionalHierarchy4",
             sample.label="Silver_concentration",
             max.label=300)
dr_heatmap3$scales$scales[[1]]$name <- "Silver concentration"
dr_heatmap3$scales$scales[[2]]$name <- "Gene"
print(dr_heatmap3)

####################################################################################
# DESeq2 ###########################################################################
####################################################################################

mainFactor="Silver_concentration"

ps <- physeq
dds <- phyloseq_to_deseq2(ps, as.formula(paste0("~",mainFactor)))

taxmat

# Change factor levels if needed
# colData(dds)[[mainFactor]] <- factor(colData(dds)[[mainFactor]], levels=myLevels)

# Geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Calculate geometric means
geoMeans = apply(counts(dds), 1, gm_mean)

# Run DESeq2
dds = estimateSizeFactors(dds)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

dds = DESeq(dds, fitType="local")

res = results(dds)
resultsNames(dds)
DESeq2::plotMA(dds)

## Custom contrasts...

res60 = results(dds, contrast=c("Silver_concentration","60","0"))
res833 = results(dds, contrast=c("Silver_concentration","833","0"))
res2000 = results(dds, contrast=c("Silver_concentration","2000","0"))

par(mfrow=c(2,2), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)

DESeq2::plotMA(res60, main="60 ppm vs. Control")
DESeq2::plotMA(res833, main="833 ppm vs. Control")
DESeq2::plotMA(res2000, main="2000 ppm vs. Control")

par(mfrow=c(1,1))

plotCounts(dds, gene=which.min(res$padj), intgroup="Silver_concentration")
plotCounts(dds, gene=804, intgroup="Silver_concentration")
plotCounts(dds, gene=15, intgroup="Silver_concentration")

par(mfrow=c(6,6), mar=c(1,1,1,1))
for (x in genesOfInterest) {
    plotCounts(dds, gene=x, intgroup="Silver_concentration")
}
par(mfrow=c(1,1))

################# Volcano ########################
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange,
     -log10(res$padj),
     col=cols,
     panel.first=grid(),
     main="Volcano plot",
     xlab="log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=20,
     cex=1,
     xlim=c(-7,8))
abline(v=0)
abline(v=c(-1.5,1.5), col="black")
abline(h=-log10(alpha), col="black")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected], cex=0.4) # taxdf$functionalHierarchy4[c(1,2)]

##################################################
# ggplot2 Volcano
ggplot(as.data.frame(res60), aes(x=log2FoldChange, y=-log10(res$padj))) +
  xlim(-10,10) +
  ylim(0,6) +
  geom_point(color="green") +
  geom_point(data=as.data.frame(res833), color="blue") +
  geom_point(data=as.data.frame(res2000), color="red")

##################################################

# Create list of DESeq results for each contrast (factor)...
resList <- list()
resList[[1]] <- res60
resList[[2]] <- res833
resList[[3]] <- res2000

# Filter results table using adjusted p-value of alpha...
# alpha = 0.01 # Also defined above for Volcano plot, so just pick once to be consistent
sigtabList <- list()

# replace "log2 fold change (MLE): Group" with "LFC"

for (i in 1:length(resList)) {
  print(i)
  sigTab <- resList[[i]]
  # Add taxonomy
  if (nrow(sigTab) == 0) {
    message("SigTab has no rows, skipping...")
    next
  } else {
    newContrast=gsub("log2 fold change \\(MLE\\)\\: Group ",
                     "",
                     sigTab@elementMetadata$description[[2]])
    sigTab <- cbind(as(sigTab, "data.frame"),
                    as(tax_table(ps)[rownames(sigTab), ], "matrix"),
                    contrast=newContrast)
    sigTab <- sigTab[!is.na(sigTab$padj) & sigTab$padj < alpha, ]
    sigtabList[[i]] <- sigTab
  }
}
sigtabList <- sigtabList[!sapply(sigtabList, is.null)] 

summaryTable <- data.frame( baseMean=resList[[1]]$baseMean )
for (i in 1:length(resList)) {
  print(i)
  n <- gsub("log2 fold change \\(MLE\\)\\: Group ",
            "",
            resList[[i]]@elementMetadata$description[[2]])
  p <- resList[[i]]@elementMetadata$description[[6]]
  message(n)
  summaryTable <- cbind(summaryTable, log2FoldChange=resList[[i]]$log2FoldChange, padj=resList[[i]]$padj)
  names(summaryTable)[[ncol(summaryTable)-1]] <- n
  names(summaryTable)[[ncol(summaryTable)]] <- p
}
summaryTable <- cbind(summaryTable, as(tax_table(ps)[rownames(summaryTable), ], "matrix"))

##############
# Plot results
theme_set(theme_bw())

# Write dataframe of all results
significantResults <- do.call(rbind, sigtabList)
#######################################
### Write results table from DESeq2
#######################################
write.table(significantResults, file="DESeq_output_significant.txt", quote=F, sep='\t', col.names=NA)
write.table(summaryTable, file="DESeq_output_all_genes.txt", quote=F, sep='\t', col.names=NA)
#######################################

# Create list of results tables where NA genera are removed
#sigtableListNA.RM <-lapply(sigtabList, subset, !is.na(Genus))

# Create list of combined results
#combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)

attach(significantResults)
significantResultsOrdered <- significantResults[order(padj),]
significantResultsOrderedTop20 <- head(significantResultsOrdered, n=20)

topResults <- NULL
numResults=100
#contrastsOfInterest <- levels(significantResultsOrdered$contrast)[c(3,4)]
contrastsOfInterest <- levels(significantResultsOrdered$contrast)
for (i in contrastsOfInterest) {
  topResults <- rbind(topResults,
                      head(significantResultsOrdered %>% 
                            dplyr::filter(contrast==i),
                           n=numResults))
}

topResults <- topResults %>%
  dplyr::mutate(functionalHierarchy4 = forcats::fct_reorder(functionalHierarchy4,
                                                            log2FoldChange)) %>%
  dplyr::arrange(log2FoldChange)

levels(topResults$contrast) <- c("60 ppm AgNP",
                                 "833 ppm AgNP",
                                 "2000 ppm AgNP")
## Genus plots, no N/A
ggplot(topResults,
       aes(x=functionalHierarchy4,
           y=log2FoldChange,
           color=functionalHierarchy1,
           size=-log(padj))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  facet_grid(contrast ~ ., scales="free_x") +
  ggtitle(paste0("Top ",
                 numResults,
                 " genes with significantly different relative abundance (adjust p-value <",
                 alpha,
                 "), grouped by treatment"))
# theme(legend.position = "none")

### Nice
ggplot(topResults,
       aes(x=log2FoldChange,
           y=functionalHierarchy4,
           color=functionalHierarchy1,
           size=-log(padj))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  facet_grid(contrast ~ ., scales="free_y") +
  ggtitle(paste0("Top ",
                 numResults,
                 " genes with significantly different relative abundance (adjust p-value <",
                 alpha,
                 "), grouped by treatment")) +
  theme(legend.position = "none") +
  geom_vline(xintercept=0,
             linetype="dashed",
             color = "black",
             size=1)

taxdf$rowname <- row.names(taxdf)
searched_names <- taxdf[taxdf$functionalHierarchy4 %in% topResults$functionalHierarchy4, ]
searched_tax_ids <- row.names(searched_names)
physeq_searched=prune_taxa(searched_tax_ids, ps_transformed)
physeq_searched_glom <- tax_glom(physeq_searched, taxrank="functionalHierarchy4")

taxdfOrderedByLog2FC <- dplyr::full_join(taxdf,
                                         topResults,
                                         by = "functionalHierarchy4") %>%
  dplyr::arrange(log2FoldChange)

taxOrder <- taxdfOrderedByLog2FC[taxdfOrderedByLog2FC$functionalHierarchy4 %in%
                                   topResults$functionalHierarchy4, ]$taxid

plot_heatmap(physeq_searched_glom,
             taxa.label="functionalHierarchy4",
             taxa.order=unique(taxOrder),
             sample.label="Silver_concentration",
             sample.order=row.names(map),
             trans=log_trans(3))

plot_heatmap(physeq_searched_glom,
             taxa.label="functionalHierarchy4",
             taxa.order=unique(taxOrder),
             sample.label="Silver_concentration",
             sample.order=row.names(map))


## Genus plots, no N/A
ggplot(topResults,
       aes(x=functionalHierarchy4,
           y=log2FoldChange,
           color=functionalHierarchy3,
           size=-log(padj))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  facet_grid(~contrast, scales="free_x") +
  ggtitle(paste0("Top ",
                 numResults,
                 " genes with significantly different relative abundance (adjust p-value <",
                 alpha,
                 "), grouped by treatment")) +
  theme(legend.position = "none")


#### FIGURE 4 #####

topResults$functionalHierarchy4_wrapped <- factor(stringr::str_wrap(topResults$functionalHierarchy4,
                                                                70),
                                              levels=stringr::str_wrap(levels(topResults$functionalHierarchy4),
                                                                           70))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

ggplot(topResults,
       aes(x=log2FoldChange,
           y=functionalHierarchy4_wrapped,
           color=functionalHierarchy1,
           size=-log(padj))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,
                                   hjust = 0,
                                   vjust=0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold")) +
  facet_grid(~contrast) +
  #labs(tag="B.") +
  ggtitle(paste0("Genes with significantly different relative \nabundance (adjusted p-value <",
                 alpha,
                 ")")) +
  #theme(legend.position = "none") +
  geom_vline(xintercept=0,
             linetype="dashed",
             color = "black",
             size=1) +
  scale_color_manual(name = "SEED subsystem 1",
                     values = getPalette(15)) +
  scale_size(name = "-log(adjusted p-value)") +
  guides(colour = guide_legend(override.aes = list(size=3)))

# Examine specific functions...
# myPS<-ps
# searchterm="Predicted cobalt transporter in Bacteroides_Porphyromonas"
# num <- grep(searchterm, myPS@tax_table@.Data[,4])
# myPS@otu_table@.Data[num,]
# dds@assays[[1]][num,]

######## BMDExpress2 #######
read.counts.sf_normalized <- counts(dds, normalized=T)
lognorm.read.counts <- log2(read.counts.sf_normalized + 1)

dds.rlog <- rlog(dds, blind=T) ### May need to set blind=F if large global changes are observed
rlog.norm.counts <- assay(dds.rlog)

par(mfrow=c(4,1))
boxplot(counts(dds, normalized=F),
        notch=T,
        main="Untransformed read counts",
        ylab="counts")

boxplot(read.counts.sf_normalized,
        notch=T,
        main="Size factor normalized read counts",
        ylab="counts")

boxplot(lognorm.read.counts,
        notch=T,
        main="log2-Transformed read counts",
        ylab="counts")

boxplot(rlog.norm.counts,
        notch=T,
        main="rlog-normalized read counts",
        ylab="counts")

par(mfrow=c(1,1))

bmdexpress <- as.data.frame(lognorm.read.counts) # log2 normalized, size factor normalized, 1 added
head(bmdexpress)
bmdexpress <- cbind(SampleID=row.names(bmdexpress), bmdexpress, stringsAsFactors=F)
head(bmdexpress)
bmdexpress <- rbind( Dose=c("Dose",
                            as.character(map$Silver_concentration)),
                     bmdexpress,
                     stringsAsFactors=F)
# bmdexpress <- rbind( Dose=c("Dose","0","1","2"), bmdexpress, stringsAsFactors=F)
head(bmdexpress)
write.table(bmdexpress, file = "bmdexpress_input.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Need two files, probe map and category map
# For probe map:

# These are the gene numbers recognized from the BMDexpress output above, aka, "Array probe"
# rownames(observation_metadata(function.biom))
# Metadata column 4 is the most detailed, so it will be the "transcript/gene" level ID, aka, "Category component"
# observation_metadata(function.biom)[[4]]

arrayProbe <- data.frame(as.numeric(rownames(observation_metadata(function.biom))),
                         as.vector(observation_metadata(function.biom)[4]))
colnames(arrayProbe) <- c("Array Probe", "Category Component")

# For category map, level 3 of the SEED subsystem metadata

categoryMap <- data.frame(as.vector(observation_metadata(function.biom)[3]),
                         as.vector(observation_metadata(function.biom)[3]),
                         as.vector(observation_metadata(function.biom)[4]))
colnames(categoryMap) <- c("Category ID", "Category Name", "Category Component")

categoryMap2 <- data.frame(as.vector(observation_metadata(function.biom)[2]),
                          as.vector(observation_metadata(function.biom)[2]),
                          as.vector(observation_metadata(function.biom)[4]))
colnames(categoryMap2) <- c("Category ID", "Category Name", "Category Component")

categoryMap1 <- data.frame(as.vector(observation_metadata(function.biom)[1]),
                           as.vector(observation_metadata(function.biom)[1]),
                           as.vector(observation_metadata(function.biom)[4]))
colnames(categoryMap1) <- c("Category ID", "Category Name", "Category Component")

write.table(arrayProbe, file = "arrayProbe.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(categoryMap, file = "categoryMap.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(categoryMap2, file = "categoryMapLevel2.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(categoryMap1, file = "categoryMapLevel1.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Overlap between DESeq2 and BMDExpress

overlap_DESeq2_BMDexpressDCA <- intersect(genesOfInterest_ordered, row.names(significantResults))
length(intersect(row.names(significantResults), genesOfInterest_ordered))

overlap_DESeq2_BMDexpressRaw <- intersect(genesOfInterest_BMD_passed_filters_ordered2, row.names(significantResults))
length(intersect(row.names(significantResults), genesOfInterest_broad))

ps_genes_overlap <- prune_taxa(intersect(row.names(significantResults),
                                         genesOfInterest_broad),
                               ps_transformed)

ps_genes_overlap2 <- prune_taxa(intersect(row.names(significantResults),
                                          genesOfInterest_ordered),
                                ps_transformed)

plot_heatmap(ps_genes_overlap,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             taxa.order=rev(overlap_DESeq2_BMDexpressRaw),
             sample.label="Silver_concentration") + 
  theme(axis.text.y = element_text(size=14))

plot_heatmap(ps_genes_overlap2,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             taxa.order=rev(overlap_DESeq2_BMDexpressDCA),
             sample.label="Silver_concentration") + 
  theme(axis.text.y = element_text(size=14))

#### FIGURE - Not used ####

fig.ord <- ordinate(ps_transformed,
                     method="PCoA",
                     distance="jaccard")

evals <- fig.ord$values$Eigenvalues

fig <- plot_ordination(ps_transformed, fig.ord,
                        color = "Silver_concentration")
#colnames(fig1$data)[1:2] <- c("PC1", "PC2")

figdata <- fig$data
# ggplot(fig3data, aes(x=PC1, y=PC2, color=Silver_Concentration)) +
fig +
  geom_point(size=2.5, shape="square") +
  scale_color_manual(values=brewer.pal(8, "Dark2")[c(1,2,5,6)]) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(color = "Silver concentration", tag="A.") +
  ggtitle("PCoA of metatranscriptome functional annotations") +
  theme_bw()

  # brewer.pal(8, "Dark2")[c(1,2,5,6)]

bmdexpress_path_results <- read.csv("./BMDexpress pathway results - levels 2 and 3.txt", sep="\t", header=T)
bmdexpress_path_results_filt <- filter(bmdexpress_path_results, Input.Genes > 3)

pathway_chart_data <- read.csv("./chartdataExport_pathway.txt", sep="\t", header=T)
ggplot(pathway_chart_data, aes(x=x,
                               y=y,
                               color=series,
                               label=components.delimited.by....)) + 
  geom_point() +
  geom_text()

summarize_phyloseq(ps)
summary(res60)
summary(res833)
summary(res2000)

rs <- rowSums(counts(dds))
rmx <- apply(counts(dds), 1, max)
plot(rs+1, rmx/rs, log="x")

BMDvsBMDL <- read.csv("./Best BMD vs Best BMDL.txt", sep="\t", header=T)
BMDvsBMDL$label <- gsub("\\s.*","",BMDvsBMDL$label)

#### Figure 5A ####

ggplot(BMDvsBMDL, aes(x=x, y=y, label=label)) + 
  geom_point(size=2) +
  theme_bw() +
  expand_limits(x = c(0,2000), y = c(0,2000)) +
  labs(tag = "A.") +
  ggtitle("BMD vs BMDL for 89 genes passing all filters") +
  scale_x_continuous(name = "Best BMD") +
  scale_y_continuous(name = "Best BMDL") + 
  coord_fixed(ratio = 1) #+
  #geom_label()

ps_genes_dose_response_filt <- prune_taxa(BMDvsBMDL$label, ps_transformed)

genesOfInterest_ordered <- names(sort(rowMeans(ps_genes_dose_response_filt@otu_table@.Data[,1:3])))
# For getting the names of the genes...
genesOfInterest_ordered_DF <- ps_genes_dose_response_filt@tax_table@.Data
write.table(file = "./genesForBMDExpressTable.txt", genesOfInterest_ordered_DF, quote = F, sep="\t")


ps_genes_dose_response_filt_BMD_gt_10xlowest <- prune_taxa(dplyr::filter(BMDvsBMDL, BMDvsBMDL$y > 1)$label,
                                                            ps_transformed)
plot_heatmap(ps_genes_dose_response_filt_BMD_gt_10xlowest,
             sample.order=row.names(map),
             taxa.order=rev(genesOfInterest_ordered),
             taxa.label="functionalHierarchy4",
             sample.label="Silver_concentration") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  labs(tag="B.")
