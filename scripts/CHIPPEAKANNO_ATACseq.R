library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(Vennerable)
library(AnnotationHub)
library(reactome.db)
library(org.Hs.eg.db)

argv<-commandArgs(T)
macs_peak <- argv[1]
macs_peak_gr <- toGRanges(macs_peak, format = "BED", header=FALSE)
anno_data <- genes(EnsDb.Hsapiens.v86)
annotated_peak4 <- annotatePeakInBatch(macs_peak_gr,
                                       AnnotationData = anno_data,
                                       output = "both")
write.csv(annotated_peak4,"diff/annotated_peak.csv",row.names=FALSE)
enriched_go <- getEnrichedGO(annotated_peak4, 
                             orgAnn = "org.Hs.eg.db", 
                             feature_id_type = "ensembl_gene_id")
p<-enrichmentPlot(enriched_go, label_wrap = 60)
pdf(file = paste0("diff/enrichmentPlot_GO.pdf"),
    width = 8, height = 8, useDingbats=F)
suppressWarnings(print(p))
invisible(dev.off())
enriched_path <- getEnrichedPATH(annotated_peak4,
                                 orgAnn = "org.Hs.eg.db",
                                 feature_id_type = "ensembl_gene_id",
                                #  pathAnn = "KEGGREST")
                                 pathAnn = "reactome.db")
p<-enrichmentPlot(enriched_path)
pdf(file = paste0("diff/enrichmentPlot_PATH.pdf"),
    width = 8, height = 8, useDingbats=F)
suppressWarnings(print(p))
invisible(dev.off())
