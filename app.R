library(Seurat)
library(ShinyCell)
#getExampleData("multi")

seu <- readRDS("readySeu_rset.rds")
scConf = createConfig(seu)
showLegend(scConf)


# Delete excessive metadata and rename some metadata
scConf = delMeta(scConf, c("orig.ident", "RNA_snn_res.0.5", "phase"))
scConf = modMetaName(scConf, 
                     meta.to.mod = c("nUMI", "nGene", "pctMT", "pctHK"), 
                     new.name = c("No. UMIs", "No. detected genes",
                                  "% MT genes", "% HK genes"))
showLegend(scConf)

# Modify colours and labels
scConf = modColours(scConf, meta.to.mod = "library", 
                    new.colours= c("black", "darkorange", "blue", "pink2"))
scConf = modLabels(scConf, meta.to.mod = "library", 
                   new.labels = c("fm", "pr", "nr", "rr"))
showLegend(scConf)

# To display the order of the meta data
showOrder(scConf)


# Add metadata back, reorder, default
scConf = addMeta(scConf, "phase", seu)
showOrder(scConf)
scConf = reorderMeta(scConf, scConf$ID[c(1:5,22,6:21)])
showOrder(scConf)
scConf = modDefault(scConf, "library", "identity")
showOrder(scConf)

# Build shiny app
checkConfig(scConf, seu)
citation = list(
  author  = "Liu X., Ouyang J.F., Rossello F.J. et al.",
  title   = "",
  journal = "Nature",
  volume  = "586",
  page    = "101-107",
  year    = "2020", 
  doi     = "10.1038/s41586-020-2734-6",
  link    = "https://www.nature.com/articles/s41586-020-2734-6")


makeShinyApp(seu, scConf, gene.mapping = TRUE, 
             gex.assay = "RNA", gex.slot = "data",
             shiny.title = "ShinyCell Tutorial",
             shiny.dir = "shinyApp/", shiny.footnotes = citation,
             default.gene1 = "NANOG", default.gene2 = "DNMT3L",
             default.multigene = c("ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
                                   "DPPA5","SLC7A2","GATA3","KRT19"))


# scConf1 = delMeta(scConf1, c("orig.ident", "RNA_snn_res.0.5"))
# scConf1 = modMetaName(scConf1, meta.to.mod = c("nUMI", "nGene", "pctMT", "pctHK"), 
#                       new.name = c("No. UMIs", "No. detected genes",
#                                    "% MT genes", "% HK genes"))
# scConf1 = modColours(scConf1, meta.to.mod = "library", 
#                      new.colours= c("black", "darkorange", "blue", "pink2"))
# makeShinyFiles(seu, scConf1, gex.assay = "RNA", gex.slot = "data",
#                gene.mapping = TRUE, shiny.prefix = "sc1",
#                shiny.dir = "shinyAppMulti/",
#                default.gene1 = "NANOG", default.gene2 = "DNMT3L",
#                default.multigene = c("ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
#                                      "DPPA5","SLC7A2","GATA3","KRT19"),
#                default.dimred = c("UMAP_1", "UMAP_2"))
