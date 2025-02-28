#load packages
library(Seurat)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)

#set organism
organism = "org.Hs.eg.db"

#read in table of genes generated as output of 'FindAllMarkers' or 'FindMarkers'
geneTable <- read.csv('filePath/geneList.csv')

#ensure the table is a data frame and arrange by descending log fold change values
geneTable <- as.data.frame(geneTable) %>% arrange(desc(`avg_log2FC`))

#Map gene IDs
geneTable$GeneID <- mapIds(org.Hs.eg.db, keys = geneTable$column_Containing_Gene_Names, column = "ENTREZID", keytype = "SYMBOL")

#Filter out NA values and create seperate tables for upregulated and downregulated genes
geneTable <- geneTable %>% filter(!is.na(`GeneID`))
UPgenes <- subset(geneTable, avg_log2FC < 0)
DOWNgenes <- subset(geneTable, avg_log2FC > 0)

#create ranks
UPgeneRanks <- UPgenes$avg_log2FC
DOWNgeneRanks <- DOWNgenes$avg_log2FC

#add names to ranks
names(UPgenesRanks) <- UPgenes$GeneID
names(DOWNgenesRanks) <- DOWNgenes$GeneID

#format the ranks as character data
geneIDs <- as.character(names(UPgenesRanks))

#run enrichment analysis
clustGO <- enrichGO(gene = geneIDs,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")

#format analysis results as a dataframe and save
df <- as.data.frame(clustGO)
write.table(df, file = "filePath/fileName.txt", sep ="\t", row.names = TRUE, col.names = TRUE)

### re-run lines 32-44 for downregulated genes ###

#Create custom bubbleplot function
make_dotplot <- function(df, title="", ylabel="Description", colour="#56B1F7", n=10){
  df$ycolour <- "black"
  if ("ONTOLOGY" %in% colnames(df)){
    df$Description <- paste(df$ONTOLOGY, df$Description, sep=' - ')
    df$ycolour <- ifelse(grepl("BP -", df$Description), 'blue', df$ycolour)
    df$ycolour <- ifelse(grepl("CC -", df$Description), 'red', df$ycolour)
    df$ycolour <- ifelse(grepl("MF -", df$Description), 'darkgreen', df$ycolour)
  }
  df <- df[order(df$p.adjust, decreasing=FALSE),]
  plt <- ggplot() +
    geom_point(data=head(df, n=n),
               aes(x = -log(p.adjust), 
                   y = reorder(Description, -p.adjust), 
                   colour = Count, 
                   size = unname(unlist(sapply(GeneRatio, function(x) eval(parse(text=x)))))*100,
               ),
    ) + 
    theme_linedraw() +
    theme(axis.text.y = element_text(size = 20)) +
    scale_color_gradient(low = "black", high = colour) +
    ggtitle(title)
  plt$labels$x <- "-log(p.adjust)"
  plt$labels$y <- ylabel
  plt$labels$size <- "GenePercentage"
  plt$labels$colour <- "GeneCount"
  return(plt)
}

#plot enriched Gene ontology terms
#make the plot
make_dotplot(df, title = paste("Upregulated Gene Ontology Enrichment Analysis"))+
  scale_color_gradient(low = "blue", high = "red") + theme(panel.grid = element_blank())
