genes <- c("PPARGC1A", "SOD1", "SOD2", "SOD3", "GPX1", "GPX2", "GPX3", 
           "GPX4,", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6", 
           "TXN1", "TRX1", "TRX", "TXN", "TRDX", "RXN", "TXN2", "TRX2", 
           "TXNRD1", "TXNRD2", "CAT", "GSR", "GSS", "NQO1")

require(xlsx2dfs)

df <- read.table("expressions_annotated_PM2.5_vs_control.csv", 
                    sep = "\t", 
                    header = T, 
                    row.names = 1, 
                    stringsAsFactors = FALSE)
df$DELTA <- df$treated_mean - df$control_mean
df$fc <- df$treated_mean / df$control_mean
df$log2FC <- log2(df$fc)
df$AveExp <- rowMeans(df[, 1:6 ])
collected.genes.df <- df[df$Gene.Symbol %in% genes, ]
collected.genes.df <- collected.genes.df[order(collected.genes.df$Gene.Symbol), ]

dfs2xlsx(withNames("point_5", collected.genes.df),
                   fpath = "collected_genes_values_20200212.xlsx")

gene_name_contained <- Vectorize(function(query_gene, subject_gene, delimiter) {
  any(query_gene == unlist(strsplit(x = subject_gene, split = delimiter)))
})

collect_rows_containing <- function(df, genes_vec, df_col, delimiter = " /// ") {
  is_containing <- Reduce(`|`, lapply(genes_vec, function(gene) {
    gene_name_contained(gene, df[, df_col], delimiter)
  }))
  df[is_containing, ]
}

collected.genes.df.corrected <- collect_rows_containing(df = df, 
                                                           genes_vec = genes, 
                                                           df_col = "Gene.Symbol", 
                                                           delimiter = " /// ")

head(collected.genes.df.corrected)
collected.genes.df.corrected <- 
  collected.genes.df.corrected[order(collected.genes.df.corrected$Gene.Symbol), ]
dfs2xlsx(withNames("table_to_point_5", 
                   collected.genes.df.corrected), 
         fpath = "collected_genes_values_200212.xlsx")


genes[genes %in% collected.genes.df.corrected$Gene.Symbol]
#  [1] "PPARGC1A" "SOD1"     "SOD2"     "SOD3"     "GPX1"     "GPX2"    
#  [7] "GPX3"     "PRDX1"    "PRDX2"    "PRDX3"    "PRDX4"    "PRDX5"   
# [13] "PRDX6"    "TXN"      "TXN2"     "TXNRD1"   "TXNRD2"   "CAT"     
# [19] "GSR"      "GSS"      "NQO1"    
# >
# 

genes[!(genes %in% collected.genes.df.corrected$Gene.Symbol)]
# [1] "GPX4," "TXN1"  "TRX1"  "TRX"   "TRDX"  "RXN"   "TRX2" 


df.1 <- collected.genes.df.corrected

genes <- c("BCL2", "BCL2L1", "XIAP", "BIRC5", "CFLAR", "MCL1", 
           "BAD", "BBC2", "BCL2L8", "BAX", "BCL2L4", "BAK1", "BAK", "BCL2L7", 
           "BIM", "BCL2L11", "BID", "HMGB1", "BECN1", "ATG5")

collected.genes.df.corrected <- collect_rows_containing(df = df, 
                                                        genes_vec = genes, 
                                                        df_col = "Gene.Symbol", 
                                                        delimiter = " /// ")
collected.genes.df.corrected <- 
  collected.genes.df.corrected[order(collected.genes.df.corrected$Gene.Symbol), ]

df.2 <- collected.genes.df.corrected

dfs2xlsx(withNames("table_to_point_5", df.1,
                   "table_to_point_6", df.2), 
         fpath = "collected_genes_values_200212_corrected.xlsx")

genes[!(genes %in% collected.genes.df.corrected$Gene.Symbol)]
# [1] "BBC2"   "BCL2L8" "BCL2L4" "BAK"    "BCL2L7" "BIM"    





















