#########################################################
# Filter expression table for gene groups of interest
#########################################################

# install the packages for Excel output if they are not installed
if ("openxlsx" %in% rownames(installed.packages()) install.packages("openxlsx")
if ("xlsx2dfs" %in% rownames(installed.packages()) install.packages("xlsx2dfs")

# genes of interest (goi)
goi.1 <- c("PPARGC1A", "SOD1", "SOD2", "SOD3", "GPX1", "GPX2", "GPX3", 
           "GPX4", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6", 
           "TXN1", "TRX1", "TRX", "TXN", "TRDX", "RXN", "TXN2", "TRX2", 
           "TXNRD1", "TXNRD2", "CAT", "GSR", "GSS", "NQO1")
goi.2 <- c("BCL2", "BCL2L1", "XIAP", "BIRC5", "CFLAR", "MCL1", 
             "BAD", "BBC2", "BCL2L8", "BAX", "BCL2L4", "BAK1", "BAK", "BCL2L7", 
             "BIM", "BCL2L11", "BID", "HMGB1", "BECN1", "ATG5")

# predicate function for filtering
gene_name_contained <- Vectorize(function(query_gene, subject_gene, delimiter) {
  any(query_gene == unlist(strsplit(x = subject_gene, split = delimiter)))
})

# filter function
filter_for_gois <- function(df, gois, col, delimiter = " /// ") {
  is_containing <- Reduce(`|`, lapply(gois, function(gene) gene_name_contained(gene, df[, col], delimiter)))
  df <- df[is_containing, ]
  # return alphabetically sorted table (by Gene.Symbol)
  df[order(df$Gene.Symbol), ]
}

# add DELTA fc log2FC AveExp
create_columns <- function(df) {
  df$DELTA <- df$treated_mean - df$control_mean
  df$fc <- df$treated_mean / df$control_mean
  df$log2FC <- log2(df$fc)
  df$AveExp <- rowMeans(df[, 1:6 ])
  df
}

# read-in expression data from excel file
df <- xlsx2dfs::xlsx2dfs("expressions_annotated_PM2.5_vs_control.xlsx")[[1]]
                                      
# eliminate uninteresting columns
cols <- c("K1_15.5_Kontrolle_200814_133Plus_2.CEL", "M1_15.5_Kontrolle_110914_133Plus_2.CEL", 
          "R1_15.5_Kontrolle_190814_133Plus_2.CEL", "K5_15.5_Filtergut_Kontrolle_220814_133Plus_2.CEL", 
          "M5_15.5_Filtergut_Kontrolle_101114_133_Plus_2.CEL", "R5_15.5_Filtergut_Kontrolle_061114_133_Plus_2.CEL", 
          "control_mean", "control_std", "treated_mean", "treated_std", "Probe.Set.ID", "Gene.Symbol")
df <- df[, cols]
                                      
# calculate and add DELTA fc log2FC AveExp columns 
df <- create_columns(df)
    
# filter for rows containing the goi gene groups
df.1 <- collect_rows_containing(df = df, genes_vec = genes, df_col = "Gene.Symbol", delimiter = " /// ")
df.2 <- collect_rows_containing(df = df, genes_vec = genes.1, df_col = "Gene.Symbol", delimiter = " /// ")

# print them out sheet-wise
xlsx2dfs::dfs2xlsx(withNames("table_to_point_5", df.1,
                             "table_to_point_6", df.2), 
                   fpath = "collected_genes_values_200212_corrected.xlsx")





















