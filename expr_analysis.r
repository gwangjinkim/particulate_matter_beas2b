################################################################################
# INSTALL r-base and the simpleaffy package into a new CONDA ENVIRONMENT
################################################################################
# $ conda update -n base -c defaults conda
# $ conda create --name expr_analysis
# $ conda activate expr_analysis
# $ conda install -c bioconda bioconductor-simpleaffy
# $ # conda install -c bioconda bioconductor-limma # not necessary
# $ # conda install -c bioconda libopenblas        # not necessary
# # but very important:
# https://support.bioconductor.org/p/122925/
# $ git clone https://github.com/bmbolstad/preprocessCore.git
# $ cd preprocessCore/
# $ conda activate exr_analysis # install into R in the conda environment!
# $ R CMD INSTALL --configure-args="--disable-threading"  .
# ## this is crucial for simpleaffy::read.affy() and simpleaffy::call.exprs()
# ## otherwise they spit out errors like:
# ## libopenblas.so.0: cannot open shared object file: No such file or directory
# ## and: return code from pthread_create() is 22 error
# ## repectively!
################################################################################

# install the packages for Excel output if they are not installed
if (!"openxlsx" %in% rownames(installed.packages())) install.packages("openxlsx")
if (!"xlsx2dfs" %in% rownames(installed.packages())) install.packages("xlsx2dfs")

################################################################################

# transfer CEL files from GEO into the same folder where this script is located.
# also make sure that the 'myphenotypes1.csv' file and the annotation file
# 'HG-U133_Plus_2.na35.annotation.mod.csv' are present in this same folder.
# Open R console in this folder (or `setwd()` to it).
# (simpleaffy unfortulately cannot handle absolute paths correctly)
dir_path <- '.' 
fname <- 'myphenotypes1.tab' # CELs must be in same folder!

# setwd(dir_path) # phentype and CEL files must be in same folder
raw.data <- simpleaffy::read.affy(fname)        
eset <- simpleaffy::call.exprs(raw.data, "rma")
data <- eset@assayData$exprs # data = rma-noramlized expression data

#####################################
# bringing annotation and expression data together
#####################################

# read-in annotation data
annot <- read.csv('HG-U133_Plus_2.na35.annotation.modified.csv',
                  header=TRUE,
                  stringsAsFactors=FALSE) ## dim(annot) # [1] 54675 41

# merge data with annotation
rownames(annot) <- as.character(annot$Probe.Set.ID)
data_annot <- merge(data, affy_ids_annot, by="row.names")

# extract only: values, affy_id, and symbols
# with group means and group standard deviations as extra columns
df <- data_annot[, c(2:8, 22)]
res_df <- cbind(df[, c(1:6)],
                control_mean = apply(df[, c(1:3)], 1, mean),
                control_std  = apply(df[, c(1:3)], 1, sd),
                treated_mean = apply(df[, c(4:6)], 1, mean),
                treated_std  = apply(df[, c(4:6)], 1, sd),
                df[, c(7:8)])

# write results into excel
xlsx2dfs::dfs2xlsx(xlsx2dfs::withNames("expr_annot_PM2.5_vs_control", res_df),
                   "expressions_annotated_PM2.5_vs_control_2008_03_17.xlsx")


################################################################################
#### the following is usually the way to translate affy_ids to symbols - but 
#### the package doesn't contain the data which have "_s_" in its affy ids
#### (non-unique data).
#### Thus, for this project, to include those, I decided to download the annotation file 
#### from the company site and to annotate manually.
#### However, in another situtation the following code might be useful:
# BiocManager::install("hgu133plus2.db")
# affy2symbol_list <- function(affy_ids) {
#    x <- hgu133plus2.db::hgu133plus2SYMBOL
#    mapping <- as.list(x[AnnotationDbi::mappedkeys(x)])
#    lapply(affy_ids, function(id) mapping[[id]])
# } # translator function affy_id to symbol
# # usage:
# affy2symbol_list(c("1007_s_at", "1053_at", "1255_g_at", "1294_at"))
################################################################################
