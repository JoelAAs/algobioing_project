library(biomaRt)
library(tidyverse)
library(magrittr)
library(KEGGREST)
library(reshape2)

# Wrapper for SNPFilter from SNPmart
# RETURN: Dataframe
snp_getBM_wrap <- function(
  my_rsid, 
  snp_mart, 
  capturing_columns
  ) {
  return(
    getBM(attributes = capturing_columns,
    filter = "snp_filter",
    values = my_rsid,
    mart = snp_mart
    )
  )
}

# Wrapper for GeneFilter from SNPmart
# RETURN: Dataframe
gene_getBM_wrap <- function(
  my_gene_ids,
  gene_mart,
  capturing_columns
  ) {
  return(
    getBM(attributes = capturing_columns,
    filter = "ensembl_gene_id",
    values = my_gene_ids,
    mart = gene_mart
    )
  )
}

# Annotate rs_ids in association, those wihtout gene annotation are dropped
# RETURN: association dataframe with gene annotation
get_gene_id <- function(
  assoc_df,
  snp_column = "SNP",
  snp_columns = c("refsnp_id", "ensembl_gene_stable_id"),
  gene_columns = c("ensembl_gene_id", "external_gene_name")
  ) {
  # Get geneame from variant  
  snp_mart <- useMart(
    biomart = "ENSEMBL_MART_SNP",
    dataset = "hsapiens_snp",
    host = "http://grch37.ensembl.org")
  
  gene_mart  <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "http://grch37.ensembl.org")
  
  snp_gene_id  <- snp_getBM_wrap(assoc_df[, snp_column], snp_mart, snp_columns) %>% 
    {.[.$ensembl_gene_stable_id != "", ]}
  gene_ids  <- snp_gene_id[, "ensembl_gene_stable_id"] %>% 
    unique()
  gene_names <- gene_getBM_wrap(gene_ids, gene_mart, gene_columns)
  
  
  # Merge and drop variants with no annotation
  snp_gene_id[, "ensembl_gene_id"] <- snp_gene_id[, "ensembl_gene_stable_id"]
  snp_gene_id[, "ensembl_gene_stable_id"] <- NULL
  refsnp_to_gene <- merge(
    x=snp_gene_id,
    y=gene_names,
    by="ensembl_gene_id"
  ) %>% 
    {.[!duplicated(.[, c("refsnp_id", "external_gene_name")]), ]}
  
  assoc_df[, "refsnp_id"] <- assoc_df[, snp_column]
  assoc_df <- merge(
    x = assoc_df,
    y = refsnp_to_gene,
    by = "refsnp_id",
    all.y = TRUE)
  
  return(assoc_df)
}

# function needed as KEGGREST only returns 
#
get_keggs_from_genes <- function(terms, organism){
  terms <- unique(terms)
  query_terms <- sapply(terms, function(x) paste("hsa", x, sep=":"))
  query_terms %>% 
    unique() %>% 
    {split(., f=1:(length(.)/10 + 1))} -> term_split
 
  query_results  <- sapply(term_split, keggGet) 
  kegg_gene_pair <- list() 
  for (query_set in query_results) {
    for (query in query_set){
      gene_names <- strsplit(query$NAME, ", ")[[1]]
      selected_terms <- terms[terms %in% gene_names][1]
      brite <- query$BRITE
      if (!is.null(brite)) {
        brite <- sapply(brite, function(x) trimws(x)) %>% unique()
      }
      print(selected_terms)
      kegg_gene_pair[[ selected_terms ]] <- brite
    }
  }
  
  
  kegg_gene_pair %>% 
    stack -> selection_query_matrix
  as.character(selection_query_matrix$ind) -> selection_query_matrix$ind
  dcast(selection_query_matrix, ind ~ values, fill = "") -> selection_query_matrix
  row.names(selection_query_matrix) <- selection_query_matrix$ind
  selection_query_matrix$ind <- NULL
  selection_query_matrix != "" -> selection_query_matrix 
  
  return(selection_query_matrix)
}


#
#
kegg_get_wrapper <- function(term, organism, selection){
  query <- tryCatch(
    keggGet(full_terms),
    error=function(e) return(NA))
  
  if (!is.na(query) & !is.null(unlist(query[[1]][selection]))) {
    query[[1]][selection] %>% 
      unlist() %>% 
      unname() -> features
    
    return(features)
  }
  return(NA)
}

# Connetect to KEGG for each gene and create sparse matrix (for grouping)  
# RETURN: Binary sparse matrix for KEGG Orthology per gene
get_KEGG_matrix_from_df <- function(
  assoc_df,
  selection,
  query_column="external_gene_name",
  organism="hsa"
  ) {

  selection_query_matrix <- assoc_df[, query_column] %>% 
    unique() %>%
    {sapply(., function(x) kegg_get_wrapper(x, organism, selection))} %>% 
    stack() %>% 
    dcast(ind ~ values) 
  row_names <- selection_query_matrix[, "ind"]
  selection_query_matrix$ind <- NULL
  selection_query_matrix %>%
    mutate_all(as.logical) -> selection_query_matrix
  row.names(selection_query_matrix) <- row_names
  return(selection_query_matrix)
  
}

# Write to stderr
write_stderr <- function(x){write(x, stderr())}

# Main
args <- commandArgs(TRUE)
parse_args <- function(x) strsplit(
  sub("^--", "", x), "=")

arg_df                <- as.data.frame(
  do.call("rbind", parse_args(args[-1])),
  stringsAsFactors=FALSE)

row.names(arg_df)     <- arg_df[, 1]
assoc_file            <- arg_df["association_file",2]
annotated_output      <- arg_df["annotated_output", 2]
feature_matrix_output <- arg_df["feature_matrix_output",2]

write_stderr(paste0("Reading rsid in: ", assoc_file))
read.csv(assoc_file, sep="\t", stringsAsFactors = F) %>% 
  get_gene_id() %T>%
  write.table(annotated_output, row.names = FALSE, sep = "\t") -> association_df
write_stderr(paste0("Gene annotated rsids can be found at: ", annotated_output))

#
association_df <- association_df[!is.na(association_df$P),]


write_stderr(paste0("Getting KEGG Onthology annotation"))
start <- Sys.time()
association_df %>% 
  get_KEGG_matrix_from_gene(selection = "BRITE") -> selection_query_matrix 
end <- Sys.time()
print(end-start)

write_stderr(paste0("Writing feature matrix to: ", feature_matrix_output))
nr_genes <- nrow(selection_query_matrix)
nr_col <- ncol(selection_query_matrix)
selection_query_matrix %>% 
  {.[, colSums(.)/nr_genes > 0.05 ]} %T>% 
  write.table(feature_matrix_output, row.names = T, sep = "\t") -> selection_query_matrix

write_stderr(paste0(
  "Keeped ",
  ncol(selection_query_matrix)/nr_col *100,
  "% of Onthologies with morethan than 5 % occurence"))


  



