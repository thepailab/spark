#!/usr/bin/env Rscript
# Author: Shaimae Elhajjajy
# Modified by Ezequiel Calvo-Roitberg
# Purpose: This R script will analyze the GENCODE annotations for the most highly expressed genes in a cell line the goal is to curate a robust, diverse, and representative gene
# set to be used in simulations.
### SET UP ###
# Import libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(parallel))
#Parse arguments
#####
option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="genome annotation", metavar="character"),
  make_option(c("-f", "--genome_fasta"), type="character", default=NULL,
              help="genome fasta", metavar="character"),
  make_option(c("-o", "--dir_out"), type="character", default=".",
              help="dir out for the bed and metric files", metavar="character"),
  make_option(c("--protein_coding_only"), action="store_true", default=FALSE,
              help="only use protein-coding genes"),
  make_option(c("-t", "--threads"), type="numeric", default="1",
              help="Number of threads", metavar="numeric"),
  make_option(c("-c", "--n_gene_clusters"), type="numeric", default="6",
              help="Number of gene clusters to get", metavar="numeric"),
  make_option(c("-n", "--number_of_genes"), type="numeric", default="120",
              help="Number of total genes to simulate", metavar="numeric"),
  make_option(c("--seed"), type="numeric", default=NULL,
              help="seed for randomization")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$gtf)){
  stop("GTF file is missing", call.=FALSE)
}
if (is.null(opt$genome_fasta)){
  stop("Genome fasta file is missing", call.=FALSE)
}
if (!is.null(opt$seed)){
  set.seed(opt$seed)
}
#Set number of threads for data table
data.table::setDTthreads(threads = opt$threads)
#Make directory to temp files
temp_dir <- paste(opt$dir_out,'/temp',sep = '')
system(paste("mkdir ", temp_dir,sep = ''))
# Read in the gencode annotation data
gencode_data = fread(opt$gtf, sep = "\t",header = F)
### FUNCTIONS ###
# Extract all GENCODE gene IDs robustly
getGencodeGeneIDs <- function(gencode_data){
  gene_ids <- sub('.*gene_id "([^"]+)".*', '\\1', gencode_data$V9)
  gene_ids[!grepl('gene_id "', gencode_data$V9)] <- NA
  return(gene_ids)
}
# Extract all GENCODE Transcript IDs robustly
getGencodeTranscriptIDs <- function(gencode_data){
  transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', gencode_data$V9)
  transcript_id[!grepl('transcript_id "', gencode_data$V9)] <- NA
  return(transcript_id)
}
# Find the type of gene securely, handling custom genes without biotypes
getGeneType <- function(gencode_data){
  gene_type <- ifelse(grepl('gene_type "', gencode_data$V9),
                      sub('.*gene_type "([^"]+)".*', '\\1', gencode_data$V9),
                      ifelse(grepl('gene_biotype "', gencode_data$V9),
                             sub('.*gene_biotype "([^"]+)".*', '\\1', gencode_data$V9),
                             "unknown"))
  return(gene_type)
}
# Find the longest transcript for the gene
getLongestTranscript <- function(gencode_data){
  longest_txpt <- gencode_data[gencode_data$V3=='transcript',] %>% dplyr::group_by(gene_id) %>% dplyr::mutate(transcript_length=V5-V4) %>% dplyr::filter(transcript_length==max(transcript_length)) %>% dplyr::slice(1)#Some genes will have transcripts from the same length, this ensures there's only one
  return(longest_txpt)

}
getfeature_metrics <- function(gencode_data,longest_txpt){
  all_features <- gencode_data[gencode_data$transcript_id%in%longest_txpt$transcript_id,] #Filter to make it faster
  all_exons <- all_features[all_features$V3=='exon',]
  all_exons$exon_length <- all_exons$V5-all_exons$V4
  #Split by strand and sort (doesn't have to be per gene)
  exons_pos_strand <- all_exons[all_exons$V7=='+',]
  exons_neg_strand <- all_exons[all_exons$V7=='-',]
  #Sort
  exons_pos_strand <- exons_pos_strand[order(exons_pos_strand$V4,decreasing = F),]
  exons_neg_strand <- exons_neg_strand[order(exons_neg_strand$V5,decreasing = T),]

  #Assign genomic order
  exons_pos_strand <- exons_pos_strand %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(genomic_order=dplyr::row_number())

  exons_neg_strand <- exons_neg_strand %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(genomic_order=dplyr::row_number())

  all_exons <- rbind(exons_pos_strand,exons_neg_strand)

  exon_coordinates <- data.frame(chr=all_exons$V1,start=all_exons$V4,end=all_exons$V5,id=all_exons$transcript_id,score=all_exons$genomic_order,strand=all_exons$V7)

  #Get exons and introns metrics
  exon_intron_metrics <- all_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(intron_length = dplyr::case_when(
      V7 == "+" ~ dplyr::lead(V4) - V5,
      V7 == "-" ~ V4 - dplyr::lead(V5),
      TRUE ~ NA_real_
    )) %>%
    dplyr::reframe(
      number_of_exons=dplyr::n(),
      number_of_introns=number_of_exons-1,
      mean_exon_length = mean(exon_length, na.rm = TRUE),
      mean_intron_length = mean(intron_length, na.rm = TRUE),
      first_exon_length = dplyr::first(exon_length),
      first_intron_length = dplyr::first(intron_length)
    )

  #Get intron coordinates
  intron_coordinates_raw <- all_exons[,c(1,4,5,14,11,7)] %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      V4_new = dplyr::case_when(
        V7 == "+" ~ V5 + 1,
        V7 == "-" ~ dplyr::lead(V5) + 1,
        TRUE ~ NA_real_
      ),
      V5_new = dplyr::case_when(
        V7 == "+" ~ dplyr::lead(V4) - 1,
        V7 == "-" ~ V4 - 1,
        TRUE ~ NA_real_
      )) %>%
    dplyr::filter(!is.na(V4_new) & !is.na(V5_new))

  intron_coordinates <- data.frame(chr=intron_coordinates_raw$V1,start=intron_coordinates_raw$V4_new,end=intron_coordinates_raw$V5_new,id=intron_coordinates_raw$transcript_id,score=intron_coordinates_raw$genomic_order,strand=intron_coordinates_raw$V7)

  #Get UTR metrics
  utrs <- all_features[all_features$V3%in%c("five_prime_utr","three_prime_utr"),]
  utrs$utr_length <- utrs$V5-utrs$V4
  #GTFs normally might have different features (rows) for utrs if they across multiple exons
  utrs <- utrs %>% dplyr::group_by(transcript_id,V3) %>% dplyr::reframe(utr_length=sum(utr_length))
  utrs_pivot <- tidyr::pivot_wider(utrs[,c('transcript_id','utr_length','V3')],id_cols = transcript_id,values_from = utr_length,names_from = V3)

  return(list(exon_intron_metrics=exon_intron_metrics,utrs_pivot=utrs_pivot,exon_coordinates=exon_coordinates,intron_coordinates=intron_coordinates))
}
getLongestTxptSequence <- function(longest_txpt,exon_coordinates,intron_coordinates,genome_fasta,temp_dir){

  coordinates_transcripts_for_bedtools <- data.frame(chr=longest_txpt$V1,start=longest_txpt$V4,end=longest_txpt$V5,id=paste(longest_txpt$transcript_id,'transcript',sep = '_'),score=longest_txpt$V6,strand=longest_txpt$V7)

  exon_coordinates$id <- paste(exon_coordinates$id,'exon',exon_coordinates$score,sep = '_')
  if (nrow(intron_coordinates)>0){
    intron_coordinates$id <- paste(intron_coordinates$id,'intron',intron_coordinates$score,sep = '_')}

  #Merge all coordinates to run bedtools only once (parsing happens downstream)
  all_coordinates <- rbind(coordinates_transcripts_for_bedtools,exon_coordinates,intron_coordinates)
  all_coordinates$start <- as.integer(all_coordinates$start)
  all_coordinates$end <- as.integer(all_coordinates$end)

  all_coordinates$start <- all_coordinates$start-1

  path_for_coordinates <- paste(temp_dir,"/longest_transcript_and_features.bed",sep = '')
  fwrite(all_coordinates,path_for_coordinates, sep = "\t", append = FALSE, col.names = FALSE,row.names = F)

  #Run nucleotide metrics
  get_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)

    if (length(file_arg) > 0) {
      # If running via Rscript
      return(normalizePath(sub("^--file=", "", file_arg)))
    } else if (interactive()) {
      # If running interactively (e.g., RStudio)
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        return(rstudioapi::getSourceEditorContext()$path)
      } else {
        stop("Cannot determine script path in interactive mode without rstudioapi.")
      }
    } else {
      stop("Script path could not be determined.")
    }
  }

  script_path <-dirname(get_script_path())

  system(paste("python ", script_path,"/nucleotide_metrics.py ", path_for_coordinates,' ',genome_fasta,' ',temp_dir,sep = ''))


  #Read nucleotide metrics
  nucleotide_metrics_path <- paste(temp_dir,'/nucleotide_metrics.tsv',sep='')
  nucleotide_metrics <- fread(nucleotide_metrics_path,sep = '\t',header = T)


  nucleotide_metrics$clean_name <- sub("::.*", "", nucleotide_metrics$Feature)
  nucleotide_metrics$clean_name <- sub("\\(.*\\)$", "", nucleotide_metrics$clean_name)


  nucleotide_metrics$transcript_id <- sub("_(transcript|exon|intron)(_\\d+)?$", "", nucleotide_metrics$clean_name)

  nucleotide_metrics$feature <- sub(".*_(transcript|exon|intron)(_\\d+)?$", "\\1", nucleotide_metrics$clean_name)

  nucleotide_metrics <- nucleotide_metrics[nucleotide_metrics$feature %in% c("transcript", "exon", "intron"), ]

  nucleotide_metrics_pivot <- nucleotide_metrics %>%
    dplyr::select(-Feature, -clean_name) %>%
    tidyr::pivot_longer(cols = A:G, names_to = "base", values_to = "value") %>%
    dplyr::group_by(transcript_id,feature,base) %>%
    dplyr::reframe(new_value=mean(value)) %>%
    tidyr::pivot_wider(names_from = c(feature, base), values_from = new_value)

  #handle GTFs entries with no introns
  expected_nuc_cols <- c("transcript_A", "transcript_C", "transcript_T", "transcript_G",
                         "exon_A", "exon_C", "exon_T", "exon_G",
                         "intron_A", "intron_C", "intron_T", "intron_G")

  missing_cols <- setdiff(expected_nuc_cols, colnames(nucleotide_metrics_pivot))
  if(length(missing_cols) > 0) {
    nucleotide_metrics_pivot[, missing_cols] <- NA
  }

  nucleotide_metrics_pivot <- nucleotide_metrics_pivot[, c("transcript_id", expected_nuc_cols)]

  return(nucleotide_metrics_pivot)

}
#Wrangle to make parsing more efficient
gencode_data$gene_id <- getGencodeGeneIDs(gencode_data)
gencode_data$transcript_id <- getGencodeTranscriptIDs(gencode_data)
gencode_data$gene_type  <- getGeneType(gencode_data)
#Optional protein coding filter
if (opt$protein_coding_only) {
  gencode_data <- gencode_data[gencode_data$gene_type=='protein_coding',]
}
#These are modular metrics instead of a unique function to make it more customizable
#Get gene length
#gene_length_ref <- gencode_data[gencode_data$V3=='gene',] %>% dplyr::group_by(gene_id) %>% dplyr::reframe(gene_length=V5-V4)
gene_rows <- gencode_data[gencode_data$V3 == 'gene', ]
if (nrow(gene_rows) > 0) {
  gene_length_ref <- gene_rows %>%
    dplyr::group_by(gene_id) %>%
    dplyr::reframe(gene_length = V5 - V4)
} else {
  gene_length_ref <- gencode_data %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::reframe(gene_length = max(V5) - min(V4))
}

#Find the longest transcript for the gene
longest_txpt <- getLongestTranscript(gencode_data)
#Add longest transcript info to gene ids to generate the final df with metrics
db_gene_transcript <- gencode_data[gencode_data$transcript_id%in%longest_txpt$transcript_id,c('gene_id','transcript_id','V3','gene_type')]
db_gene_transcript <- db_gene_transcript[db_gene_transcript$V3=='transcript',]
gene_length_ref <- merge(gene_length_ref,db_gene_transcript[,-'V3'],by='gene_id')
#Get the number of exons and introns per transcript
feature_metrics <- getfeature_metrics(gencode_data,longest_txpt)
exon_introns_metrics <- feature_metrics[["exon_intron_metrics"]]
#Get length of UTRs
utr_metrics <- feature_metrics[["utrs_pivot"]]
utr_gene_has <- colnames(utr_metrics)[-1]
exon_coordinates <- feature_metrics[["exon_coordinates"]]
intron_coordinates <- feature_metrics[["intron_coordinates"]]

#Get longest transcript sequence and metrics then parse
#Use cached metrics file if available to skip the expensive bedtools call
path_metrics <- paste(opt$dir_out,'/metrics_per_transcript_df.tsv',sep='')
if (file.exists(path_metrics)) {
  cat("[Info] Found cached metrics_per_transcript_df.tsv — skipping bedtools/nucleotide computation.\n")
  metric_df <- fread(path_metrics, sep='\t', header=T)
} else {
  nucleotide_metrics_pivot <- getLongestTxptSequence(longest_txpt,exon_coordinates,intron_coordinates,opt$genome_fasta,temp_dir)
  #Build the output
  metric_df <- merge(merge(merge(gene_length_ref,merge(longest_txpt[,c(1,4,5,11,6,7,13)],exon_introns_metrics,by='transcript_id'),by='transcript_id'),utr_metrics,by='transcript_id',all.x=T),nucleotide_metrics_pivot,by='transcript_id')
  colnames(metric_df) = c('transcript_id','gene_id','gene_length','gene_type','chr','start','end','score','strand','transcript_length','number_of_exons','number_of_introns','mean_exon_length','mean_intron_length','first_exon_length', 'first_intron_length', utr_gene_has, 'transcript_A', 'transcript_C', 'transcript_T','transcript_G', 'exon_A', 'exon_C', 'exon_T', 'exon_G', 'intron_A', 'intron_C', 'intron_T','intron_G')
  #Write metrics
  fwrite(metric_df, path_metrics, sep = "\t")
}

######
#HC
data_scaling <- function(metric_df){
  #Create a separate dataframe with only the variables that will be used in the hiearchical clustering
  safe_cols <- c("transcript_id", "gene_length", "transcript_length", "number_of_exons", "number_of_introns",
                 "mean_exon_length", "mean_intron_length", "first_exon_length", "first_intron_length",
                 "five_prime_utr", "three_prime_utr",
                 "transcript_A", "transcript_C", "transcript_T", "transcript_G",
                 "exon_A", "exon_C", "exon_T", "exon_G",
                 "intron_A", "intron_C", "intron_T", "intron_G")
  safe_cols_present <- intersect(safe_cols, colnames(metric_df))
  data_to_scale <- metric_df[, ..safe_cols_present]


  #Replace NAs by 0 for genes that don't have a given feature (for example intronless genes)
  data_to_scale <- replace(data_to_scale, is.na(data_to_scale), 0)
  # Standardize the data by subtracting the mean and dividing by the standard deviation
  data_scaled <- scale(data_to_scale[, 2:ncol(data_to_scale)], center = TRUE, scale = TRUE)

  #For positive scaled values, take the log; for negative scaled values, take the absolute value and then the negative log
  log10_transform <- function(x) {
    sign(x) * log10(pmax(abs(x), .Machine$double.eps))
  }

  data_scaled <- log10_transform(data_scaled)
  data_scaled <- cbind(data_to_scale[,c("transcript_id")],data_scaled)

  return(data_scaled)
}
cluster_making <- function(data_scaled,gene_clusters){
  #Run HC
  rownames_data_matrix <- data_scaled$transcript_id
  data_matrix <- as.matrix(data_scaled[, -c(1)])
  rownames(data_matrix) <- rownames_data_matrix
  data.clust <- hclust(d = dist(x = data_matrix), method = "complete")
  #Get a given number of clusters
  gene_clusters_df = cutree(data.clust, k = gene_clusters)
  #Wrangle
  gene_clusters_df <- as.data.frame(gene_clusters_df)
  gene_clusters_df$transcript_id <- rownames(gene_clusters_df)
  rownames(gene_clusters_df) <- NULL
  return(gene_clusters_df)
}
gene_sampling <- function(gene_clusters_df, number_of_genes, n_gene_clusters) {
  # Calculate the initial target per cluster
  n_genes_per_cluster_to_sample <- floor(number_of_genes / n_gene_clusters)
  n_genes_per_cluster <- table(gene_clusters_df$gene_clusters) %>% as.data.frame()
  colnames(n_genes_per_cluster) <- c('gene_clusters','total_genes')

  #find clusters with fewer genes than the target and calculate leftovers
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      to_sample = pmin(total_genes, n_genes_per_cluster_to_sample),  
      leftover = pmax(0, n_genes_per_cluster_to_sample - total_genes)
    )

  #calculate total leftover genes
  total_leftover <- sum(n_genes_per_cluster$leftover)

  #redistribute leftover genes to clusters that can accommodate more
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      can_take_more = total_genes - to_sample,
      additional = pmin(total_leftover, can_take_more))%>%
    dplyr::mutate(additional=floor(pmin(total_leftover, can_take_more)/length(additional>0)))

  #adjust to_sample based on redistribution
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      to_sample = to_sample + additional
    )

  #sample based on adjusted to_sample per cluster
  sampled_genes <- merge(gene_clusters_df, n_genes_per_cluster, by = "gene_clusters") %>%
    dplyr::group_by(gene_clusters) %>%
    dplyr::group_map(~ dplyr::slice_sample(.x, n = .x$to_sample[1]), .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::select(gene_clusters, transcript_id)

  return(sampled_genes)
}
gtf_generator <- function(path_to_fastatsv,sampled_genes,exon_coordinates,intron_coordinates,db_gene_transcript,path_for_gtf,gene_length_ref,threads=1){

  path_to_tsv <- paste(path_to_fastatsv, 'longest_transcript_and_features.tsv', sep='')
  seq_per_feature <- fread(path_to_tsv, header = F)

  seq_per_feature <- seq_per_feature[!grepl("_transcript(\\(.*\\))?$", seq_per_feature$V1), ]

  #extract transcript_id, feature, and position even if IDs contain underscores ()
  seq_per_feature <- seq_per_feature %>%
    dplyr::mutate(
      transcript_id = sub("_(exon|intron)_\\d+.*$", "", V1),
      feature = sub(".*_(exon|intron)_\\d+.*$", "\\1", V1),
      position = as.numeric(sub(".*_(exon|intron)_(\\d+).*$", "\\2", V1))
    )

  seq_per_feature <- seq_per_feature[, c("transcript_id", "feature", "position", "V2")]
  colnames(seq_per_feature)[4] <- 'sequence'

  #filter
  seq_per_feature <- seq_per_feature[seq_per_feature$transcript_id%in%sampled_genes$transcript_id,]

  #add coordinate info per feature
  exon_coordinates$feature <- 'exon'
  if (nrow(intron_coordinates)>0){intron_coordinates$feature <- 'intron'}

  all_features <- rbind(exon_coordinates,intron_coordinates)
  colnames(all_features)[c(4,5)] <- c('transcript_id','position')
  #add gene id for export
  all_features <- merge(all_features,db_gene_transcript[,c('transcript_id','gene_id')],by='transcript_id')

  seq_per_feature_full <- merge(seq_per_feature,all_features,by=c('transcript_id','feature','position'))

  #include PAS coords
  seq_per_feature_full <- merge(seq_per_feature_full, gene_length_ref[, c('gene_id', 'gene_length')], by='gene_id')
  seq_per_feature_full$pas_coordinate <- seq_per_feature_full$gene_length

  seq_per_feature_full_for_export <- seq_per_feature_full[,c('chr','start','end','gene_id','feature','position','strand','sequence','pas_coordinate')]

  #Sort
  seq_per_feature_full_for_export_pos <- seq_per_feature_full_for_export[seq_per_feature_full_for_export$strand=='+',]
  seq_per_feature_full_for_export_neg <- seq_per_feature_full_for_export[seq_per_feature_full_for_export$strand=='-',]
  seq_per_feature_full_for_export_pos <- seq_per_feature_full_for_export_pos[order(seq_per_feature_full_for_export_pos$start,decreasing = F),]
  seq_per_feature_full_for_export_neg <- seq_per_feature_full_for_export_neg[order(seq_per_feature_full_for_export_neg$end,decreasing = T),]
  seq_per_feature_full_for_export <- rbind(seq_per_feature_full_for_export_pos,seq_per_feature_full_for_export_neg)


  seq_per_feature_full_for_export <- seq_per_feature_full_for_export %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(position = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    as.data.frame()


  path_for_gtf <- paste(path_for_gtf,'gtf/',sep='')
  system(paste('mkdir ',path_for_gtf,sep=''))

  n_cores <- max(1L, as.integer(threads))
  gene_list <- split(seq_per_feature_full_for_export, seq_per_feature_full_for_export$gene_id)
  parallel::mclapply(gene_list, function(gene_df) {
    out_path <- paste0(path_for_gtf, gene_df$gene_id[1], ".tsv.gz")
    data.table::fwrite(gene_df, file=out_path, sep="\t", quote=FALSE, compress="gzip")
  }, mc.cores=n_cores)

  system(paste('rm ',path_to_fastatsv,'longest_transcript_and_features.tsv',sep='')) #Remove temp file

}

metric_df <- fread(path_metrics,sep='\t',header = T)
data_scaled <- data_scaling(metric_df)
if (opt$number_of_genes == 0) {
  #cat("[Info] --n set to 0. Simulating ALL genes found in the GTF...\n")
  sampled_genes <- data_scaled[, "transcript_id", drop = FALSE]

} else {

  if (nrow(data_scaled) >= opt$n_gene_clusters & nrow(data_scaled) > 1){
    gene_clusters_df <- cluster_making(data_scaled, opt$n_gene_clusters)
    colnames(gene_clusters_df)[1] <- 'gene_clusters'
    # Sample clusters
    sampled_genes <- gene_sampling(gene_clusters_df, opt$number_of_genes, opt$n_gene_clusters)
  } else {
    sampled_genes <- dplyr::sample_n(data_scaled, opt$number_of_genes)
  }
}
#############
#GTF for mRNA generation
gtf_generator(paste(temp_dir,'/',sep=''),sampled_genes,exon_coordinates,intron_coordinates,db_gene_transcript,opt$dir_out,gene_length_ref,threads=opt$threads)

