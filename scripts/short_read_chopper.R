suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("--tsv"), type="character", default=NULL, 
              help="path to tsv with simulated mRNAs", metavar="character"),
  make_option(c("--insert_size"), type="character", default="200,300", 
              help="library insert size range, e.g. 200,300", metavar="character"),
  make_option(c("--read_length"), type="numeric", default=100, 
              help="length of each individual read", metavar="numeric"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="Number of threads", metavar="numeric"),
  make_option(c("--seq_depth"), type="numeric", default=20000000, 
              help="sequencing depth", metavar="numeric"),
  make_option(c("--tpm_lower_limit"), type="numeric", default=5, 
              help="TPM lower limit", metavar="numeric"),
  make_option(c("--tpm_upper_limit"), type="numeric", default=200, 
              help="TPM upper limit", metavar="numeric"),
  make_option(c("-o", "--dir_out"), type="character", default=".", 
              help="dir out to create temp files", metavar="character"),
  make_option(c("--seed"), type="numeric", default=NULL,
              help="seed for randomization"),
  make_option(c("--fragments"), type="character", default='no', 
              help="to export ground truth", metavar="character"),
  make_option(c("--sizeselectiontype"), type="character", default="probabilistic",
              help="Type of size selection: none, hardcut, or probabilistic"),
  make_option(c("--no_fragmentation"), action="store_true", default=FALSE,
              help="If specified, do not fragment; keep full molecule length")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tsv)){
  stop("path to tsv with simulated mRNAs is missing", call.=FALSE)
}

if (!is.null(opt$seed)){
  set.seed(opt$seed)
}

setDTthreads(threads = opt$threads)
insert_size <- opt$insert_size
insert_size_split <- strsplit(insert_size, ",")[[1]]
insert_size <- c(as.numeric(insert_size_split[1]), as.numeric(insert_size_split[2]))

file <- opt$tsv

file_importer <- function(file){
  # Use data.table to quickly load, intentionally dropping strings if present
  full_reads <- data.table::fread(file=file, sep = '\t', header = T, stringsAsFactors = FALSE)
  full_reads$transcript_id <- as.numeric(1:nrow(full_reads))
  
  if (!"sequence_length" %in% names(full_reads)) {
      stop("sequence_length column missing. Ensure upstream scripts are updated.")
  }
  
  # Clean up memory immediately if strings are present
  if ("full_molecule_sequence" %in% names(full_reads)) {
      full_reads[, full_molecule_sequence := NULL]
  }
  
  return(full_reads)
}

get_reads <- function(lengths, eta_val = 200, insert_size, transcript_id, no_fragmentation=FALSE, sizeselectiontype="hardcut") {
  
  if (no_fragmentation) {
    fragments <- data.frame(
      transcript_id = transcript_id,
      read_start = 1,
      read_end = lengths,
      length = lengths, 
      size_selection = 'filtered_out' 
    )
    
    if (sizeselectiontype == 'none') {
      if (fragments$length >= 1) fragments$size_selection <- 'passed_size_selection'
    } else if (sizeselectiontype == 'probabilistic') {
      p_left <- 1 / (1 + exp(-0.1 * (fragments$length - insert_size[1])))
      p_right <- 1 / (1 + exp(0.1 * (fragments$length - insert_size[2])))
      prob_keep <- p_left * p_right
      if (runif(1) <= prob_keep) fragments$size_selection <- 'passed_size_selection'
    } else {
      if (fragments$length >= insert_size[1] && fragments$length <= insert_size[2]) {
        fragments$size_selection <- 'passed_size_selection'
      }
    }
    fragments$length <- NULL 
    return(fragments)
  }
  
  eta_val <- mean(c(insert_size[1],insert_size[2]))
  
  deltas <- log10(lengths)
  ns_minus_1 <- pmax(round(lengths / eta_val / gamma(1 / deltas + 1)) - 1, 0)
  xis <- lapply(ns_minus_1, function(n) { diff(sort(c(runif(n), 0, 1))) })
  xis_transformed <- mapply(function(x, d) { x^(1 / d) }, xis, deltas, SIMPLIFY = FALSE)
  delta_is <- mapply(function(len, x_t) {
    round(len * x_t / sum(x_t))
  }, lengths, xis_transformed, SIMPLIFY = FALSE)
  
  starts <- lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(sample(min(insert_size[1], d[1]), 1), cumsum(d[1:(length(d) - 1)]))
    } else {
      sample(min(insert_size[1], d), 1)
    }
  })
  
  ends <- lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(cumsum(d[1:(length(d) - 1)]), sum(d) - sample(min(insert_size[1], sum(d) - d[length(d)]), 1))
    } else {
      d
    }
  })
  
  fragments <- data.frame(
    transcript_id = rep(transcript_id, length(delta_is[[1]])),
    read_start = unlist(starts),
    read_end = unlist(ends)
  )
  
  first_row <- data.frame(transcript_id = fragments$transcript_id[1], read_start = 1, read_end = fragments$read_start[1])
  last_row <- data.frame(transcript_id = fragments$transcript_id[nrow(fragments)], read_start = fragments$read_end[nrow(fragments)], read_end = lengths)
  
  fragments <- rbind(first_row, fragments, last_row)
  fragments$length <- fragments$read_end - fragments$read_start
  
  fragments$size_selection <- 'filtered_out'
  
  if (sizeselectiontype == 'none') {
    fragments[fragments$length >= 1, 'size_selection'] <- 'passed_size_selection'
  } else if (sizeselectiontype == 'probabilistic') {
    p_left <- 1 / (1 + exp(-0.1 * (fragments$length - insert_size[1])))
    p_right <- 1 / (1 + exp(0.1 * (fragments$length - insert_size[2])))
    prob_keep <- p_left * p_right
    keep_roll <- runif(nrow(fragments))
    fragments[keep_roll <= prob_keep, 'size_selection'] <- 'passed_size_selection'
  } else {
    fragments[fragments$length >= insert_size[1] & fragments$length <= insert_size[2], 'size_selection'] <- 'passed_size_selection'
  }
  
  fragments$length <- NULL
  return(fragments)
}

# --- Fast Vectorized Coverage Function ---
calculate_coverage_fast <- function(starts, ends, multiplier = 1) {
  max_pos <- max(ends)
  change_vec <- numeric(max_pos + 2)
  
  start_counts <- table(starts)
  end_counts <- table(ends + 1)
  
  change_vec[as.integer(names(start_counts))] <- as.integer(start_counts)
  change_vec[as.integer(names(end_counts))] <- change_vec[as.integer(names(end_counts))] - as.integer(end_counts)
  
  freq_vec <- cumsum(change_vec)[1:max_pos]
  
  freq_table <- data.table(position = 1:max_pos, frequency = freq_vec * multiplier)
  return(freq_table[frequency > 0])
}
# -----------------------------------------

full_reads <- file_importer(file)
gene_id <- sub(".*/(ENSG[0-9]+(?:_background)?)(?:\\.tsv\\.gz)$", "\\1", file)

list_with_all_fragments <- lapply(seq_len(nrow(full_reads)), function(i) {
  get_reads(lengths = full_reads$sequence_length[i],
            insert_size = insert_size,
            transcript_id = full_reads$transcript_id[i],
            no_fragmentation = opt$no_fragmentation,
            sizeselectiontype = opt$sizeselectiontype)
})

reads_list <- dplyr::bind_rows(list_with_all_fragments, .id = 'transcript_id')
reads_list_all_fragments <- reads_list
reads_list <- reads_list[reads_list$size_selection == 'passed_size_selection', ]

if (nrow(reads_list) > 10){ 
  reads_list$transcript_id <- as.numeric(reads_list$transcript_id)
  gene_length <- mean(full_reads$sequence_length)/1000
  gene_tpm <- runif(n=1, min=as.integer(opt$tpm_lower_limit), max=as.integer(opt$tpm_upper_limit))
  seq_depth <- opt$seq_depth / 10^6
  
  reads_to_get <- gene_length * gene_tpm * seq_depth
  
  reads_list$read_start <- as.integer(reads_list$read_start)
  reads_list$read_end <- as.integer(reads_list$read_end)
  
  # --- Post-Selection Ground Truth ---
  if (opt$fragments == 'with_ground_truth'){
    # Use fast vectorized math instead of sequence expansion
    freq_table_selected <- calculate_coverage_fast(reads_list$read_start, reads_list$read_end)
    
    dir_ground_truth <- paste(opt$dir_out, '/ground_truth_after_size_selection', sep = '')
    dir.create(dir_ground_truth, showWarnings = FALSE, recursive = TRUE)
    path_for_file <- paste(dir_ground_truth, '/', gene_id, '.tsv.gz', sep = '')
    fwrite(freq_table_selected, path_for_file, sep = '\t', col.names = T, row.names = F, quote = F, compress = 'gzip')
  }
  
  reads_list_collapsed <- reads_list %>%
    dplyr::mutate(read_coordinate = paste0(read_start, "-", read_end)) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(read_coordinates = paste(read_coordinate, collapse = ","), .groups = "drop")
  
  # --- Pre-Selection Ground Truth ---
  if (opt$fragments == 'with_ground_truth'){
    
    n_passed <- nrow(reads_list_all_fragments[reads_list_all_fragments$size_selection == 'passed_size_selection',])
    if(n_passed == 0) {
      scaling_factor <- 0
    } else {
      # Instead of looping 100 times to approximate expectation, scale mathematically:
      # (Target fragments / Total available pre-selection fragments)
      scaling_factor <- reads_to_get / n_passed
    }
    
    final_freqs <- calculate_coverage_fast(reads_list_all_fragments$read_start, reads_list_all_fragments$read_end, multiplier = scaling_factor)
    
    dir_ground_truth <- paste(opt$dir_out, '/ground_truth_pre_size_selection', sep = '')
    dir.create(dir_ground_truth, showWarnings = FALSE, recursive = TRUE)
    path_for_file <- paste(dir_ground_truth, '/', gene_id, '.tsv.gz', sep = '')
    fwrite(final_freqs, path_for_file, sep = '\t', col.names = T, row.names = F, quote = F, compress = 'gzip')
  }
  
  temp_dir <- paste(opt$dir_out, '/temp/mRNAs_with_fragments', sep = '')
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  path_for_file <- paste(temp_dir, '/', gene_id, '_fragments.tsv', sep = '')
  
  fwrite(reads_list_collapsed, path_for_file, sep = '\t', col.names = T, row.names = F, quote = F)
}
