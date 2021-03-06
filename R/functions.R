# Functions for use in main plan ----

# When run, this function will concatenate all files in 02_clustering
# ending with 'allbyall.blast.outfmt6' into a single file called 'all.rawblast'
cat_cdhitest_blast_files <- function (...) {
  system("cat 02_clustering/*.allbyall.blast.outfmt6 > 02_clustering/all.rawblast")
}

# Split a vector into a list of n chunks
# Helper function for sort_clusters_into_subfolders.
# 144k views on SO and counting
# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# Given a large number of clusters (fasta files) in `main_dir`,
# copy them into subfolders within in `subfolders_dir`. Set the number
# subfolders with `num_subfolders`. If `overwrite` is set to TRUE, ALL content
# in `subfolders_dir` will be erased first.
sort_clusters_into_subfolders <- function (main_dir, subfolders_dir, num_subfolders, overwrite = FALSE, ...) {
  
  assertthat::assert_that(assertthat::is.dir(main_dir))
  
  assertthat::assert_that(assertthat::is.dir(subfolders_dir))
  
  assertthat::assert_that(assertthat::is.number(num_subfolders))
  
  assertthat::assert_that(is.logical(overwrite))
  
  # Delete old output if overwrite is TRUE
  old_output <- list.files(subfolders_dir)
  if(isTRUE(overwrite) & length(old_output) > 0) {
    fs::dir_delete(subfolders_dir)
    fs::dir_create(subfolders_dir)
  }
  
  # Create subfolders in specified directory
  fs::dir_create(fs::path(subfolders_dir, 1:num_subfolders))
  
  # Get list of fasta files in main directory, split into groups
  fastas_to_copy <- list.files(main_dir, pattern = "*.fa$", full.names = TRUE) %>%
    chunk2(num_subfolders)
  
  # Copy fastas into subfolders
  
  purrr::walk2(
    fastas_to_copy, 1:num_subfolders,
    ~fs::file_copy(.x, fs::path(subfolders_dir, .y), overwrite = TRUE)
  )
  
  # Return info about the subfolders directory for tracking
  fs::dir_info(subfolders_dir)
  
}

# Copy all of the files in subfolders to a single main folder.
# `subfolders_dir` contains all the subfolders.
aggregate_files_in_subfolders <- function (subfolders_dir, main_dir, ...) {
  
  assertthat::assert_that(assertthat::is.dir(main_dir))
  
  assertthat::assert_that(assertthat::is.dir(subfolders_dir))
  
  files_to_copy <- list.files(subfolders_dir, recursive = TRUE, full.names = TRUE)
  
  purrr::walk(
    files_to_copy,
    ~fs::file_copy(., main_dir, overwrite = TRUE)
  )
  
  # Return info about the main_dir for tracking
  fs::dir_info(main_dir)
}

# Functions for processing example data -----

#' Downsize a file from the top
#' 
#' Only the lines specified from the first to the fraction
#' of the file specified by `keep_frac` will be kept.
#'
#' @param full_file Path to input file
#' @param path Path to write out downsized file
#' @param keep_frac Fraction of file to retain
#'
#' @return NULL
downsize_simple <- function(full_file, path, keep_frac) {
  
  full <- readr::read_lines(full_file)
  short <- full[1:round(length(full)*keep_frac, 0)]
  readr::write_lines(short, path)
  
}

# Make a table of oneKp sample codes and URLs to download transcriptome assemblies.
make_download_table <- function(url = "http://www.onekp.com/public_data.html") {
  
  # Scrape oneKp links for assemblies.
  # Use Selectorgadget to identify the CSS component with links.
  # Note: run vignette(“selectorgadget”) for a refresher on how to do this.
  links <- xml2::read_html(url) %>% rvest::html_nodes("td:nth-child(5) a")
  
  # Next, read in the main table (just plain text, without links).
  tbls_xml <- XML::readHTMLTable(url)
  
  # Extract the table and add links to assemblies,
  # keeping only transcriptomes in baitfindR example dataset .
  tbls_xml[[1]] %>% 
    tibble::as_tibble() %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate(assembly_link = map_chr(links, xml2::xml_attrs)) %>%
    dplyr::select(code = `1kP_Code`, assembly_link)
}

# Downsize a transcriptome file
downsize_transcriptome <- function (file, keep_frac) {
  # Read in file. Use bzfile() because it's compressed
  seq <- ape::read.FASTA(bzfile(file))
  # Randomly downsize to specified fraction of transcripts
  seq <- seq[sample(1:length(seq), keep_frac*length(seq))]
}

# Untar with tracking
untar_tracked <- function (tarfile, exdir, ...) {
  untar(tarfile = tarfile, exdir = exdir)
}

# For some reason there are a bunch of duplicate names with
# different sequences in the lygodium data. Exclude these.
trim_proteome <- function (proteome) {
  dup_names <- names(proteome)[duplicated(names(proteome))]
  trimmed_proteome <- proteome[!(names(proteome) %in% dup_names)]
  return(trimmed_proteome)
}

# Functions to help standardize names in plan 
get_gff3 <- function (genome) {
  switch(genome,
         arabidopsis = here::here("data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz"),
         azolla = here::here("data_raw/Azolla_filiculoides.gene_models.highconfidence_v1.1.gff"),
         salvinia = here::here("data_raw/Salvinia_cucullata.gene_models.highconfidence_v1.2.gff"))
}

get_genome <- function (genome) {
  switch(genome,
         arabidopsis = here::here("data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta"),
         azolla = here::here("data_raw/Azolla_filiculoides.genome_v1.2.fasta"),
         salvinia = here::here("data_raw/Salvinia_cucullata.genome_v1.2.fasta")
  )
}

get_sources <- function (genome) {
  switch(genome,
         arabidopsis = "araport11",
         azolla = "maker",
         salvinia = "maker"
  )
}

# Rename arabidopsis genomes so that 
# names match between the fasta file
# and the bed file
rename_arabidopsis_genome <- function (fasta_in, fasta_out, ...) {
  arabidopsis_genome <- read.FASTA(fasta_in)
  
  new_names <- case_when(
    grepl("1 dna", names(arabidopsis_genome)) ~ "1",
    grepl("2 dna", names(arabidopsis_genome)) ~ "2",
    grepl("3 dna", names(arabidopsis_genome)) ~ "3",
    grepl("4 dna", names(arabidopsis_genome)) ~ "4",
    grepl("5 dna", names(arabidopsis_genome)) ~ "5",
    grepl("Mt dna", names(arabidopsis_genome)) ~ "Mt",
    grepl("Pt dna", names(arabidopsis_genome)) ~ "Pt")
  
  if(anyNA(new_names)) {stop("Replacement names not provided for all fasta headers")}
  
  names(arabidopsis_genome) <- new_names
  
  write.FASTA(arabidopsis_genome, fasta_out)
}

# Additional functions ----

#' Write out a list of DNA sequences in chunks
#' 
#' The list will be split into a list of lists, to
#' save memory when writing out sequences.
#'
#' @param fasta_list The list of DNA sequences
#' @param chunk_size Size of chunks to split up list into.
#' @param out_dir Directory to write to
#'
#' @return List of hashes; digest of each chunk. Externally, the 
#' sequences will be written to out_dir
#' 
write_fasta_files_chunked <- function(fasta_list, chunk_size = 100, out_dir) {
# Split list of fasta files to write out into chunks
n <- length(fasta_list)
r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
fasta_list_split <- split(fasta_list, r)

# Write out in chunks
purrr::map(fasta_list_split, ~baitfindR::write_fasta_files(., out_dir = out_dir))
}

#' Helper function to subset an alignment
#'
#' The alignment will be subset to the bases between the start 
#' and stop positions, inclusive
#'
#' @param alignment Object of class "DNAbin"; input alignment
#' @param start Integer; start position
#' @param stop Integer; stop position
#'
#' @return Object of class "DNAbin"; subset alignment
#' 
subset_alignment <- function (alignment, start, stop) {
  assertthat::assert_that(assertthat::is.number(start))
  assertthat::assert_that(assertthat::is.number(stop))
  assertthat::assert_that(is.matrix(alignment))
  assertthat::assert_that(inherits(alignment, "DNAbin"))
  assertthat::assert_that(start <= stop)
  assertthat::assert_that(stop >= start)
  assertthat::assert_that(start >= 1)
  assertthat::assert_that(stop <= ncol(alignment))
  
  alignment[,start:stop]
}

#' Extract exons from an alignment
#' 
#' Introns are assumed to be coded as "n" for all sequences at a given site.
#'
#' @param alignment 
#'
#' @return Dataframe with exons as a list-column
#' 
extract_exons <- function (alignment) {
  # Get the column positions of all exons
  exon_cols <- apply(as.character(alignment), 2, function(x) all(x =="n")) %>%
    magrittr::not() %>%
    which
  
  tibble(exon_cols = exon_cols) %>%
    # Group the exons by consecutive column positions
    mutate(exon = cumsum(c(1, diff(exon_cols) != 1))) %>%
    group_by(exon) %>%
    # Find the start and stop of each exon
    summarize(
      exon_start = min(exon_cols),
      exon_stop = max(exon_cols)
    ) %>%
    ungroup %>%
    # Extract each exon from the alignment
    mutate(exon_alignment = map2(.x = exon_start, .y = exon_stop, ~subset_alignment(alignment, .x, .y))) %>%
    select(exon, exon_alignment)
}

#' Convert a dataframe of combined alignments to exons
#'
#' @param combined_alignments_data Dataframe; candidate baits (alignments),
#' with columns `bait_id` and `alignment`. Each row is one gene, including
#' both introns and exons (introns are regions where all sequences are "n").
#'
#' @return Dataframe, with each row an exon from the input alignments.
#' 
combined_alignments_to_exons <- function (combined_alignments_data) {
  
  combined_alignments_data %>%
    # Extract exons from each alignment
    select(bait_id, alignment) %>%
    mutate(exons = map(alignment, extract_exons)) %>%
    select(bait_id, exons) %>%
    # Now "alignment" referers to each exon
    unnest(exons) %>%
    rename(alignment = exon_alignment) %>%
    # Remove empty sequences from alignments
    mutate(alignment = map(alignment, exclude_empty_taxa)) %>%
    # Calculate summary stats for each exon
    mutate(
      length = purrr::map_dbl(alignment, ncol),
      ntaxa = purrr::map_dbl(alignment, nrow),
      GC_content = purrr::map_dbl(alignment, GC.content),
      mean_dist = purrr::map_dbl(alignment, ~ape::dist.dna(., model = "raw", pairwise.deletion = TRUE) %>% mean(na.rm = TRUE)),
      max_dist = purrr::map_dbl(alignment, ~ape::dist.dna(., model = "raw", pairwise.deletion = TRUE) %>% max(na.rm = TRUE)),
      pars_inf = purrr::map_dbl(alignment, ~ips::pis(., "fraction")))
  
}

#' Calculate base frequences *per sample* in an alignment
#'
#' @param align DNA alignment (list of class "DNAbin")
#'
#' @return Tibble with frequencies of each base as columns and samples as rows
base_freq_per_taxon <- function (align) {
  
  purrr::map_df(
    1:nrow(align), 
    ~ape::base.freq(align[.,], freq = TRUE, all = TRUE) %>% 
      dplyr::bind_rows()) %>%
    dplyr::rename(gap = `-`, missing = `?`) %>%
    dplyr::mutate(sample = rownames(align))
}

#' Exclude empty sequences from an alignment
#' 
#' (Empty sequences are those with zero A, C, T, or G)
#'
#' @param align DNA alignment (list of class "DNAbin")
#'
#' @return DNA alignment (list of class "DNAbin") with empty sequences removed
#' 
exclude_empty_taxa <- function (align) {
  
  to_keep <- 
    # Calculate base frequences for each sample (row of alignment)
    base_freq_per_taxon(align) %>%
    # Exclude sample if frequences of all actg bases are zero
    dplyr::mutate(to_exclude = ifelse(a == 0 & c == 0 & g == 0 & t == 0, TRUE, FALSE)) %>%
    dplyr::filter(to_exclude == FALSE) %>%
    dplyr::pull(sample)
  
  # Subset to only samples flagged to keep
  align[rownames(align) %in% to_keep,]
  
}

#' Filter a set of exons to a bait set of a specified size
#' 
#' After filtering by `min_exon_length`, `max_exon_length`, `min_gc`, `max_gc` and `min_ntaxa`,
#' exons will be arranged in decreasing order of parsimony informativeness (PI), and
#' one exon per gene with the greatest PI will be selected. The cumulative number of
#' baits will then be calculated, and the final set of baits truncated so it does not
#' exceed `bait_kit_size`.
#'
#' @param exon_alignments Dataframe; exons in a tibble. Output of combined_alignments_to_exons()
#' @param bait_length Length of the synthesized baits in bp (typically 105 or 120)
#' @param bait_coverage Mean depth of tiled baits (typically in the range of 1-3)
#' @param bait_kit_size Total number of baits in the kit (typically in the range of 20000-60000)
#' @param min_exon_length Minimum length of exons to include in the baits; 
#' exons shorter than this will be excluded
#' @param max_exon_length Maximum length of exons to include in the baits; 
#' exons longer than this will be excluded
#' @param min_gc Minimum GC content allowed to include in the baits; 
#' exons with lower % GC than this will be excluded
#' @param max_gc Maximum GC content allowed to include in the baits; 
#' exons with higher % GC than this will be excluded
#' @param min_ntaxa Minimum number of taxa required in each exon alignment to be included in the baits;
#' exon alignments with fewer taxa will be excluded
#'
#' @return Tibble
#' 
filter_baits <- function(
  exon_alignments, 
  bait_length = 105, 
  bait_coverage = 2.5, 
  bait_kit_size = 20000,
  min_exon_length = 120,
  max_exon_length = 300,
  min_gc = 0.25,
  max_gc = 0.55,
  min_ntaxa) {
  
  exon_alignments %>%
    # Filter exon alignments by length, number of taxa, and GC content
    dplyr::filter(
      length > min_exon_length, 
      length < max_exon_length, 
      GC_content > min_gc,
      GC_content < max_gc,
      ntaxa > min_ntaxa) %>%
    # Select one exon per gene with the greatest pars. inform.
    dplyr::group_by(bait_id) %>%
    dplyr::arrange(dplyr::desc(pars_inf)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # Arrange by descending parsimony informativeness
    dplyr::arrange(dplyr::desc(pars_inf)) %>%
    # Calculate number of baits required
    dplyr::mutate(
      # - Calculate total number of non-missing bases in each alignment
      total_bp = purrr::map_dbl(alignment, ~ape::base.freq(x = ., all = FALSE, freq = TRUE) %>% sum),
      # - Calculate the number of baits needed to cover each alignment
      num_baits = (total_bp / bait_length) * bait_coverage,
      # - Calculate running total of the number of baits needed as the number of 
      # included exons increases
      num_baits_cum = cumsum(num_baits)
    ) %>%
    # Cap the number of baits so it does not exceed the bait kit size
    dplyr::filter(num_baits_cum < bait_kit_size) %>%
    dplyr::select(-num_baits_cum)
}

