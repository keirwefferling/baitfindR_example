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
