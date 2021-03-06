# baitfindR_simple

# An example project showing how baitfindR can be combined with drake to
# automate the Yang and Smith (2014) (hereafter, Y&S) workflow to find
# candidate baits for sequence capture from a set of transcriptomes and
# reference genomes.

# Set working directory
setwd(here::here())

# Load packages
source("R/packages.R")

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Update baitfindR settigs
baitfindR::set_ys_path("/home/phylogenomic_dataset_construction")

# Load functions
source("R/functions.R")

### Define basic input values

# Set seed for reproducibility
set.seed(5394)

# Vector of transcriptome samples. These are all cheilanthoid ferns.
codes <- c("ACAS", "ADCA", "ADCS", "ADJO", "ALCS", "ANSS", "BOHS", "CHCS", "COFS", "GAAG", "HAAS", "HAHE", "POPJ", "SKYV", "UJTT", "WQML", "YCKE", "ACSS", "ADCF", "ADHI", "ADMA", "ANCS", "ASCG", "CERI", "CHNI", "FLTD", "GSXD", "HAES", "ONJA", "PTSG", "TABS", "VATS", "XDDT", "ZXJO")

# Specify taxonomy table to use for filtering
cheilanthoid_taxonomy <- tibble(
  code = codes,
  group = rep("in", length(codes)),
  family = rep("pteridaceae", length(codes))
)

# Set fraction for randomly downsizing transcriptomes (e.g., 0.10 = 10%)
keep_frac = 1

# Values to use for mcl I value and Y&S hit-frac-cutoff
my_hit_frac <- 0.4
my_i_value <- 2

# Vector of reference genomes to use for making masked blast db
genomes <- c("arabidopsis", "azolla", "salvinia")

# Set number of subfolders for splitting fasta clusters
n_cluster_subfolders <- 1000
cluster_subfolders <- glue::glue("03_clusters_split/{1:n_cluster_subfolders}")

### Load plans
source("R/example_data.R")
source("R/main_plan.R")

### Download and pre-process example data
# make(example_data_plan)

### Run analyses
#
# `jobs` is set to 1, but may be increased up to the number of
# cores (CPUs) available.
#
# Check the number of cores available with this command:
# future::availableCores()
future::plan(future::multiprocess)
make(main_plan, parallelism = "future", jobs = 32)
