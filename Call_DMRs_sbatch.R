#!/usr/bin/env Rscript
# The aim of this script is to non-interactively call DMRs through Slurm

# Get variables passed from Shell
{
  args = commandArgs (trailingOnly = TRUE)
  Yaml_filepath <- args [1]
  Control_group <- args [2]
  Test_group <- args [3]
  
  print (paste ("The Yaml_filepath is", Yaml_filepath))
  print (paste ("The Control_group is", Control_group))
  print (paste ("The Test_group is", Test_group))
}

# Set paths to load/download R packages taking into account of Bluebear architecture
{
  libdir <- paste (Sys.getenv ('HOME'), "R/library", getRversion (), Sys.getenv ('BB_APPS_BASE'), sep = "/")
  if (!file.exists (libdir)) {
    dir.create (libdir, recursive = TRUE)
  }
  .libPaths (c (libdir, .libPaths ()))
}

# Define variables and load libraries
{
  ## Get parameters from yaml file and command line
  {
    library ("yaml") # read yaml files
    library ("rlang")
    library ("optparse") # parse arguments from command line
    
    ### Params
    {
      Params <- read_yaml (Yaml_filepath)
      list2env (Params, envir = globalenv ())
      rm ("Params")
      Specondition <- paste0 (Species, Condition)
      Scripts_dir <- file.path (dirname (Scripts_dir), "Data", Specondition, "Analysis_records")
      setwd (Scripts_dir)
    }
  }
  
  ## Load libraries
  {
    ### Standard libraries
    {
      CRAN_mirror <- "https://www.stats.bris.ac.uk/R/"
      library ("GenomicRanges") # Various GRanges support (this needs to go first, I think)
      library ("tools") # file_path_sans_ext function
      library ("ggplot2") # plot things
      library ("readr") # read tsv files
      library ("qs") # Quick saving and loading of R objects
      library ("GenomeInfoDb") # seqlengths helper function for GRanges
      library ("tidyr") # various tools to deal with r objects
      library ("multcompView") # Calculate significance letters from Tukeys HSD
      library ("tibble") # various dataframe tools
    }
    
    ### Bioconductor libraries
    {
      bioconductor_libraries <- c ("Biostrings", 
                                   "rtracklayer", 
                                   "DMRcaller", 
                                   "GenomicFeatures", 
                                   "GenomicRanges", 
                                   "DESeq2",
                                   "methylKit")
      if (!require ("BiocManager", quietly = TRUE)) {
        install.packages ("BiocManager")
      }
      BiocManager::install (version = "3.19", update = FALSE) # Can't write to centrally managed R package locations anyway
      for (lib in bioconductor_libraries) {
        BiocManager::install (lib, update = FALSE)
      }
    }
  }
  
  ## Set variables
  {
    ### Simple variables
    {
      # Genome_annotation_filename <- file_path_sans_ext (basename (Gene_annotation_url)) # needs fixing
      Reference_genome_filename <- file_path_sans_ext (basename (Genome_filepath))
      Contexts <- c ("CG", "CHG", "CHH")
      # condition_labels <- c (rep ("Control", length (Control_individuals)), rep ("Test", length (Test_individuals)))
      Proportion_methylated_hist_breaks <- 10 # Number of histogram blocks in proportions methylated plots
    }
    
    ### Get chromosome names and their lengths
    {
      Chromosome_lengths <- read.table (file.path (Working_dir, Specondition, "annotation/bismark_index", paste0 (Reference_genome_filename, ".txt")), header = FALSE)
      colnames (Chromosome_lengths) <- c ("Chromosome", "Length")
      if (length (DMR_calling_chromosomes) == 1 && DMR_calling_chromosomes == TRUE) {
        DMR_calling_chromosomes <- sort (unique (Chromosome_lengths$Chromosome))
        DMR_calling_chromosomes <- GenomicRanges::GRanges (seqnames = Rle (Chromosome_lengths$Chromosome), ranges = IRanges (1, end = Chromosome_lengths$Length, names = Chromosome_lengths$Chromosome), strand = "*")
      } else {
        DMR_calling_chromosomes <- Chromosome_lengths [Chromosome_lengths$Chromosome %in% DMR_calling_chromosomes, ]
        DMR_calling_chromosomes <- GenomicRanges::GRanges (seqnames = Rle (DMR_calling_chromosomes$Chromosome), ranges = IRanges (1, end = DMR_calling_chromosomes$Length, names = DMR_calling_chromosomes$Chromosome), strand = "*")
      }
    }
  }
}

# Load data
{
  Control_group <- sub ("^All_analysis_groups_", "", Control_group)
  Test_group <- sub ("^All_analysis_groups_", "", Test_group)
  Control_individuals <- All_analysis_groups [Control_group]
  Test_individuals <- All_analysis_groups [Test_group]
  
  All_CX_reports <- lapply (c (Control_individuals [[1]], Test_individuals [[1]]), function (individual) {
    print (paste ("Loading CX report for", individual))
    qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
    })
  names (All_CX_reports) <- c (Control_individuals [[1]], Test_individuals [[1]])
}

# Call DMRs using DMRCaller (based on Marco's method)
{
  dir.create (file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons")), showWarnings = FALSE, recursive = TRUE)
  chrs <- DMR_calling_chromosomes
  
  condition_labels <- c (rep (Control_group, length (Control_individuals [[1]])), rep (Test_group, length (Test_individuals [[1]])))
  for (context in Contexts) {
    if (exists ("Combined_report")) {rm ("Combined_report")}
    DMRsbinsReplicates <- GRanges ()
    for (individual in c (Control_individuals [[1]], Test_individuals [[1]])) {
      CX_report <- All_CX_reports [[which (names (All_CX_reports) == individual)]]
      CX_report <- CX_report [IRanges::overlapsAny (CX_report, chrs)]
      CX_report <- CX_report [CX_report$context == context]
      
      if (!exists ("Combined_report")) {
        Combined_report <- CX_report
      } else {
        print (paste ("Combining", individual, "with others in the", context))
        Combined_report <- DMRcaller::joinReplicates (Combined_report, CX_report, usecomplete = FALSE)
      }
    }
    
    ## Call DMRs
    print (paste ("Calculating DMRs in", context, "context"))
    DMRsbinsReplicates <- DMRcaller::computeDMRsReplicates (Combined_report,
                                                            condition = condition_labels,
                                                            regions = chrs,
                                                            context = context,
                                                            method = "bins",
                                                            binSize = binSize,
                                                            test = "betareg",
                                                            pValueThreshold = pValueThreshold,
                                                            minCytosinesCount = minCytosinesCount,
                                                            minProportionDifference = methylationDiff,
                                                            minGap = 0, # do not merge
                                                            minSize = minSize,
                                                            minReadsPerCytosine = minReadsPerCytosine,
                                                            cores = cores)
    print (paste ("Saving DMRs in the", context, "context to", file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", Control_group, "_vs_", Test_group, ".qs"))))
    qsave (DMRsbinsReplicates, file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", Control_group, "_vs_", Test_group, ".qs")))
    rtracklayer::export (DMRsbinsReplicates [which (DMRsbinsReplicates$regionType == "loss")], file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRsbinsReplicates_", context, "_", Control_group, "_vs_", Test_group, "_loss.bed")))
    rtracklayer::export (DMRsbinsReplicates [which (DMRsbinsReplicates$regionType == "gain")], file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRsbinsReplicates_", context, "_", Control_group, "_vs_", Test_group, "_gain.bed")))
  }
}
