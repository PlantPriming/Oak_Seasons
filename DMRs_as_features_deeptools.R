#!/usr/bin/env Rscript
rm (list = ls ())



# Get the DMRs from oak and use them as a feature in deeptools and compare differences between seasons and offspring
# Also get particular subsections of the CHH DMRs to use as a deeptools feature

# Get variables passed from Shell
{
  Yaml_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Control/Oak_eCO2/Oak_All_eCO2.yaml"
  print (Yaml_filepath)
}

# Set paths to load/download R packages taking into account of Bluebear architecture
{
  libdir <- paste (Sys.getenv ('HOME'), "R/library", getRversion (), Sys.getenv ('BB_APPS_BASE'), sep = "/")
  if (!file.exists (libdir)) {
    dir.create (libdir, recursive = TRUE)
  }
  .libPaths (c (libdir, .libPaths ()))
  options(bitmapType='cairo') # this is needed for plotting?
}

# Get parameters from yaml file
{
  library ("yaml")
  library ("rlang")
  
  Params <- read_yaml (Yaml_filepath)
  list2env (Params, envir = globalenv ())
  rm ("Params")
  Specondition <- paste0 (Species, Condition)
  Scripts_dir <- file.path (Working_dir, Specondition, "Analysis_records")
  setwd (Scripts_dir)
}

# Load libraries
{
  ### Standard libraries
  {
    CRAN_mirror <- "https://www.stats.bris.ac.uk/R/"
    library ("GenomeInfoDb") # seqlengths helper function for GRanges
    library ("tools") # file_path_sans_ext function
    library ("ggplot2") # plot things
    library ("dplyr") # various tools for manipulating R objects
    library ("tidyverse") # various tools for manipulating R objects
    library ("readr") # read tsv files
    library ("qs") # Quick saving and loading of R objects
    library ("multcompView") # Calculate significance letters from Tukeys HSD
    library ("viridis") # Colour-blind friendly palette
    library ("pheatmap") # Plot nice heatmaps
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
    BiocManager::install (version = "3.19", update = FALSE, force = TRUE) # Can't write to centrally managed R package locations anyway
    for (lib in bioconductor_libraries) {
      BiocManager::install (lib, update = FALSE)
    }
  }
}

# Set variables
{
  ### Simple variables
  {
    # Genome_annotation_filename <- file_path_sans_ext (basename (Gene_annotation_url)) # needs fixing
    Reference_genome_filename <- file_path_sans_ext (basename (Genome_filepath))
    CX_filepath <- file.path (Working_dir, Specondition, "meth_call", "CX_reports")
    Contexts <- c ("CG", "CHG", "CHH")
    # condition_labels <- c (rep ("Control", length (Control_individuals)), rep ("Test", length (Test_individuals)))
    Out_dir <- "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/"
  }
  
  ### Get chromosome names and their lengths
  {
    Chromosome_lengths <- read.table (file.path (Working_dir, Specondition, "annotation/bismark_index", paste0 (Reference_genome_filename, ".txt")), header = FALSE)
    colnames (Chromosome_lengths) <- c ("Chromosome", "Length")
    if (length (DMR_calling_chromosomes) == 1 && DMR_calling_chromosomes == TRUE) {
      DMR_calling_chromosomes <- sort (unique (Chromosome_lengths$Chromosome))
      DMR_calling_chromosomes <- GenomicRanges::GRanges (seqnames = Rle (Chromosome_lengths$Chromosome), ranges = IRanges::IRanges (1, end = Chromosome_lengths$Length, names = Chromosome_lengths$Chromosome), strand = "*")
    } else {
      DMR_calling_chromosomes <- Chromosome_lengths [Chromosome_lengths$Chromosome %in% DMR_calling_chromosomes, ]
      DMR_calling_chromosomes <- GenomicRanges::GRanges (seqnames = Rle (DMR_calling_chromosomes$Chromosome), ranges = IRanges::IRanges (1, end = DMR_calling_chromosomes$Length, names = DMR_calling_chromosomes$Chromosome), strand = "*")
    }
  }
}

# Load Heatmap data and make dataset for plotting
{
  All_individuals <- c (SpringAmbientParents, SummerAmbientParents, AutumnAmbientParents, AmbientOffspring)
  All_heatmap_data <- lapply (All_individuals, function (individual) {
    qread (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps", paste0 ("Heatmap_recomputed_ambient_par_and_ambient_off_", individual, ".qs")))
  })
  names (All_heatmap_data) <- All_individuals
  
  
  Pairwise_combinations <- combn (names (All_analysis_groups), 2, paste, collapse = "_vs_")  
  Pairwise_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
  
  ## Get all DMRs and write to file
  {
    All_DMRs <- lapply (Contexts, function (context) {
      Pairwise_comparison_DMRs <- lapply (Pairwise_combinations, function (combination) {
        if (file.exists (file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", combination, ".qs")))) {
          qread (file.path (Working_dir, Specondition, paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", combination, ".qs")))
        }
      })
      names (Pairwise_comparison_DMRs) <- Pairwise_combinations
      return (Pairwise_comparison_DMRs)
    })
    names (All_DMRs) <- Contexts
    
    ## Generate list of all DMRs found in multiple comparisons
    {
      All_conditions_DMRs <- lapply (All_DMRs, function (context) {
        ### Maybe need the unique function, but needs testing to ensure that only ranges (context) are considered for unique 
        DMRs <- unique (Reduce (append, context))
        mcols (DMRs) <- NULL
        return (DMRs)
      })
    }
  }
}

# Get Overlapping DMRs and make into bed files for deeptools
{
  Contexts <- c ("CG", "CHG", "CHH")
  for (context in Contexts) {
    ## For the intersect of Spring vs Summer and Spring vs Autumn
    {
      Intersect_DMRs <- All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents [overlapsAny (All_DMRs[[context]]$SummerAmbientParents_vs_SpringAmbientParents, All_DMRs[[context]]$AutumnAmbientParents_vs_SpringAmbientParents)]
      Out_data_filepath <- paste0 ("/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/DMR_bed_files/Spring_Summer_Spring_Autumn_intersect_", context, ".bed")
      rtracklayer::export (Intersect_DMRs, Out_data_filepath, format = "bed")
    }
    
    ## For Spring vs Summer only
    {
      ### Get the DMRs
      {
        Spring_summer_only_DMRs <- GenomicRanges::setdiff (All_DMRs[[context]]$SummerAmbientParents_vs_SpringAmbientParents, All_DMRs[[context]]$AutumnAmbientParents_vs_SpringAmbientParents)
        Out_data_filepath <- paste0 ("/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/DMR_bed_files/Spring_Summer_only_", context, ".bed")
        rtracklayer::export (Spring_summer_only_DMRs, Out_data_filepath, format = "bed")
      }
    }
    
    ## For Spring vs Autumn only
    {
      Spring_autumn_only_DMRs <- GenomicRanges::setdiff (All_DMRs[[context]]$AutumnAmbientParents_vs_SpringAmbientParents, All_DMRs[[context]]$SummerAmbientParents_vs_SpringAmbientParents)
      Out_data_filepath <- paste0 ("/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/DMR_bed_files/Spring_Autumn_only_", context, ".bed")
      rtracklayer::export (Spring_autumn_only_DMRs, Out_data_filepath, format = "bed")
    }
  }
}

# Write to bed file
{
  for (context in Contexts) {
    Out_data <- All_conditions_DMRs [[context]]
    if (length (Out_data) > 10000){
      Out_data <- Out_data [sample (length (Out_data), 10000, replace = FALSE)]
    }
    Out_data_filepath <- file.path ("/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/DMR_bed_files", paste0 ("Ambient_par_offspring_10000_rand_", context, ".bed"))
    rtracklayer::export (Out_data, Out_data_filepath, format = "bed")
  }
}