#!/usr/bin/env Rscript
# This script plots methylation boxplots in a (hopefully) standardised format for MEMBRA
# It should run non-interactively as part of the MEMBRA pipeline

# Get variables passed from Shell
{
  args = commandArgs (trailingOnly = TRUE)
  Yaml_filepath <- args [1]
  
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
  paste ("Scripts_dir is", Scripts_dir)
  # setwd (Scripts_dir)
}

# Load libraries
{
  ### Standard libraries
  {
    CRAN_mirror <- "https://www.stats.bris.ac.uk/R/"
    library ("tools") # file_path_sans_ext function
    library ("ggplot2") # plot things
    library ("dplyr") # various tools for manipulating R objects
    library ("tidyverse") # various tools for manipulating R objects
    library ("readr") # read tsv files
    library ("qs") # Quick saving and loading of R objects
    library ("GenomeInfoDb") # seqlengths helper function for GRanges
    library ("multcompView") # Calculate significance letters from Tukeys HSD
    library ("viridis") # Colour-blind friendly palette
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

# Load plotting data
{
  ### Create directory for storing CX report in GRanges
  {
    if (!file.exists (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports"))) {
      dir.create (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports"), recursive = TRUE)
    }
  }
  
  ### Create CX report in GRanges format and save
  {
    for (individual in All_individuals) {
      if (!file.exists (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))) {
        print (paste ("Creating CX report R object for", individual))
        CX_report <- DMRcaller::readBismark (paste0 (file.path (Working_dir, Specondition, "meth_call/CX_reports", paste0 (individual, "_trimmed_1_bismark_bt2_pe.deduplicated.CX_report.txt"))))
        qsave (CX_report, file = file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
      }
    }
    gc ()
  }
  
  ### Load and compute methylation proportions
  {
    if (file.exists (file.path (Working_dir, Specondition, "analysis/objects", paste0 ("Methylation_proportions.qs")))) {
      Methylation_proportions <- qread (file.path (Working_dir, Specondition, "analysis/objects", paste0 ("Methylation_proportions.qs")))
    } else {
      #### Create empty dataframe
      Methylation_proportions <- data.frame (
        matrix (NA, 
                nrow = length (All_individuals), 
                ncol = length (Contexts)),
        row.names = All_individuals
      )
      colnames(Methylation_proportions) <- Contexts
      
      for (individual in All_individuals) {
        print (paste ("loading CX file for", individual))
        CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
        for (context in Contexts) {
          print (paste ("computing methylation proportions in", context))
          CX_report_context <- CX_report [CX_report$context == context, ]
          CX_report_context <- CX_report_context [CX_report_context$readsN >= Bigwig_min_cov, ]
          methylation_proportion <- sum (CX_report_context$readsM) / sum (CX_report_context$readsN)
          Methylation_proportions [individual, context] <- methylation_proportion
        }
      }
      qsave (Methylation_proportions, file.path (Working_dir, Specondition, "analysis/objects", paste0 ("Methylation_proportions.qs")))
    }
    
  }
  
  ### Convert into long format suitable for ggplot
  {
    Methylation_proportions_backup <- Methylation_proportions
    Methylation_proportions <- Methylation_proportions_backup
    Methylation_proportions <- rownames_to_column (Methylation_proportions, var = "individual")
    Methylation_proportions <- pivot_longer (Methylation_proportions,
                                             cols = -individual,
                                             names_to = "context",
                                             values_to = "methylation_proportion"
    )
    All_analysis_groups_flat <- tibble::enframe (All_analysis_groups, name = "group", value = "individual") %>%
      unnest (individual)
    Methylation_proportions <- left_join (Methylation_proportions, All_analysis_groups_flat, by = "individual")
  }
}

# Make bigwig files (for IGV)
{
  function_make_bigwig <- function (CX_report, Out_file, Context, Bigwig_min_cov) {
    CX_report <- CX_report [which (CX_report$context == Context)]
    CX_report <- CX_report [which (CX_report$readsN >= Bigwig_min_cov)]
    score_CX_report <- mcols (CX_report)$readsM / mcols (CX_report)$readsN
    mcols (CX_report) <- data.frame (score = score_CX_report)
    seqlengths (CX_report) <- Chromosome_lengths$Length
    rtracklayer::export.bw (CX_report, Out_file)
  }
  
  dir.create (file.path (Working_dir, Specondition, "analysis/objects/bedgraph_files"), recursive = TRUE, showWarnings = FALSE)
  for (individual in All_individuals) {
    print (paste ("Generating BigWig files for", individual))
    for (context in Contexts) {
      Out_file <- file.path (Working_dir, Specondition, "analysis/objects/bedgraph_files", paste0 (individual, "_", context, "_min_cov_", Bigwig_min_cov, ".bw"))
      if (!file.exists (Out_file)) {
        CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
        function_make_bigwig (CX_report, Out_file, context, Bigwig_min_cov)
      }
    }
  }
}

# Create boxplot of methylation proportions
{
  dir.create (file.path (Working_dir, Specondition, "analysis/plots/Methylation_distribution"), recursive = TRUE, showWarnings = FALSE)
  for (context in Contexts) {
    ### Subset the appropriate context
    Plot_data <- Methylation_proportions [Methylation_proportions$context == context, ]
    Plot_data <- as.data.frame (Plot_data)
    if (Do_plot_order) {
      Plot_data$group <- factor(
        Plot_data$group, 
        levels = Plot_order)
    }
    
    ### Compute ANOVA and Tukey's HSD for plotting
    anova_result <- aov (methylation_proportion ~ group, data = Plot_data)
    tukey_result <- TukeyHSD (anova_result, conf.level = 0.95)
    letters <- multcompLetters4 (anova_result, tukey_result)
    letters_df <- as.data.frame (letters$group$Letters)
    letters_df <- rownames_to_column (letters_df)
    colnames (letters_df) <- c ("name", "label")
    
    ### Make the plot
    Outplot <- ggplot (Plot_data, aes (x = group, y = methylation_proportion, fill = group)) +
      geom_boxplot (outliers = FALSE) + 
      geom_point () + 
      scale_fill_viridis (discrete = TRUE, option = "D") +
      geom_text (data = letters_df, aes (x = name, y = max (Plot_data$methylation_proportion) * 1.01, label = label), inherit.aes = FALSE) +
      labs (title = paste ("Methylation Proportions in the", context, "Context for", Species, "under", Condition, "Condition"), 
            x = "Context", 
            y = "Methylation Proportion") +
      theme_bw () +
      theme (axis.text.x = element_blank (),
             axis.ticks.x = element_blank ())
    ggsave (Outplot, height = 10, width = 15 + length (unique (Plot_data$group)) * 2, units = "cm", 
            filename = file.path (Working_dir, Specondition, "analysis/plots/Methylation_distribution", paste0 ("Methylation_distribution_", context, "_boxplot.png")))
  }
}

# Plot methylation "histogram" (adapted from Marco's methods)
{
  dir.create (file.path (Working_dir, Specondition, "analysis/plots/Methylation_histogram"), recursive = TRUE, showWarnings = FALSE)
  CG_freq <- data.frame(freq = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  CHG_freq <- data.frame(freq = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  CHH_freq <- data.frame(freq = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  
  # Create table of methyaltion "histogram"
  if (!file.exists (file.path (Working_dir, Specondition, "analysis/objects/Methylation_histogram.qs"))) {
    print ("Calculating histogram of methylation")
    for (individual in All_individuals) {
      CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
      CX_report_covered <- CX_report [which (CX_report$readsN >= Methylation_min_cov)]
      
      print (paste ("Calculation of frequencies per context for", individual))
      freq_summary <- data.frame(freq = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
      for (context in Contexts) {
        freq <- c()
        percontext <- CX_report_covered[which(CX_report_covered$context==context)]
        score_percontext <- mcols(percontext)$readsM / mcols(percontext)$readsN
        
        freq[1] <- length(which(score_percontext <= 0.1)) / length(score_percontext)
        freq[2] <- length(which(score_percontext > 0.1 & score_percontext <= 0.2)) / length(score_percontext)
        freq[3] <- length(which(score_percontext > 0.2 & score_percontext <= 0.3)) / length(score_percontext)
        freq[4] <- length(which(score_percontext > 0.3 & score_percontext <= 0.4)) / length(score_percontext)
        freq[5] <- length(which(score_percontext > 0.4 & score_percontext <= 0.5)) / length(score_percontext)
        freq[6] <- length(which(score_percontext > 0.5 & score_percontext <= 0.6)) / length(score_percontext)
        freq[7] <- length(which(score_percontext > 0.6 & score_percontext <= 0.7)) / length(score_percontext)
        freq[8] <- length(which(score_percontext > 0.7 & score_percontext <= 0.8)) / length(score_percontext)
        freq[9] <- length(which(score_percontext > 0.8 & score_percontext <= 0.9)) / length(score_percontext)
        freq[10] <- length(which(score_percontext > 0.9)) / length(score_percontext)
        
        freq_summary[,context] <- freq
      }
      CG_freq[,individual] <- freq_summary$CG
      CHG_freq[,individual] <- freq_summary$CHG
      CHH_freq[,individual] <- freq_summary$CHH
      list_frequency <- list(CG_freq,CHG_freq,CHH_freq)
      names(list_frequency) <- Contexts
    }
    qsave (list_frequency, file = file.path (Working_dir, Specondition, "analysis/objects/Methylation_histogram.qs"))
  }
  else {
    list_frequency <- qread (file.path (Working_dir, Specondition, "analysis/objects/Methylation_histogram.qs"))
  }
  
  
  #plot frequencies
  for (context in Contexts) {
    ## Create colour-blind friendly colours for each group
    {
      group_colors <- viridis(length(names (All_analysis_groups)))
      # Assign colors to individuals based on their groups
      colours_vector <- unlist(lapply(seq_along(All_analysis_groups), function(i) {
        group <- names(All_analysis_groups)[i]  # Get group name
        individuals <- All_analysis_groups[[i]]  # Get individuals in the group
        setNames(rep(group_colors[i], length(individuals)), individuals)  # Assign group color
      }))
    }
    
    pdf (file.path (Working_dir, Specondition, "analysis/plots/Methylation_histogram", paste0 ("Methylation_freq_by_context_", context, ".pdf")), width=6, height=6, onefile=F)
    plot(c(0,100),c(0,0.8),type="n", 
         xlab="% of methylation", 
         ylab=paste0("m",context, "frequency"), 
         main = paste (Specondition, "Methylation proportions in", context),
         las=1)
    for (individual in All_individuals){
      lines(list_frequency[[context]]$freq*100, list_frequency[[context]][,individual], type="l",lwd=2,col=colours_vector[individual])
    }
    legend("topright",  # Position (e.g., "topleft", "topright", etc.)
           legend = names (All_analysis_groups),  # Names of the individuals
           # There is a hidden assumption of colour groups here - may need to change
           col = unique (colours_vector),  # Corresponding colors 
           lwd = 2,  # Line width in the legend
           cex = 0.8,  # Adjust text size
           title = "Groups")  # Optional title for the legend
    dev.off ()
  }
}

# Plot methylation per chromosome (adapted from Marco's methods)
{
  ## Load/create the data object for plotting
  {
    #Generate tiles
    only_chrs <- DMR_calling_chromosomes
    chr.size <- end (only_chrs)
    names (chr.size) <- DMR_calling_chromosomes@seqnames
    
    
    size <- 1000000
    overlap <- 0
    
    
    tiles <- data.frame ()
    for (i in c(1:length(chr.size))) { 
      start <- seq(1,chr.size[i] - (size-1) , by=size - overlap)
      end <- seq(size,chr.size[i], by=size - overlap)
      d <- data.frame(chr=rep(names(chr.size[i]), length(start)),start=start,end=end)
      tiles <- rbind(tiles,d)
    }
    rm(size,overlap,i,start,end,d)
    
    tiles <- format(tiles, scientific = FALSE)
    tiles$chr <- as.factor(tiles$chr)
    tiles.gr <- GenomicRanges::makeGRangesFromDataFrame(tiles)
    
    if (!file.exists (file.path (Working_dir, Specondition, "analysis/objects/tiles.qs"))) {
      #calculate methylation per tile per condition,
      for (individual in All_individuals) {
        CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports", paste0 ("CX_report_", individual, ".qs")))
        for (context in Contexts) {
          tiles.gr <- DMRcaller::analyseReadsInsideRegionsForCondition (tiles.gr, CX_report, context = context, label = individual)
        }
        tiles <- as.data.frame (tiles.gr)
      }
      qsave (tiles, file = file.path (Working_dir, Specondition, "analysis/objects/tiles.qs"))
    }
    else {
      tiles <- qread (file.path (Working_dir, Specondition, "analysis/objects/tiles.qs"))
    }
  }
  
  ## Create colour-blind friendly colours for each group
  {
    group_colors <- viridis(length(names (All_analysis_groups)))
    # Assign colors to individuals based on their groups
    colours_vector <- unlist(lapply(seq_along(All_analysis_groups), function(i) {
      group <- names(All_analysis_groups)[i]  # Get group name
      individuals <- All_analysis_groups[[i]]  # Get individuals in the group
      setNames(rep(group_colors[i], length(individuals)), individuals)  # Assign group color
    }))
  }
  
  ## Make the plot
  {
    dir.create (file.path (Working_dir, Specondition, "analysis/plots/Methylation_by_chromosome"), recursive = TRUE, showWarnings = FALSE)
    
    for (context in Contexts) {
      for (chromosome in as.character (seqnames (DMR_calling_chromosomes))) {# this bit needs fixing
        plot_data <- dplyr::filter (tiles, seqnames == chromosome)
        
        pdf (file.path (Working_dir, Specondition, paste0 ("analysis/plots/Methylation_by_chromosome/Chrom_plots_", chromosome, "_m", context, ".pdf")),
             width=(1.5*max(plot_data$end)/1000000), height=6, onefile=F)
        plot(c(0,max(plot_data$end)/1000000),c(0,1),type="n", 
             xlab="position", 
             ylab=paste0("m",context), 
             main=paste (Specondition, chromosome, "Methylation in the", context, "Context"), 
             las=1)
        
        for (individual in All_individuals){
          lines(plot_data$end/1000000, plot_data[,paste0("proportion",individual,context)], type="l",lwd=1,col=colours_vector[individual])
        }
        legend("topright",  # Position (e.g., "topleft", "topright", etc.)
               legend = names (All_analysis_groups),  # Names of the individuals
               # There is a hidden assumption of colour groups here - may need to change
               col = unique (colours_vector),  # Corresponding colors 
               lwd = 2,  # Line width in the legend
               cex = 0.8,  # Adjust text size
               title = "Groups")  # Optional title for the legend
        dev.off()
      }
    }
  }
}

# Plot PCA, scree plot and clustering of methylation (WIP)
{
  # Need to look into extending the PCA plots into multiple groups
  # Probably can be done, but harder than just extending Marco's codes...
  # {
  #   conditions_inverted <- c (Test_individuals, Control_individuals)
  #   treatments <- c (rep (1, length (Test_individuals)), rep (0, length (Control_individuals)))
  #   assembly1 <- "Hazel"
  #   file.list <- as.list (paste0 (file.path (Working_dir, Specondition, "meth_call/CX_reports", paste0 (conditions_inverted, "_trimmed_1_bismark_bt2_pe.deduplicated.CX_report.txt"))))
  #   
  #   myobj_cg <- methylKit::methRead(file.list,
  #                                   sample.id=as.list(conditions_inverted),
  #                                   assembly=assembly1,
  #                                   treatment=treatments,
  #                                   pipeline="bismarkCytosineReport",
  #                                   context="CpG")
  #   myobj_chg <- methylKit::methRead(file.list,
  #                                    sample.id=as.list(conditions_inverted),
  #                                    assembly=assembly1,
  #                                    treatment=treatments,
  #                                    pipeline="bismarkCytosineReport",
  #                                    context="CHG")
  #   myobj_chh <- methylKit::methRead(file.list,
  #                                    sample.id=as.list(conditions_inverted),
  #                                    assembly=assembly1,
  #                                    treatment=treatments,
  #                                    pipeline="bismarkCytosineReport",
  #                                    context="CHH")
  #   meth_cg <- methylKit::unite(myobj_cg, destrand=FALSE)
  #   meth_chg <- methylKit::unite(myobj_chg, destrand=FALSE)
  #   meth_chh <- methylKit::unite(myobj_chh, destrand=FALSE)
  #   
  #   pdf (file.path (Working_dir, Specondition, "analysis/plots/Methylkit_clusters.pdf"), width=18, height=6, onefile=F)
  #   par(mfrow=c(1,3))
  #   methylKit::clusterSamples(meth_cg, dist="correlation", method="ward", plot=TRUE)
  #   methylKit::clusterSamples(meth_chg, dist="correlation", method="ward", plot=TRUE)
  #   methylKit::clusterSamples(meth_chh, dist="correlation", method="ward", plot=TRUE)
  #   dev.off()
  #   
  #   pdf (file.path (Working_dir, Specondition, "analysis/plots/Methylkit_PCAbar.pdf"), width=18, height=6, onefile=F)
  #   par(mfrow=c(1,3))
  #   methylKit::PCASamples(meth_cg, screeplot=TRUE)
  #   methylKit::PCASamples(meth_chg, screeplot=TRUE)
  #   methylKit::PCASamples(meth_chh, screeplot=TRUE)
  #   dev.off()
  #   
  #   pdf (file.path (Working_dir, Specondition, "analysis/plots/Methylkit_PCA_plot.pdf"), width=18, height=6, onefile=F)
  #   par(mfrow=c(1,3))
  #   pc <- methylKit::PCASamples(meth_cg,obj.return = TRUE, adj.lim=c(1,1))
  #   pc <- methylKit::PCASamples(meth_chg,obj.return = TRUE, adj.lim=c(1,1))
  #   pc <- methylKit::PCASamples(meth_chh,obj.return = TRUE, adj.lim=c(1,1))
  #   dev.off()
  # }
}



