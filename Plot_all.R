# Try to plot all figures in a publishable way
# 1. Repeat stacked bar charts
# 2. Heatmap of Spring, and Progeny in CHH in white/red
# 3. Ontology of CDS and Promoters overlapping with DMRs
# 4. Venn Diagrams move overlapping labels
# 5. Deeptools 


# Get parameters from yaml file, load libraries, set variables and load data
{
  ## Get parameters from yaml file
  {
    library ("yaml")
    library ("rlang")
    
    Yaml_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Control/Oak_eCO2/Oak_All_eCO2.yaml"
    Params <- read_yaml (Yaml_filepath)
    
    list2env (Params, envir = globalenv ())
    rm ("Params")
    Specondition <- paste0 (Species, Condition)
    Scripts_dir <- file.path (Scripts_dir, "Data", Specondition, "Analysis_records")
    Scripts_dir <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/Analysis_records"
    Out_plot_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/Plots"
    setwd (Scripts_dir)
  }
  
  ## Load libraries
  {
    ### Standard libraries
    {
      CRAN_mirror <- "https://www.stats.bris.ac.uk/R/"
      library ("tools") # file_path_sans_ext function
      library ("ggplot2") # plot things
      library ("readr") # read tsv files
      library ("qs") # Quick saving and loading of R objects
      library ("GenomeInfoDb") # seqlengths helper function for GRanges
      library ("GenomicRanges") # Various GRanges support
      library ("dplyr") # manipulation of dataframes
      library ("viridis") # Colour-blind friendly palette
      library ("VennDiagram") # Venn Diagrams
      library ("grid") # For putting venn diagram into pdf
      library ("gridExtra")
      library ("Rgraphviz") # needed for Jack's functions
      library ("topGO") # GO
      library ("tidyr") # Various useful functions
      library ("tibble")
      library ("pheatmap")
      library ("stringr")
      library("tidyverse")
      library("scales") # make pretty breaks in scales
      options(bitmapType='cairo') # Needed for plotting apparently
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
      promoter_width <- 1000 # Assume this many bp upstream of CDS represents the promoter region
      GFF_filepath_TEs <- "/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/Qrob_PM1N_refTEs.gff"
      GFF_filepath_genes <-"/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/Qrob_PM1N_genes_20161004.gff"
      
      # Get this number of random regions to be the overlap controls
      sample_number <- 10000
      # with this region width. Should set equal to DMR size for control comparison
      region_width <- 400
      
      # Output GO analysis here
      GO_out_folder <- "/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/"
      GO_database <- "/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/trinotate_GO_terms.tsv"
      gff_file <- "/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/Qrob_PM1N_genes_20161004.gff"
      GFF_filepath_genes <-"/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/Qrob_PM1N_genes_20161004.gff"
      
      # Mappability and TE analysis
      Oak_mapability_filepath <- "/rds/projects/c/catonim-minicircles/mappability/results/Qrob_PM1N.genmap.txt"
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

## Load DMRs
{
  Out_dir <- "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots"
  Pairwise_comparison_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/"
  Groups <- c ("SpringElevatedParents", "SummerElevatedParents", "AutumnElevatedParents", "ElevatedOffspring", "SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")
  
  dir.create (Out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Load DMRs
  Pairwise_combinations <- combn (Groups, 2, paste, collapse = "_vs_")
  Reversed_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
}

## Define functions
{
  ### Make random control
  {
    #### Function to generate random regions for a given set of chromosomes and their lengths with a specific width
    #### Used as a control when overlapped with gene features to determine if DMR are overlapping disproportionately with features
    function_generate_random_regions <- function (sample_number, Chromosome_lengths, region_width = 400) {
      random_chromosomes <- sample (x = Chromosome_lengths$Chromosome, 
                                    size = sample_number, 
                                    replace = TRUE, 
                                    prob = ifelse ((Chromosome_lengths$Length - region_width) < 0, 0, (Chromosome_lengths$Length - region_width)))
      random_positions <- sapply (random_chromosomes, 
                                  function (chromosome) 
                                  {sample ((Chromosome_lengths [Chromosome_lengths$Chromosome == chromosome, "Length"] - region_width), 1)})
      GRanges (seqnames = names (random_positions),
               ranges = IRanges (start = random_positions, end = random_positions + region_width))
    }
    
    ## Function to generate random regions
    function_generate_random_regions_2 <- function (sample_number, Chromosome_lengths, region_width = 400, DMR_calling_chromosomes) {
      Chromosome_lengths <- Chromosome_lengths [Chromosome_lengths$Chromosome %in% DMR_calling_chromosomes, ]
      random_chromosomes <- sample (x = Chromosome_lengths$Chromosome, 
                                    size = sample_number, 
                                    replace = TRUE, 
                                    prob = ifelse ((Chromosome_lengths$Length - region_width) < 0, 0, (Chromosome_lengths$Length - region_width)))
      random_positions <- sapply (random_chromosomes, 
                                  function (chromosome) 
                                  {sample ((Chromosome_lengths [Chromosome_lengths$Chromosome == chromosome, "Length"] - region_width), 1)})
      GenomicRanges::GRanges (seqnames = names (random_positions),
                              ranges = IRanges (start = random_positions, end = random_positions + region_width - 1))
    }
  }
  
  ### Get ranges
  {
    #### Get the ranges from a GenomicRanges object in a format suitable for VennDiagram plotting
    get_ranges <- function (Granges) {
      paste (seqnames (Granges), start (Granges), end (Granges), sep = ":")
    }
  }
  
  ### Find groups
  {
    #### For an individual, check to see if the individual is in a analysis group (e.g. season)
    #### Allows labelling of heatmap plot by season automatically for individuals
    find_group <- function (individual, All_analysis_groups = All_analysis_groups) {
      # Iterate through each group and check if the individual is in that group
      group <- NULL
      for (name in names (All_analysis_groups)) {
        if (individual %in% All_analysis_groups [[name]]) {
          group <- name
          break
        }
      }
      return(group)
    }
  }
  
  ### Gene Ontology
  {
    GO_enrichment <- function(genes_of_interest, 
                              gff_file, 
                              GO_database, 
                              out_folder,
                              fisher_threshold = 0.05,
                              min_sig_genes = 2,
                              useInfo = "all"){
      
      # Load and prep data #########
      
      dir.create (out_folder, recursive = TRUE, showWarnings = FALSE)
      
      # Load oak GO database
      GO_database <- readMappings(file = GO_database, 
                                  sep = "\t", IDsep = ",")
      
      # Load gff file
      gff <- rtracklayer::import(gff_file)
      
      #count number of unique genes in gff
      number_of_genes <- 
        gff$ID %>%
        str_extract(., "[^_]*_[^_]*") %>%
        unique()%>%
        length()
      
      #extract unique gene names from gff
      all_gene_names <-
        gff$ID %>%
        str_extract(., "[^_]*_[^_]*") %>%
        unique()
      
      # load DMRs
      DMRs <- load(file = genes_of_interest)
      DMRs <- get(DMRs)
      
      #extract file name for saving data
      out_file_ID = basename(genes_of_interest) 
      out_file_ID = tools::file_path_sans_ext(out_file_ID)
      
      # extract gene names from DMRs
      test_gene_names <- stringr::str_extract(DMRs$ID, "[^_]*_[^_]*")
      #POSSIBLY ADD A WAY TO EXTRACT JUST CDS SEQUENCES HERE
      test_gene_names <- unique(test_gene_names) #select only unique entries
      
      #set up vector to store genes of interest for GO analysis
      test_genes <- vector(length = number_of_genes)
      
      #set genes of interest to 1 and others to 0
      test_genes <- 
        factor(as.integer(all_gene_names %in% test_gene_names))
      
      #add genes names to vector
      names(test_genes) <- all_gene_names
      
      # BP (biological process) ##############
      
      #create topGOdata object to be used for enrichment analysis
      BP_topGOdata <- new("topGOdata", 
                          ontology = "BP", 
                          allGenes = test_genes, 
                          annot = annFUN.gene2GO, 
                          gene2GO = GO_database)
      
      #perform enrichment analysis
      BP_fisher_test <- runTest(BP_topGOdata, 
                                algorithm = "weight01", 
                                statistic = "fisher") 
      BP_fisher_test
      
      #visualise enrichment analysis results (ordered by fisher p-value)
      BP_GO_results <- GenTable(BP_topGOdata, 
                                fisher = BP_fisher_test, 
                                orderBy = "fisher", 
                                topNodes = length(score(BP_fisher_test)), 
                                numChar=1000)
      
      #plot and save GO hierarchy graph
      
      #select significant GOs for hierarchy graph
      sig_GOs <- BP_GO_results %>%
        filter(Significant >= min_sig_genes) %>%
        filter(fisher <= fisher_threshold) %>%
        pull(GO.ID)
      
      #only plot hierarchy graph if at least 1 sig_GO found
      if(length(sig_GOs) == 0){
        cat("No BP GO terms with at least", min_sig_genes, "significant genes at a fisher threshold of", fisher_threshold, "found.\n")
      } else {
        
        #filter scores to only significant GOs for hierarchy graph
        filtered_scores <- score(BP_fisher_test)[sig_GOs]
        
        pdf(file = paste0(out_folder, "BP_GO_hierarchy_", out_file_ID, ".pdf"))
        
        showSigOfNodes(BP_topGOdata, 
                       filtered_scores, 
                       firstSigNodes = length(sig_GOs),
                       useInfo = useInfo ) 
        
        dev.off()
        
      }
      
      # MF (molecular function) ##############
      
      #create topGOdata object to be used for enrichment analysis
      MF_topGOdata <- new("topGOdata", 
                          ontology = "MF", 
                          allGenes = test_genes, 
                          annot = annFUN.gene2GO, 
                          gene2GO = GO_database)
      
      #perform enrichment analysis
      MF_fisher_test <- runTest(MF_topGOdata, 
                                algorithm = "weight01", 
                                statistic = "fisher") 
      MF_fisher_test
      
      #visualise enrichment analysis results (ordered by fisher p-value)
      MF_GO_results <- GenTable(MF_topGOdata, 
                                fisher = MF_fisher_test, 
                                orderBy = "fisher", 
                                topNodes = length(score(MF_fisher_test)), 
                                numChar=1000)
      
      #plot and save GO hierarchy graph
      
      #select significant GOs for hierarchy graph
      sig_GOs <- MF_GO_results %>%
        filter(Significant >= min_sig_genes) %>%
        filter(fisher <= fisher_threshold) %>%
        pull(GO.ID)
      
      #only plot hierarchy graph if at least 1 sig_GO found
      if(length(sig_GOs) == 0){
        cat("No MF GO terms with at least", min_sig_genes, "significant genes at a fisher threshold of", fisher_threshold, "found.\n")
      } else {
        
        #filter scores to only significant GOs for hierarchy graph
        filtered_scores <- score(MF_fisher_test)[sig_GOs]
        
        pdf(file = paste0(out_folder, "MF_GO_hierarchy_", out_file_ID, ".pdf"))
        
        showSigOfNodes(MF_topGOdata, 
                       filtered_scores, 
                       firstSigNodes = length(sig_GOs),
                       useInfo = useInfo ) 
        
        dev.off()
        
      }
      
      # CC (cellular compartment) ##############
      
      #create topGOdata object to be used for enrichment analysis
      CC_topGOdata <- new("topGOdata", 
                          ontology = "CC", 
                          allGenes = test_genes, 
                          annot = annFUN.gene2GO, 
                          gene2GO = GO_database)
      
      #perform enrichment analysis
      CC_fisher_test <- runTest(CC_topGOdata, 
                                algorithm = "weight01", 
                                statistic = "fisher") 
      CC_fisher_test
      
      #visualise enrichment analysis results (ordered by fisher p-value)
      CC_GO_results <- GenTable(CC_topGOdata, 
                                fisher = CC_fisher_test, 
                                orderBy = "fisher", 
                                topNodes = length(score(CC_fisher_test)), 
                                numChar=1000)
      
      #plot and save GO hierarchy graph
      
      #select significant GOs for hierarchy graph
      sig_GOs <- CC_GO_results %>%
        filter(Significant >= min_sig_genes) %>%
        filter(fisher <= fisher_threshold) %>%
        pull(GO.ID)
      
      #only plot hierarchy graph if at least 1 sig_GO found
      if(length(sig_GOs) == 0){
        cat("No CC GO terms with at least", min_sig_genes, "significant genes at a fisher threshold of", fisher_threshold, "found.\n")
      } else {
        
        #filter scores to only significant GOs for hierarchy graph
        filtered_scores <- score(CC_fisher_test)[sig_GOs]
        
        pdf(file = paste0(out_folder, "CC_GO_hierarchy_", out_file_ID, ".pdf"))
        
        showSigOfNodes(CC_topGOdata, 
                       filtered_scores, 
                       firstSigNodes = length(sig_GOs),
                       useInfo = useInfo ) 
        
        dev.off()
        
      }
      
      # Save final enrichment results #########
      
      # save results to list
      enrichment_results <- list(BP_GO_results, 
                                 MF_GO_results, 
                                 CC_GO_results,
                                 BP_topGOdata,
                                 MF_topGOdata,
                                 CC_topGOdata,
                                 BP_fisher_test,
                                 MF_fisher_test,
                                 CC_fisher_test)
      
      # add names to list to make results clearer
      names(enrichment_results) <- c("BP_GO_results", 
                                     "MF_GO_results", 
                                     "CC_GO_results",
                                     "BP_topGOdata",
                                     "MF_topGOdata",
                                     "CC_topGOdata",
                                     "BP_fisher_test",
                                     "MF_fisher_test",
                                     "CC_fisher_test")
      
      #save results as rdata file
      save(enrichment_results, file = paste0(out_folder, "Enrichment_results_", out_file_ID, ".RData"))
      
      #return results
      return(enrichment_results)
      
    }
  }
  
  ### GO Ontology Plotting
  {
    GO_plots <- function(enrichment_RData, 
                         min_sig_genes = 2, 
                         fisher_threshold = 0.05,
                         plot_width = 30,
                         plot_height = NULL,
                         max_y_lim = 1000, 
                         plot_title,
                         show_ontology = TRUE,
                         out_folder){
      
      # load enrichment
      enrichment <- load(file = enrichment_RData)
      enrichment <- get(enrichment)
      
      #extract file name for saving data
      out_file_ID = basename(enrichment_RData) 
      out_file_ID = tools::file_path_sans_ext(out_file_ID)
      
      #convert data frames for plotting to list
      BP_GO_results <- enrichment[["BP_GO_results"]]
      MF_GO_results <- enrichment[["MF_GO_results"]]
      CC_GO_results <- enrichment[["CC_GO_results"]]
      GO_results_list <- list(BP = BP_GO_results,
                              MF = MF_GO_results,
                              CC = CC_GO_results)
      
      # Initialize an empty list to store the plots
      GO_graphs <- list()
      
      # Filter data and determine the maximum significant value across all categories
      #this will be used to keep a consistent y-axis scale across all graphs
      max_significant <- max(
        sapply(GO_results_list, function(df) {
          max(df %>% filter(Significant > 2, fisher <= 0.05) %>% pull(Significant), na.rm = TRUE)
        }),
        na.rm = TRUE
      ) + 10 #add 10 to give a bit of leeway at top of graph
      
      # Loop over each data frame
      for (category in names(GO_results_list)) {
        
        #convert fisher value to numeric
        GO_results_list[[category]]$fisher <- as.numeric(GO_results_list[[category]]$fisher)
        
        #filter data
        filtered_data <- GO_results_list[[category]] %>%
          filter(Significant >= min_sig_genes) %>%
          filter(fisher <= fisher_threshold) 
        
        #set up max y-axis value
        max_filtered_significant <- ifelse(nrow(filtered_data) > 0, 
                                           max(filtered_data$Significant, na.rm = TRUE), 10)
        
        
        custom_theme <- theme (
          axis.text.y = element_text (size = 18),
          # legend.position = "none"
        )
        if (show_ontology) {
          custom_theme <- custom_theme + theme (axis.ticks.y = element_blank ())
        } else {
          custom_theme <- custom_theme + theme (axis.text.y = element_blank (), , axis.title.y = element_blank())
        }
        
        #plot data
        GO_graphs[[category]] <- filtered_data %>%
          ggplot(aes(x = fct_reorder(Term, fisher), y = Significant, fill = fisher)) +
          geom_col() +
          ggtitle (plot_title) +
          xlab(paste(category, "GO Term")) +
          ylab("Number of Significant Genes") +
          #scale_y_continuous(limits = c(0, max_significant)) + #make it so all plots have same max-value on y-axis
          scale_y_log10(limits = c(1, max_y_lim)) +
          theme_bw () + 
          custom_theme +
          guides(fill = guide_colorbar(title = "Fisher P-value")) +
          coord_flip () +
          scale_fill_viridis(limits = c(0, fisher_threshold))
        
        #save graphs
        if (is.null (plot_height)) {
          plot_height <- 1 * nrow (filtered_data) + 2
        }
        
        ggsave(filename = paste0(out_folder, category, "_GO_graph_", out_file_ID, show_ontology, ".png"), 
               plot = GO_graphs[[category]], 
               width = plot_width, 
               height = plot_height, units = "cm")
        ggsave(filename = paste0(out_folder, category, "_GO_graph_", out_file_ID, show_ontology, ".pdf"), 
               plot = GO_graphs[[category]], 
               width = plot_width, 
               height = plot_height, units = "cm")
      }
      
      return(GO_graphs)
    }
  }
  
  ### function_get_TE_from_granges
  {
    # Use an input matrix with names in format: chromosome_start_end and granges object to get only the subbset of granges object that match
    function_get_TE_from_granges <- function (input_matrix, granges_object) {
      split_rn <- do.call (rbind, strsplit (row.names (input_matrix), "_"))
      colnames (split_rn) <- c("prefix", "chr", "start", "end")
      rn_info <- data.frame(
        chr = paste0(split_rn[,1], "_", split_rn[,2]),
        start = as.numeric(split_rn[,3]),
        end = as.numeric(split_rn[,4]),
        row.names = paste(split_rn[,1], split_rn[,2], split_rn[,3], split_rn[,4], sep = "_")
      )
      
      rn_info <- GRanges (rn_info)
      
      # Then subset
      matches <- findOverlaps (rn_info, granges_object, type = "equal")
      GRanges_subset <- granges_object [subjectHits (matches)]
      return (GRanges_subset)
    }
  }
  
  ### second_to_last_char
  {
    second_to_last_char <- function(x) {
      # Length of the string
      len <- nchar(x)
      # Extract the second to last character
      substr(x, 2, len)
    }
  }
}

# Figure 1
{
  # Plot methylation per chromosome (adapted from Marco's methods)
  {
    
    # Old code
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
        All_analysis_groups <- All_analysis_groups [c("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")]
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
        dir.create ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/Methylation_by_chrom", recursive = TRUE, showWarnings = FALSE)
        
        All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
        for (context in Contexts) {
          for (chromosome in as.character (seqnames (DMR_calling_chromosomes))) {# this bit needs fixing
            plot_data <- dplyr::filter (tiles, seqnames == chromosome)
            
            pdf (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/Methylation_by_chrom", paste0 ("Chrom_plots_", chromosome, "_m", context, ".pdf")),
                 width=10, height=6, onefile=F)
            plot(c(0,max(plot_data$end)/1000000),c(0,ifelse (context == "CHH", 0.15, 1)),type="n", 
                 xlab="position/Mbp", 
                 ylab=paste0("m",context), 
                 main=paste ("Methylation Profile of", chromosome, "in the", context, "Context"), 
                 las=1)
            
            for (individual in All_individuals){
              lines(plot_data$end/1000000, plot_data[,paste0("proportion",individual,context)], type="l",lwd=1,col=colours_vector[individual])
            }
            legend("topright",  # Position (e.g., "topleft", "topright", etc.)
                   legend = c ("Spring", "Summer", "Autumn", "Offspring"),  # Names of the individuals
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
    
    
    ## Make the plot (new code in ggplot)
    {
      dir.create ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/Methylation_by_chrom", recursive = TRUE, showWarnings = FALSE)
      All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
      for (context in Contexts) {
        for (chromosome in as.vector (seqnames (DMR_calling_chromosomes))) {
          plot_data <- dplyr::filter (tiles, seqnames == chromosome)
          plot_data <- plot_data [, paste0 ("proportion", All_individuals, context)]
          
          season_info <- data.frame (
            Sample = paste0 ("proportion", All_individuals, context),
            Season = c ("Spring", "Spring", "Spring", "Spring", "Spring", "Spring", "Summer", "Summer", "Autumn", "Autumn", "Autumn", "Offspring", "Offspring", "Offspring")  # Assign seasons
          )
          
          plot_data <- mutate (plot_data, Index = row_number())
          plot_data <- as.data.frame (pivot_longer (plot_data, cols = -Index))
          plot_data <- left_join (plot_data, season_info, by = c ("name" = "Sample"))
          
          plot_order <- c ("Spring","Summer", "Autumn", "Offspring")
          plot_data$Season <- factor (plot_data$Season, levels = plot_order)
          
          Out_plot <- ggplot (data = plot_data, aes (x = Index, y = value, group = name, colour = Season)) + 
            geom_line (aes ()) +
            scale_color_viridis_d (option = "viridis") + 
            labs (title = paste (chromosome, context)) + 
            xlab ("Position/Mbp") +
            ylab ("Methylation Proportion") +
            labs (color = "Group") +
            ylim (c (0, ifelse (context == "CHH", 0.15, 1))) +
            theme_bw () + 
            theme (plot.title = element_text (hjust = 0.5),
                   legend.position = "none"
            )
          ggsave (filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/Methylation_by_chrom/", chromosome, context, "_methylation.png"),
                  height = 10, width = 7, units = "cm")
        }
      }
    }
  }
}

# Figure 2
{
  ## Bar chart of DMR counts
  {
    ### Load data of DMRs
    {
      Pairwise_combinations <- combn (names (All_analysis_groups), 2, paste, collapse = "_vs_")
      Pairwise_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
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
      
      DMRs_barplot_data <- lapply (Contexts, function (context) {
        barplot_data <- All_DMRs [[context]]
        Gain_loss_DMRs <- lapply (barplot_data, function (pairwise_comparison) {
          table (pairwise_comparison$regionType)
        })
        flattened_df <- do.call (rbind, lapply (names (Gain_loss_DMRs), function (name) {
          data.frame(
            Comparison = name,
            Context = context,
            Gain = as.numeric (ifelse (!is.na (Gain_loss_DMRs[[name]]["gain"]), Gain_loss_DMRs[[name]]["gain"], 0)),
            Loss = as.numeric (ifelse (!is.na (Gain_loss_DMRs[[name]]["loss"]), Gain_loss_DMRs[[name]]["loss"], 0))
          )
        }))
        flattened_df$Loss <- flattened_df$Loss * -1
        flattened_df <- as.data.frame (pivot_longer (flattened_df, cols = c ("Gain", "Loss")))
      })
      names (DMRs_barplot_data) <- Contexts
    }
    
    ### Make plots
    {
      ### Make custom dataset to plot
      {
        #### Get Ambient Datasets only
        {
          DMRs_barplot_data_par <- list ()
          DMRs_barplot_data_par$CG <- DMRs_barplot_data$CG [DMRs_barplot_data$CG$Comparison %in% 
                                                              c ("SummerAmbientParents_vs_SpringAmbientParents", 
                                                                 "AutumnAmbientParents_vs_SpringAmbientParents", 
                                                                 "AutumnAmbientParents_vs_SummerAmbientParents"), ]
          DMRs_barplot_data_par$CHG <- DMRs_barplot_data$CHG [DMRs_barplot_data$CHG$Comparison %in% 
                                                                c ("SummerAmbientParents_vs_SpringAmbientParents", 
                                                                   "AutumnAmbientParents_vs_SpringAmbientParents", 
                                                                   "AutumnAmbientParents_vs_SummerAmbientParents"), ]
          DMRs_barplot_data_par$CHH <- DMRs_barplot_data$CHH [DMRs_barplot_data$CHH$Comparison %in% 
                                                                c ("SummerAmbientParents_vs_SpringAmbientParents", 
                                                                   "AutumnAmbientParents_vs_SpringAmbientParents", 
                                                                   "AutumnAmbientParents_vs_SummerAmbientParents"), ]
          
          DMRs_barplot_data_off <- list ()
          DMRs_barplot_data_off$CG <- DMRs_barplot_data$CG [DMRs_barplot_data$CG$Comparison %in% 
                                                              c ("AmbientOffspring_vs_SpringAmbientParents", 
                                                                 "AmbientOffspring_vs_SummerAmbientParents", 
                                                                 "AmbientOffspring_vs_AutumnAmbientParents"), ]
          DMRs_barplot_data_off$CHG <- DMRs_barplot_data$CHG [DMRs_barplot_data$CHG$Comparison %in% 
                                                                c ("AmbientOffspring_vs_SpringAmbientParents", 
                                                                   "AmbientOffspring_vs_SummerAmbientParents", 
                                                                   "AmbientOffspring_vs_AutumnAmbientParents"), ]
          DMRs_barplot_data_off$CHH <- DMRs_barplot_data$CHH [DMRs_barplot_data$CHH$Comparison %in% 
                                                                c ("AmbientOffspring_vs_SpringAmbientParents", 
                                                                   "AmbientOffspring_vs_SummerAmbientParents", 
                                                                   "AmbientOffspring_vs_AutumnAmbientParents"), ]
        }
      }
      
      for (context in Contexts) {
        #### Reverse some Comparisons Made
        {
          Plot_dataframe <- DMRs_barplot_data_par [[context]]
          Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "SummerAmbientParents_vs_SpringAmbientParents", "Comparison"] <- "SpringAmbientParents_vs_SummerAmbientParents"
          
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SpringAmbientParents", "Comparison"] <- "SpringAmbientParents_vs_AutumnAmbientParents"
          
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "AutumnAmbientParents_vs_SummerAmbientParents", "Comparison"] <- "SummerAmbientParents_vs_AutumnAmbientParents"
          
          Plot_dataframe$Comparison <- factor (Plot_dataframe$Comparison, levels = c ("SpringAmbientParents_vs_SummerAmbientParents", 
                                                                                      "SpringAmbientParents_vs_AutumnAmbientParents", 
                                                                                      "SummerAmbientParents_vs_AutumnAmbientParents"))
        }
        
        Out_plot <- ggplot (data = Plot_dataframe, aes (x = Comparison, y = value, fill = name)) +
          geom_bar (stat = "identity", position = "stack") + 
          labs (title = NULL,
                x = NULL, 
                y = paste (context, "DMR count")) + 
          lims (y = c (-100* 1.1, max (abs (DMRs_barplot_data_off$CG$value)) * 1.1)) + 
          scale_fill_viridis (discrete = TRUE, option = "inferno", begin = 0.2, end = 0.8, direction = -1) + 
          scale_x_discrete (labels = c ("Spring vs\nSummer", "Spring vs\nAutumn", "Summer vs\nAutumn")) + 
          labs (fill = "DMR") +
          theme_bw () + 
          theme (axis.text.x = element_text (angle = 90, hjust = 0.5, vjust = 0.5),
                 axis.ticks.x = element_blank ()
          )
        ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_par", "_", context, ".png"),  height = 10, width = 7, units = "cm")
        ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_par", "_", context, ".pdf"),  height = 10, width = 7, units = "cm")
        
        
        #### Reverse some Comparisons Made
        {
          Plot_dataframe <- DMRs_barplot_data_off [[context]]
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SpringAmbientParents", "Comparison"] <- "SpringAmbientParents_vs_AmbientOffspring"
          
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_SummerAmbientParents", "Comparison"] <- "SummerAmbientParents_vs_AmbientOffspring"
          
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents" & Plot_dataframe$name == "Gain", ]$name <- "Loss_2"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents" & Plot_dataframe$name == "Loss", ]$name <- "Gain"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents" & Plot_dataframe$name == "Loss_2", ]$name <- "Loss"
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents", "value" ] <- Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents", "value" ] * -1
          Plot_dataframe [Plot_dataframe$Comparison == "AmbientOffspring_vs_AutumnAmbientParents", "Comparison"] <- "AutumnAmbientParents_vs_AmbientOffspring"
          
          Plot_dataframe$Comparison <- factor (Plot_dataframe$Comparison, levels = c ("SpringAmbientParents_vs_AmbientOffspring", 
                                                                                      "SummerAmbientParents_vs_AmbientOffspring", 
                                                                                      "AutumnAmbientParents_vs_AmbientOffspring"))
        }
        
        Out_plot <- ggplot (data = Plot_dataframe, aes (x = Comparison, y = value, fill = name)) +
          geom_bar (stat = "identity", position = "stack") + 
          labs (title = NULL,
                x = NULL, 
                y = paste (context, "DMR count")) + 
          lims (y = c (-max (abs (DMRs_barplot_data_off$CG$value)) * 1.1, max (abs (DMRs_barplot_data_off$CG$value)) * 1.1)) + 
          scale_fill_viridis (discrete = TRUE, option = "inferno", begin = 0.2, end = 0.8, direction = -1) + 
          scale_x_discrete (labels = c ("Spring vs\nOffspring", "Summer vs\nOffspring", "Autumn vs\nOffspring")) + 
          labs (fill = "DMR") +
          theme_bw () + 
          theme (axis.text.x = element_text (angle = 90, hjust = 0.5, vjust = 0.5),
                 axis.ticks.x = element_blank ()
          )
        ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_off","_", context, ".png"),  height = 10, width = 7, units = "cm")
        ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_off","_", context, ".pdf"),  height = 10, width = 7, units = "cm")
        
      }
    }
  }
  
  ## Venn Diagrams
  {
    Contexts <- c ("CG", "CHG", "CHH")
    Pairwise_comparison_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/"
    
    for (context in Contexts) {
      print (paste ("Working on", context))
      All_context_DMRs <- lapply (Reversed_combinations, function (DMR_pair) {
        if (file.exists (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
          qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
        } else if (context == "Random") {
          function_generate_random_regions (sample_number, Chromosome_lengths, region_width)
        }
      })
      names (All_context_DMRs) <- Reversed_combinations
      
      Spring_summer_DMRs <- All_context_DMRs [["SummerAmbientParents_vs_SpringAmbientParents"]]
      Spring_summer_DMRs <- get_ranges (Spring_summer_DMRs)
      
      Spring_autumn_DMRs <- All_context_DMRs [["AutumnAmbientParents_vs_SpringAmbientParents"]]
      Spring_autumn_DMRs <- get_ranges (Spring_autumn_DMRs)
      
      Summer_autumn_DMRs <- All_context_DMRs [["AutumnAmbientParents_vs_SummerAmbientParents"]]
      Summer_autumn_DMRs <- get_ranges (Summer_autumn_DMRs)
      
      Venn_data <- list(
        "Spring_summer_DMRs" = unique (Spring_summer_DMRs),
        "Spring_autumn_DMRs" = unique (Spring_autumn_DMRs),
        "Summer_autumn_DMRs" = unique (Summer_autumn_DMRs)
      )
      
      output_path <- file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all")
      venn_filename_base <- paste0 ("Venn_Seasons_", context)
      
      venn_diagram <- venn.diagram (
        x = Venn_data,
        category.names = c ("Spring\nvs\nSummer", "Spring\nvs\nAutumn", "Summer\nvs\nAutumn"),
        fill = viridis (3, option = "turbo"),
        alpha = 0.2,
        print.mode = c ("raw", "percent"),
        main = context,
        main.cex = 1.5, 
        main.fontface = "bold", 
        filename = file.path (output_path, paste0 (venn_filename_base, ".png")),
        width = 15,
        height = 15,
        units = "cm",
        margin = 0.05,
        output = TRUE,
        cex = 1,
        cat.cex = 1,
        cat.dist = c (0.11, 0.11, 0.11)
      )
      
      pdf_file <- file.path(output_path, paste0(venn_filename_base, ".pdf"))
      venn_diagram <- venn.diagram(
        x = Venn_data,
        category.names = c ("Spring\nvs\nSummer", "Spring\nvs\nAutumn", "Summer\nvs\nAutumn"),
        fill = viridis (3, option = "turbo"),
        alpha = 0.2,
        print.mode = c ("raw", "percent"),
        main = context,
        filename = NULL,
        output = FALSE,
        width = 10,
        height = 10,
        units = "cm",
        margin = 0.1,
        output = TRUE,
        cex = 1.5,
        cat.cex = 1.5,
        cat.dist = c (0.13, 0.13, 0.13)
      )
      pdf (pdf_file)
      grid.draw(grobTree(venn_diagram))
      dev.off()
    }
  }
  
  ## DMR overlaps stacked bar chart
  {
    ### Load gene features files
    {
      print ("Loading gene feature files")
      #### Load TEs, remove partial matches, unknown TEs and filter for >500 bp
      {
        GFF_data_TEs <- rtracklayer::import (GFF_filepath_TEs, format = "gff")
        GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$type == "match"]
        GFF_data_TEs <- GFF_data_TEs [!is.na (GFF_data_TEs$Order)]
        GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$Order != "NA"]
        GFF_data_TEs <- GFF_data_TEs [seqnames (GFF_data_TEs) %in% seqnames (DMR_calling_chromosomes)]
        GFF_data_TEs <- unique (GFF_data_TEs)
        GFF_data_TEs_500bp <- GFF_data_TEs [width (GFF_data_TEs) > 500] # length (GFF_data_TEs) = 129220
        
        ##### Group TEs with < 10% of total number of TEs to "Other_TEs"
        non_na_idx <- !(GFF_data_TEs_500bp$Superfamily == "NA")
        GFF_data_TEs_500bp$Order [non_na_idx] <- paste (GFF_data_TEs_500bp$Order [non_na_idx], GFF_data_TEs_500bp$Superfamily [non_na_idx], sep = "_")
        Cutoff_threshold <- 0.1 * length (GFF_data_TEs_500bp)
        Cutoff_TE_elements <- table (GFF_data_TEs_500bp$Order) < Cutoff_threshold
        mcols (GFF_data_TEs_500bp)$Order [mcols (GFF_data_TEs_500bp)$Order %in% names (Cutoff_TE_elements [Cutoff_TE_elements])] <- "Other_TE"
        mcols (GFF_data_TEs_500bp)$type <- mcols (GFF_data_TEs_500bp)$Order
      }
      
      #### Load genic features and promoters
      {
        GFF_data_genes <- rtracklayer::import (GFF_filepath_genes, format = "gff")
        GFF_data_promoters <- GenomicRanges::flank (GFF_data_genes [GFF_data_genes$type == "mRNA"], width = promoter_width, start = TRUE)
        GFF_data_promoters$type <- "promoter"
      }
      
      GFF_all_gene_features <- unique (c (GFF_data_genes, GFF_data_promoters, GFF_data_TEs_500bp))
      #### Remove three_prime_UTR, five_prime_UTR and CDS
      GFF_all_gene_features <- GFF_all_gene_features [!GFF_all_gene_features$type %in% c ("three_prime_UTR", "five_prime_UTR", "CDS")]
      mcols (GFF_all_gene_features)$type <- droplevels (mcols (GFF_all_gene_features)$type)
    }
    
    Comparison_DMRs <- c ("SummerAmbientParents_vs_SpringAmbientParents",
                          "AutumnAmbientParents_vs_SpringAmbientParents",
                          "AutumnAmbientParents_vs_SummerAmbientParents")
    Contexts <- c ("CG", "CHG", "CHH", "Random")
    
    Outer_list <- list ()
    
    for (context in Contexts) {
      print (paste ("Working on", context))
      All_context_DMRs <- lapply (Reversed_combinations, function (DMR_pair) {
        if (file.exists (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
          qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
        } else if (context == "Random") {
          function_generate_random_regions (sample_number, Chromosome_lengths, region_width)
        }
      })
      names (All_context_DMRs) <- Reversed_combinations
      
      Plot_data <- lapply (Comparison_DMRs, function (comparison_DMRs) {
        DMR_overlaps <- IRanges::overlapsAny (query = GFF_all_gene_features, subject = All_context_DMRs [[comparison_DMRs]])
        DMR_no_overlaps <- countOverlaps (query = All_context_DMRs [[comparison_DMRs]], subject = GFF_all_gene_features) > 0
        All_overlaps <- table (GFF_all_gene_features  [DMR_overlaps]$type)
        All_overlaps <- as.data.frame (All_overlaps)
        All_overlaps <- rbind (All_overlaps, data.frame (Var1 = "Intergenic", Freq = sum (DMR_no_overlaps)))
        colnames (All_overlaps) <- c ("Feature", "Freq")
        return (All_overlaps)
      })
      names (Plot_data) <- Comparison_DMRs
      Outer_list [[context]] <- Plot_data
    }
    
    Contexts <- c ("CG", "CHG", "CHH")
    for (context in Contexts) {
      Outer_list [[context]]$Random <- Outer_list$Random [["SummerAmbientParents_vs_SpringAmbientParents"]]
      Out_data <- do.call (rbind, lapply (names(Outer_list[[context]]), function(name) {
        df <- Outer_list[[context]][[name]]
        df$Group <- name
        return(df)
      }))
      
      Out_data$Feature <- factor (Out_data$Feature,
                                  levels = c ("mRNA", "promoter", "LTR_Gypsy", "LTR_Copia", "TIR", "LINE", "Other_TE","Intergenic"))
      Feature_custom_labels <- c (
        "mRNA" = "Gene",
        "promoter" = "Promoter Region",
        "LTR_Gypsy" = "Gypsy LTR",
        "LTR_Copia" = "Copia LTR",
        "Other_TE" = "Other TE",
        "TIR" = "TIR",
        "LINE" = "LINE",
        "Intergenic" = "Intergenic Region"
      )
      Out_data$Group <- factor (Out_data$Group,
                                levels = c ( "SummerAmbientParents_vs_SpringAmbientParents", "AutumnAmbientParents_vs_SpringAmbientParents", "Random"))
      Group_custom_labels <- c (
        "Random" = "Random",
        "SummerAmbientParents_vs_SpringAmbientParents" = "Spring vs\nSummer", 
        "AutumnAmbientParents_vs_SpringAmbientParents" = "Spring vs\nAutumn"
      )
      Out_data <- Out_data [!is.na (Out_data$Group), ]
      
      #### Fisher test 
      {
        function_do_fisher_test <- function (frequency_table, group, feature) {
          print (feature)
          print (group)
          feature_random <- frequency_table [frequency_table$Group == "Random" & frequency_table$Feature == feature, "Freq"]
          not_feature_random <- sum (frequency_table [frequency_table$Group == "Random", "Freq"]) - feature_random
          feature_spring_summer <- frequency_table [frequency_table$Group == group & frequency_table$Feature == feature, "Freq"]
          not_feature_spring_summer <- sum (frequency_table [frequency_table$Group == group, "Freq"]) - feature_spring_summer
          
          counts <- matrix (c (feature_random, not_feature_random, feature_spring_summer, not_feature_spring_summer),
                            nrow = 2,
                            byrow = FALSE,
                            dimnames = list (
                              Feature = c (feature, paste0 ("Not", feature)),
                              comparison = c ("Random", group)
                            ))
          print (counts)
          fisher.test (counts)
        }
        
        function_do_fisher_test (Out_data, "SummerAmbientParents_vs_SpringAmbientParents", "TIR")
        function_do_fisher_test (Out_data, "AutumnAmbientParents_vs_SpringAmbientParents", "TIR")
        function_do_fisher_test (Out_data, "SummerAmbientParents_vs_SpringAmbientParents", "promoter")
        function_do_fisher_test (Out_data, "AutumnAmbientParents_vs_SpringAmbientParents", "promoter")
        
        function_do_fisher_test (Out_data, "SummerAmbientParents_vs_SpringAmbientParents", "LTR_Gypsy")
      }
      
      Out_plot <- ggplot (data = as.data.frame (Out_data), aes (x = Group, y = Freq, fill = Feature)) + 
        geom_bar (position = "fill", stat = "identity") + 
        scale_fill_viridis (discrete = TRUE, option = "F", labels = Feature_custom_labels) +
        scale_x_discrete (labels = Group_custom_labels) +
        labs (title = context, x = "DMR Comparisons", y = "Proportion") + 
        coord_flip () + 
        theme_bw () + 
        theme (axis.ticks.x = element_blank (),
               plot.title = element_text (size = 14, face = "bold", hjust = 0.5),
               axis.text.x = element_text (size = 10, hjust = 0.5, vjust = 0.5),
               axis.text.y = element_text (size = 10)
               #legend.position = "none
        )
      ggsave (plot = Out_plot, 
              filename = file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", paste0 ("DMR_GFF_overlap_parents_", context, ".png")), 
              height = 8, width = 16, units = "cm")
      ggsave (plot = Out_plot, 
              filename = file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", paste0 ("DMR_GFF_overlap_parents_", context, ".pdf")), 
              height = 8, width = 16, units = "cm")
    }
  }
  
  ## Heatmap of Parents at DMRs
  {
    All_individuals <- c (SpringAmbientParents, SummerAmbientParents, AutumnAmbientParents, AmbientOffspring)
    All_heatmap_data <- lapply (All_individuals, function (individual) {
      qread (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps", paste0 ("Heatmap_recomputed_ambient_par_and_ambient_off_", individual, ".qs")))
    })
    names (All_heatmap_data) <- All_individuals
    
    Contexts <- c ("CG", "CHG", "CHH")
    for (context in Contexts) {
      print (paste ("Working on", context))
      if (exists ("Heatmap_matrix")) {
        rm ("Heatmap_matrix")
      }
      
      for (individual in All_individuals) {
        if (!exists ("Heatmap_matrix")) {
          Heatmap_matrix_data <- cbind (as.character (seqnames (All_heatmap_data [[individual]] [[context]])), 
                                        start = start(All_heatmap_data [[individual]] [[context]]),                               # Start positions
                                        end = end(All_heatmap_data [[individual]] [[context]]),                                   # End positions
                                        mcols (All_heatmap_data [[individual]] [[context]] [,paste0 ("proportion", context)]))
          Heatmap_matrix_data$ID <- paste (Heatmap_matrix_data$`as.character(seqnames(All_heatmap_data[[individual]][[context]]))`, Heatmap_matrix_data$start, Heatmap_matrix_data$end, sep = "_")
          Heatmap_matrix_data <- Heatmap_matrix_data [, c ("ID", paste0 ("proportion", context))]
          Group <- find_group (individual, All_analysis_groups)
          colnames (Heatmap_matrix_data) <- c ("ID", paste (individual, Group, sep = "_"))
          Heatmap_matrix <- as.data.frame (Heatmap_matrix_data)
        } else {
          Heatmap_matrix_data <- cbind (as.character (seqnames (All_heatmap_data [[individual]] [[context]])), 
                                        start = start(All_heatmap_data [[individual]] [[context]]),                               # Start positions
                                        end = end(All_heatmap_data [[individual]] [[context]]),                                   # End positions
                                        mcols (All_heatmap_data [[individual]] [[context]] [,paste0 ("proportion", context)]))
          Heatmap_matrix_data$ID <- paste (Heatmap_matrix_data$`as.character(seqnames(All_heatmap_data[[individual]][[context]]))`, Heatmap_matrix_data$start, Heatmap_matrix_data$end, sep = "_")
          Heatmap_matrix_data <- Heatmap_matrix_data [, c ("ID", paste0 ("proportion", context))]
          Group <- find_group (individual, All_analysis_groups)
          colnames (Heatmap_matrix_data) <- c ("ID", paste (individual, Group, sep = "_"))
          Heatmap_matrix <- left_join (Heatmap_matrix, as.data.frame (Heatmap_matrix_data), by = "ID")
        }
      }
      
      Heatmap_matrix <- column_to_rownames (Heatmap_matrix, var = "ID")
      ### This step may not be "technically" proper way of handeling, but not sure what else to do
      Heatmap_matrix [is.na (Heatmap_matrix)] <- 0
      if (max (nrow (Heatmap_matrix)) > 10000) {
        Heatmap_matrix <- Heatmap_matrix [sample (nrow (Heatmap_matrix), 10000, replace = FALSE), ]
      }
      
      colnames (Heatmap_matrix) <- gsub ("Ambient", "", colnames (Heatmap_matrix))
      colnames (Heatmap_matrix) <- gsub ("Parents", "", colnames (Heatmap_matrix))
      
      colnames (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3",
                                      "Spring6","Summer1", "Summer2", "Autumn1", 
                                      "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
      
      # Heatmap_matrix <- Heatmap_matrix [, colnames (Heatmap_matrix) [grepl ("[12]", colnames (Heatmap_matrix))]] # Remove the progeny from this plot
      png (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("ALL_DMR_heatmap_recomputed_ambient_par_and_ambient_offspring_", context, ".png")), width = 400, height = 400)
      a <- pheatmap (as.matrix (Heatmap_matrix), 
                     main = context,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     legend = TRUE
      )
      dev.off()
      pdf (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("ALL_DMR_heatmap_recomputed_ambient_par_and_ambient_offspring_", context, ".pdf")))
      a <- pheatmap (as.matrix (Heatmap_matrix), 
                     main = context,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     legend = TRUE
      )
      dev.off()
    }
  }
}

# Figure 3
{
  ## Ontology Plots
  {
    ### Get Spring DMRs
    {
      context <- "CHH"
      Comparison_DMRs <- c ("SummerAmbientParents_vs_SpringAmbientParents",
                            "AutumnAmbientParents_vs_SpringAmbientParents",
                            "AutumnAmbientParents_vs_SummerAmbientParents")
      Overlapping_genes_file <- paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/Seasons_", context, "_overlapping_genes.RData")
      genes_of_interest <- Overlapping_genes_file
      enrichment_RData <- paste0 (GO_out_folder, "Enrichment_results_", basename (genes_of_interest))
      
      All_context_DMRs <- lapply (Comparison_DMRs, function (DMR_pair) {
        if (file.exists (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/", paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
          qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/", paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
        }
      })
      names (All_context_DMRs) <- Comparison_DMRs
      
      Spring_DMRs <- c (All_context_DMRs$SummerAmbientParents_vs_SpringAmbientParents, All_context_DMRs$AutumnAmbientParents_vs_SpringAmbientParents)
      Spring_DMRs <- unique (Spring_DMRs)
    }
    
    ### Overlap with Genes and Promotors
    {
      Overlapping_promotors_file <- paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/Seasons_", context, "_overlapping_promoters.RData")
      Overlaps_promotor <- GFF_data_promoters [overlapsAny (GFF_data_promoters, Spring_DMRs)]
      save (Overlaps_promotor, file = Overlapping_promotors_file)
      
      Overlapping_genes_file <- paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/Seasons_", context, "_overlapping_genes.RData")
      Overlaps_genes <- GFF_data_genes [overlapsAny (GFF_data_genes, Spring_DMRs)]
      save (Overlaps_genes, file = Overlapping_genes_file)
      
      genes_of_interest <- Overlapping_promotors_file
      enrichment_RData <- paste0 (GO_out_folder, "Enrichment_results_", basename (genes_of_interest))
      enrichment <- GO_enrichment(
        genes_of_interest = Overlapping_promotors_file,
        GO_database = GO_database, 
        gff_file = gff_file,
        out_folder = GO_out_folder,
        useInfo = "all",
        min_sig_genes = 10,
        fisher_threshold = 0.05)
      
      GO_plots(enrichment_RData = enrichment_RData,
               min_sig_genes = 10,
               fisher_threshold = 0.02,
               plot_title = paste (context, "Promoter Ontology"),
               max_y_lim = 150,
               plot_height = 16,
               plot_width = 43.5,
               out_folder = "/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe")
      
      genes_of_interest <- Overlapping_genes_file
      enrichment_RData <- paste0 (GO_out_folder, "Enrichment_results_", basename (genes_of_interest))
      enrichment <- GO_enrichment(
        genes_of_interest = Overlapping_genes_file,
        GO_database = GO_database, 
        gff_file = gff_file,
        out_folder = GO_out_folder,
        useInfo = "all",
        min_sig_genes = 10,
        fisher_threshold = 0.05)
      
      GO_plots(enrichment_RData = enrichment_RData,
               min_sig_genes = 10,
               fisher_threshold = 0.02,
               plot_width = 43.5,
               plot_title = paste (context, "Gene Ontology"),
               out_folder = "/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe")
    }
  }
  
  ## TE PCA, Kmeans and heatmap
  {
    All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
    Heatmap_matrix <- qread (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/TE_500bp_Ambient_par_offspring_", context, ".qs"))
    Heatmap_matrix <- Heatmap_matrix [!duplicated (Heatmap_matrix$ID), ] # remove duplicated data
    rownames (Heatmap_matrix) <- NULL
    Heatmap_matrix <- column_to_rownames (Heatmap_matrix, var = "ID")
    context <- "CHH"
    Heatmap_matrix <- Heatmap_matrix [, paste (All_individuals, context, sep = "_")] # Keep only datasets with the individuals of interest
    Heatmap_matrix [is.na (Heatmap_matrix)] <- 0 # Set NA values in heatmap arbitrarily to 0
    names (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3", "Spring6", "Summer1", "Summer2", "Autumn1", "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
    
    Kmean <- kmeans (Heatmap_matrix, centers = 2)
    table (Kmean$cluster)
    
    ## Make PCA (need to load GFF_data_TEs_500bp)
    {
      Out_data <- prcomp (t(Heatmap_matrix))
      Prop_var <- summary(Out_data)$importance[2, ]
      Out_data <- Out_data$rotation [, 1:2] # get only first 2 principle components
      TE_classes <- data.frame (TE = paste (seqnames (GFF_data_TEs_500bp), start (GFF_data_TEs_500bp), end (GFF_data_TEs_500bp),sep = "_"),
                                Class = paste (GFF_data_TEs_500bp$Order, GFF_data_TEs_500bp$Superfamily, sep = "_"))
      TE_classes <- unique (TE_classes)
      TE_classes <- TE_classes [!duplicated (TE_classes$TE), ]
      ### Get classes of the TEs
      Out_data <- as.data.frame (Out_data)
      Out_data <- rownames_to_column (Out_data, var = "TE")
      Out_data <- left_join (Out_data, TE_classes, by = "TE")
      Out_data <- column_to_rownames (Out_data, var = "TE")
      
      Out_data$Class <- factor (Out_data$Class,
                                levels = c ("Other_TE_NA", "LTR_Copia_Copia", "LTR_Gypsy_Gypsy", "TIR_NA", "LINE_NA"),
                                labels = c ("Other_TE", "LTR_Copia", "LTR_Gypsy", "TIR", "LINE"))
      
      Outplot <- ggplot (data = Out_data, aes (x = PC1, y = PC2, color = Class)) +
        xlab (paste0 ("PC1 (", signif (Prop_var [1], 2), ")")) + 
        ylab (paste0 ("PC2 (", signif (Prop_var [2], 2), ")")) + 
        scale_color_viridis_d (option = "rocket") +
        geom_point (size = 0.01) + 
        # geom_hline(yintercept = - 0.03, linetype = "dashed") + 
        geom_density_2d (alpha = 1) + 
        theme_bw () +
        theme (plot.title = element_text (size = 16, face = "bold", hjust = 0.5),
               axis.text.x = element_text (size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
               axis.text.y = element_text (size = 12)
        )
      ggsave (plot = Outplot, filename = "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/TEs_500bp_PCA_ind.png", 
              height = 10, width = 10, units = "cm")
      ggsave (plot = Outplot, filename = "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/TEs_500bp_PCA_ind.pdf", 
              height = 10, width = 10, units = "cm")
      
      #### Plot using kmeans clusters
      {
        Out_data <- merge (Out_data, as.data.frame(Kmean$cluster), by = "row.names", all.x = TRUE)
        colnames (Out_data) [ncol(Out_data)] <- "Kmean_cluster"
        Out_data$Kmean_cluster <- as.character (Out_data$Kmean_cluster)
        
        Outplot <- ggplot (data = Out_data, aes (x = PC1, y = PC2, color = Kmean_cluster)) +
          xlab (paste0 ("PC1 (", signif (Prop_var [1], 2), ")")) + 
          ylab (paste0 ("PC2 (", signif (Prop_var [2], 2), ")")) + 
          scale_color_viridis_d (option = "rocket") +
          geom_point (size = 0.01) + 
          # geom_hline(yintercept = - 0.03, linetype = "dashed") + 
          geom_density_2d (alpha = 1) + 
          theme_bw () +
          theme (plot.title = element_text (size = 16, face = "bold", hjust = 0.5),
                 axis.text.x = element_text (size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
                 axis.text.y = element_text (size = 12)
          )
        
        ggsave (plot = Outplot, filename = "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/TEs_500bp_PCA_ind_kmeans.png", 
                height = 10, width = 10, units = "cm")
        ggsave (plot = Outplot, filename = "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/TEs_500bp_PCA_ind_kmeans.pdf", 
                height = 10, width = 10, units = "cm")
      }
    }
    
    ## Make main heatmap
    {
      if (nrow (Heatmap_matrix) > 50000) {
        Heatmap_matrix_sample <- Heatmap_matrix [sample (1:nrow (Heatmap_matrix), 50000, replace = FALSE),]
      } else {
        Heatmap_matrix_sample <- Heatmap_matrix
      }
      
      split_rn <- do.call (rbind, strsplit (row.names (Heatmap_matrix_sample), "_"))
      colnames (split_rn) <- c("prefix", "chr", "start", "end")
      rn_info <- data.frame(
        chr = paste0(split_rn[,1], "_", split_rn[,2]),
        start = as.numeric(split_rn[,3]),
        end = as.numeric(split_rn[,4]),
        row.names = paste(split_rn[,1], split_rn[,2], split_rn[,3], split_rn[,4], sep = "_")
      )
      
      ## Slow function ()
      matches <- mapply (function(chr, start, end) {
        which (seqnames(GFF_data_TEs_500bp) == chr & start(GFF_data_TEs_500bp) == start & end(GFF_data_TEs_500bp) == end)
      }, 
      rn_info$chr, rn_info$start, rn_info$end)
      
      # Filter(function(x) length(x) > 1, matches) # sometimes >1 matches can be found - check what they are and why
      matches <- lapply(matches, function(x) x[1]) # get rid of instances where there are more than one matches
      match_indices <- unlist (matches)
      match_indices <- unique (match_indices)
      order_values <- mcols (GFF_data_TEs_500bp)$Order [match_indices]
      unique_orders <- c ("Other_TE", "LTR_Copia", "LTR_Gypsy", "TIR", "LINE")
      custom_labels <- c(
        "Other_TE"   = "Other TEs",
        "LTR_Copia"  = "Copia LTR",
        "LTR_Gypsy"  = "Gypsy LTR",
        "TIR"        = "TIR",
        "LINE"       = "LINE"
      )
      
      annotation_row <- data.frame("TE Order" = factor(custom_labels[as.character(order_values)],
                                                       levels = custom_labels))
      rownames(annotation_row) <- rownames(Heatmap_matrix_sample)
      
      # Create colors using the custom labels
      ann_colors <- list(
        "TE Order" = setNames(
          viridisLite::viridis(length(custom_labels), option = "rocket"),
          custom_labels
        )
      )
      rownames(annotation_row) <- rownames(Heatmap_matrix_sample)
      png (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/TE_500bp_Ambient_par_offspring_50000_", context,".png"), width = 500, height = 500)
      a <- pheatmap (Heatmap_matrix_sample, 
                     main = paste (context, "Methylation in TEs"),
                     annotation_row = annotation_row,
                     annotation_colors = ann_colors,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     legend = TRUE
      )
      dev.off ()
      pdf (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/TE_500bp_Ambient_par_offspring_50000_", context,".pdf"))
      a <- pheatmap (Heatmap_matrix_sample, 
                     main = paste (context, "Methylation in TEs"),
                     annotation_row = annotation_row,
                     annotation_colors = ann_colors,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     legend = TRUE
      )
      dev.off ()
    }
    
    ## Make small cluster heatmap
    {
      All_individuals <- c (All_analysis_groups$SpringAmbientParents, All_analysis_groups$SummerAmbientParents, All_analysis_groups$AutumnAmbientParents, All_analysis_groups$AmbientOffspring)
      Heatmap_matrix <- qread (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/TE_500bp_Ambient_par_offspring_", context, ".qs"))
      Heatmap_matrix <- Heatmap_matrix [!duplicated (Heatmap_matrix$ID), ] # remove duplicated data
      rownames (Heatmap_matrix) <- NULL
      Heatmap_matrix <- column_to_rownames (Heatmap_matrix, var = "ID")
      Heatmap_matrix <- Heatmap_matrix [, paste (All_individuals, context, sep = "_")] # Keep only datasets with the individuals of interest
      Heatmap_matrix [is.na (Heatmap_matrix)] <- 0 # Set NA values in heatmap arbitrarily to 0
      names (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3", "Spring6", "Summer1", "Summer2", "Autumn1", "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
      
      Out_data <- prcomp (t(Heatmap_matrix))
      Prop_var <- summary(Out_data)$importance[2, ]
      Out_data <- Out_data$rotation [, 1:2] # get only first 2 principle components
      TE_classes <- data.frame (TE = paste (seqnames (GFF_data_TEs_500bp), start (GFF_data_TEs_500bp), end (GFF_data_TEs_500bp),sep = "_"),
                                Class = paste (GFF_data_TEs_500bp$Order, GFF_data_TEs_500bp$Superfamily, sep = "_"))
      TE_classes <- unique (TE_classes)
      TE_classes <- TE_classes [!duplicated (TE_classes$TE), ]
      ### Get classes of the TEs
      Out_data <- as.data.frame (Out_data)
      Out_data <- rownames_to_column (Out_data, var = "TE")
      Out_data <- left_join (Out_data, TE_classes, by = "TE")
      Out_data <- column_to_rownames (Out_data, var = "TE")
      
      Out_data$Class <- factor (Out_data$Class,
                                levels = c ("Other_TE_NA", "LTR_Copia_Copia", "LTR_Gypsy_Gypsy", "TIR_NA", "LINE_NA"),
                                labels = c ("Other_TE", "LTR_Copia", "LTR_Gypsy", "TIR", "LINE"))
      
      #### Full matrix
      {
        
      }
      
      #### Small cluster 
      {
        Small_cluster <- Out_data [Out_data$PC2 < -0.03, ]
        Counts <- cbind (as.matrix (table (Out_data$Class)), as.matrix (table (Small_cluster$Class)))
        colnames (Counts) <- c ("All_TEs", "Cluster_TEs")
        chisq.test (Counts) # X-squared = 11.921, df = 4, p-value = 0.01795
        
        Small_cluster_heatmap <- as.matrix (Heatmap_matrix [row.names (Small_cluster),])
        annotation_df <- data.frame (Class = Small_cluster$Class)
        rownames(annotation_df) <- rownames(Small_cluster_heatmap)
        png ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots/Small_cluster_heatmap.png", width = 1000, height = 800, res = 150) 
        d <- pheatmap (Small_cluster_heatmap, 
                       main = "Small Cluster",
                       color = colorRampPalette (c ("white", "red"))(100),
                       annotation_row = annotation_df ,
                       show_rownames = FALSE,
                       show_colnames = TRUE,
                       legend = TRUE
        )
        dev.off()
      }
    }
  }
  
  ## Most/least methylated TEs cluster mapability (violin plot)
  {
    ## Read first 24 lines (corresponding to first 12 chromosomes)
    
    context <- "CHH"
    All_individuals <- c (All_analysis_groups$SpringAmbientParents, All_analysis_groups$SummerAmbientParents, All_analysis_groups$AutumnAmbientParents, All_analysis_groups$AmbientOffspring)
    All_TE_500bp_mapability <- qread ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/All_TE_500bp_mappability.qs")
    Heatmap_matrix <- qread (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/TE_500bp_Ambient_par_offspring_CHH.qs"))
    Heatmap_matrix <- Heatmap_matrix [!duplicated (Heatmap_matrix$ID), ] # remove duplicated data
    rownames (Heatmap_matrix) <- NULL
    Heatmap_matrix <- column_to_rownames (Heatmap_matrix, var = "ID")
    Heatmap_matrix <- Heatmap_matrix [, paste (All_individuals, context, sep = "_")] # Keep only datasets with the individuals of interest
    Heatmap_matrix [is.na (Heatmap_matrix)] <- 0 # Set NA values in heatmap arbitrarily to 0
    names (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3", "Spring6", "Summer1", "Summer2", "Autumn1", "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
    
    Heatmap_matrix_top <- Heatmap_matrix
    Heatmap_matrix_top$Spring_mean <- rowMeans (Heatmap_matrix_top [grep ("^Spring", colnames (Heatmap_matrix_top))])
    Heatmap_matrix_top$Season_mean <- rowMeans (Heatmap_matrix_top [grep ("^(Summer|Autumn)", colnames (Heatmap_matrix_top))])
    Heatmap_matrix_top$Season_difference <- Heatmap_matrix_top$Spring_mean - Heatmap_matrix_top$Season_mean
    Heatmap_matrix_top <- Heatmap_matrix_top [order (-Heatmap_matrix_top$Season_difference), ]
    Heatmap_matrix_top <- Heatmap_matrix_top [1:round (nrow (Heatmap_matrix_top) / 10), ]
    
    ### Generate Mappability_list - mappability of all loci in genome
    {
      ## Read first 24 lines (corresponding to first 12 chromosomes)
      Mapability <- readLines (Oak_mapability_filepath, n = 24) # Don't print to screen - probably will crash R
      
      Mapability_list <- list ()
      for (chromosome in seq (2,24,2)) {
        print (paste ("Working on Chrom", Mapability [[chromosome - 1]]))
        Mapability_list [Mapability [[chromosome - 1]]] <- strsplit (Mapability [[chromosome]], split = " ")
      }
      names (Mapability_list) <- sapply (names (Mapability_list), second_to_last_char)
      Mapability_list <- lapply (Mapability_list, as.numeric)
      qsave (Mapability_list, file = "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Mapability_list_whole_genome.qs")
      Mapability_list <- qread ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Mapability_list_whole_genome.qs")
    }
    
    ### Get mappability of most methylated TEs
    {
      Methylated_TEs <- list ()
      for (analysis in c ("Top", "Bottom")) {
        print (paste ("Working on ", analysis))
        Heatmap_matrix_subset <- Heatmap_matrix
        Heatmap_matrix_subset$Spring_mean <- rowMeans (Heatmap_matrix_subset [grep ("^Spring", colnames (Heatmap_matrix_subset))])
        Heatmap_matrix_subset$Season_mean <- rowMeans (Heatmap_matrix_subset [grep ("^(Summer|Autumn)", colnames (Heatmap_matrix_subset))])
        Heatmap_matrix_subset$Season_difference <- Heatmap_matrix_subset$Spring_mean - Heatmap_matrix_subset$Season_mean
        if (analysis == "Top") {
          Heatmap_matrix_subset <- Heatmap_matrix_subset [order (-Heatmap_matrix_subset$Season_difference), ]
        } else if (analysis == "Bottom") {
          Heatmap_matrix_subset <- Heatmap_matrix_subset [order (Heatmap_matrix_subset$Season_difference), ]
        }
        Heatmap_matrix_subset <- Heatmap_matrix_subset [1:round (nrow (Heatmap_matrix_subset) / 10), ]
        TE_data <- function_get_TE_from_granges (Heatmap_matrix_subset, GFF_data_TEs_500bp)
        
        TE_500bp_mappability <- list ()
        for (TE in 1:length (TE_data)) {
          print (paste (TE, "of", length (TE_data)))
          TE_500bp_mappability [[TE]] <- Mapability_list [as.character (seqnames (TE_data [TE]))][[1]] [start (TE_data [TE]) : end (TE_data [TE])]
        }
        means_list <- lapply (TE_500bp_mappability, function (x) median (round (1/x), na.rm = TRUE))
        means_list <- unlist (means_list)
        means_list <- data.frame (Frequency = means_list)
        Methylated_TEs [analysis] <- means_list
        
        Out_plot <- ggplot (data = means_list, aes (x = Frequency)) + 
          geom_histogram () +
          geom_vline (xintercept = median (means_list$Frequency)) + 
          annotate ("text", x = median (means_list$Frequency), y = 1000, hjust = 0,
                    label = paste ("Median:", round (median (means_list$Frequency), 2))) +
          ggtitle (paste (analysis, "Differentially Methylated TEs")) +
          theme_bw () + 
          scale_x_continuous (limits = c (0, 500)) +
          scale_y_continuous (limits = c (0, 100000)) +
          ylab (NULL) +
          theme (axis.ticks.x = element_blank (),
                 plot.title = element_text (size = 14, face = "bold", hjust = 0.5),
                 axis.text.x = element_text (size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
                 axis.text.y = element_text (size = 10))
        
        ggsave (paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all/", analysis, "_methylated_Spring_TEs.png"), 
                height = 15, width = 15, units = "cm")
      }
      
      #### Plot violin plot of Methylated TEs
      {
        Violin_plot_data <- data.frame(
          value = c (Methylated_TEs$Top, Methylated_TEs$Bottom),
          group = factor (c(
            rep ("Top", length(Methylated_TEs$Top)),
            rep ("Bottom", length(Methylated_TEs$Bottom))
          ))
        )
        
        Violin_plot <- ggplot (data = Violin_plot_data, aes (x = group, y = value, fill = group)) +
          geom_violin () +
          scale_x_discrete (labels = c ("Least Methylated 10%", "Most Methylated 10%")) +
          scale_y_continuous (limits = c (0, 100)) +
          scale_fill_manual(
            values = c("Top" = "#1f77b4", "Bottom" = "#ff7f0e"),
            labels = c("Top" = "Most Methylated 10%", "Bottom" = "Least Methylated 10%")
          ) +
          stat_summary (fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
          stat_summary (fun = median,geom = "text",aes(label = round(..y.., 2)),
                        vjust = -1.5,fontface = "bold",color = "black") +
          labs (x = NULL, y = "Frequency", fill = "TE Methylation") +
          theme_bw () + 
          theme (axis.text.x = element_text (angle = 90, hjust = 0.5, vjust = 0.5),
                 axis.ticks.x = element_blank ()
          )
        ggsave ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all/TE_Frequency.png", 
                height = 15, width = 15, units = "cm")
        ggsave ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all/TE_Frequency.pdf",
                height = 20, width = 20, units = "cm")
      }
    }
  }
}

# Figure 4
{
  ## Parents Offspring DMRs count bar chart 
  {
    Par_vs_Off <- "AmbientParents_vs_AmbientOffspring"
    Par_vs_Off_DMRsCG <- qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakProgeny_vs_par/analysis/objects/Pairwise_Comparisons",
                                           paste0 ("DMRs_CG", "_", Par_vs_Off, ".qs")))
    Par_vs_Off_DMRsCHG <- qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakProgeny_vs_par/analysis/objects/Pairwise_Comparisons",
                                            paste0 ("DMRs_CHG", "_", Par_vs_Off, ".qs")))
    Par_vs_Off_DMRsCHH <- qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakProgeny_vs_par/analysis/objects/Pairwise_Comparisons",
                                            paste0 ("DMRs_CHH", "_", Par_vs_Off, ".qs")))
    
    Plot_dataframe <- data.frame (CG = table (Par_vs_Off_DMRsCG$regionType),
                                  CHG = table (Par_vs_Off_DMRsCHG$regionType),
                                  CHH = table (Par_vs_Off_DMRsCHH$regionType))
    
    rownames (Plot_dataframe) <- c ("gain", "loss")
    Plot_dataframe <- Plot_dataframe [, c ("CG.Freq", "CHG.Freq", "CHH.Freq")]
    colnames (Plot_dataframe) <- c ("CG", "CHG", "CHH")
    Plot_dataframe$category <- rownames(Plot_dataframe)
    
    Plot_dataframe <- pivot_longer(Plot_dataframe, cols = CG:CHH, names_to = "context", values_to = "count")
    Plot_dataframe [Plot_dataframe$category == "loss", "count"] <- Plot_dataframe [Plot_dataframe$category == "loss", "count"] * -1
    Plot_dataframe$count <- Plot_dataframe$count * -1
    
    Out_plot <- ggplot (data = Plot_dataframe, aes (x = context, y = count, fill = category)) +
      geom_bar (stat = "identity", position = "stack") + 
      labs (title = "Progeny vs Parents",
            x = NULL,
            y = NULL) + 
      lims (y = c (-max (abs (Plot_dataframe$count)) * 1.1, max (abs (Plot_dataframe$count)) * 1.1)) + 
      scale_fill_viridis (discrete = TRUE, option = "inferno", begin = 0.2, end = 0.8) + 
      # scale_x_discrete (labels = c ("Spring\nvs\nSummer", "Spring\nvs\nAutumn", "Summer\nvs\nAutumn")) + 
      theme_bw () + 
      theme (# axis.text.x = element_blank (),
        axis.ticks.x = element_blank (),
        legend.position = "none"
      )
    
    ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_par_offs", ".png"),  height = 10, width = 7, units = "cm")
    ggsave (Out_plot, filename = paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/", "DMRs_barplot_data_par_offs", ".pdf"),  height = 10, width = 7, units = "cm")
  }
  
  ## DMR overlaps stacked bar chart
  {
    ### Load GFF file as GRanges object and filter
    {
      GFF_data_genes <- rtracklayer::import (Gene_annotation_file, format = "gff")
      GFF_data_mRNA <- GFF_data_genes [GFF_data_genes$type == "mRNA"]
      GFF_data_mRNA$score [is.na (GFF_data_mRNA$score)] <- 0 # for testing writing to bed file
      rtracklayer::export (GFF_data_mRNA, format = "BED", "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Par_offspring_deeptools/allmRNA.bed")
      GFF_data_genes <- GFF_data_genes [GFF_data_genes$type == "mRNA"]
      
      ## Load TEs and remove partial matches and unknown TEs
      GFF_data_TEs <- rtracklayer::import (GFF_filepath_TEs, format = "gff")
      GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$type == "match"]
      GFF_data_TEs <- GFF_data_TEs [!is.na (GFF_data_TEs$Order)]
      GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$Order != "NA"]
      GFF_data_TEs <- GFF_data_TEs [seqnames (GFF_data_TEs) %in% seqnames (DMR_calling_chromosomes)]
      GFF_data_TEs <- unique (GFF_data_TEs)
      # Filter for TEs >500bp in length and where order 
      GFF_data_TEs_500bp <- GFF_data_TEs [width (GFF_data_TEs) > 500] # length (GFF_data_TEs) = 129220
      non_na_idx <- !(GFF_data_TEs_500bp$Superfamily == "NA")
      GFF_data_TEs_500bp$Order [non_na_idx] <- paste (GFF_data_TEs_500bp$Order [non_na_idx], GFF_data_TEs_500bp$Superfamily [non_na_idx], sep = "_")
      Cutoff_threshold <- 0.1 * length (GFF_data_TEs_500bp)
      Cutoff_TE_elements <- table (GFF_data_TEs_500bp$Order) < Cutoff_threshold
      mcols (GFF_data_TEs_500bp)$Order [mcols (GFF_data_TEs_500bp)$Order %in% names (Cutoff_TE_elements [Cutoff_TE_elements])] <- "Other_TE"
      mcols (GFF_data_TEs_500bp)$type <- mcols (GFF_data_TEs_500bp)$Order
      
      GFF_data_promoters <- GenomicRanges::flank (GFF_data_genes [GFF_data_genes$type == "mRNA"], width = promoter_width, start = TRUE)
      GFF_data_promoters$type <- "promoter"
      
      GFF_combined <- c (GFF_data_genes, GFF_data_TEs_500bp, GFF_data_promoters)
    }
    
    ### Load DMR data comparing parent and offspring
    {
      Contexts <- c ("CG", "CHG", "CHH")
      gain_loss <- "loss"
      All_DMRs <- lapply (Contexts, function (context) {
        if (file.exists (file.path (Working_dir, "OakProgeny_vs_par", paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", "AmbientParents_vs_AmbientOffspring", ".qs")))) {
          a <- qread (file.path (Working_dir, "OakProgeny_vs_par", paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", "AmbientParents_vs_AmbientOffspring", ".qs")))
          a <- a [a$regionType == gain_loss]
        }
      })
      names (All_DMRs) <- Contexts
    }
    
    ### Tabulate overlaps between genes, TEs and non-overlapping regions
    {
      Inner_list <- list ()
      Contexts <- c (Contexts, "Random")
      Random_overlaps_list <- function_generate_random_regions_2 (sample_number, Chromosome_lengths, binSize, as.character (seqnames (DMR_calling_chromosomes)))
      All_DMRs [["Random"]] <- Random_overlaps_list
      for (context in Contexts) {
        print (paste ("working on:", context))
        gene_overlaps <- subsetByOverlaps (GFF_combined, All_DMRs [[context]], ignore.strand = TRUE)
        gene_overlaps_table <- table (gene_overlaps$type)
        gene_overlaps_table <- as.data.frame (gene_overlaps_table)
        names (gene_overlaps_table) <- c ("Feature", "Freq")
        
        # Does not overlap with gene features
        no_overlaps <- All_DMRs [[context]] [!overlapsAny (All_DMRs [[context]], GFF_combined, ignore.strand = TRUE)]
        no_overlaps_table <- data.frame (Feature = "non_overlapping", Freq = length (no_overlaps))
        gene_overlaps_table <- rbind (gene_overlaps_table, no_overlaps_table)
        
        Inner_list [[context]] <- gene_overlaps_table
      }
    }
    
    ### Plot
    {
      dir.create ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots", recursive = TRUE, showWarnings = FALSE)
      if (exists ("Out_data")) {
        rm (Out_data)
      }
      Out_data2 <- rbind( # Need to separately do Out_data and Out_data2
        cbind (Inner_list$Random, Context = paste0 ("Random_", gain_loss)),
        cbind (Inner_list$CG, Context = paste0 ("CG_", gain_loss)),
        cbind (Inner_list$CHG, Context = paste0 ("CHG_", gain_loss)),
        cbind (Inner_list$CHH, Context = paste0 ("CHH_", gain_loss))
      )
      Out_data <- rbind (Out_data1, Out_data2) # Need to separately do Out_data and Out_data2
      Out_data <- Out_data [!Out_data$Feature %in% c ("CDS", "five_prime_UTR", "three_prime_UTR"), ]
      Out_data$Context [Out_data$Context == "Random_gain"] <- "Random"
      Out_data$Context <- factor (Out_data$Context, levels = c("CG_gain", "CG_loss", "CHG_gain", "CHG_loss", "Random"))
      Out_data <- droplevels(Out_data)
      Out_data$Feature <- factor (Out_data$Feature, levels = c ("mRNA", "promoter", "LTR_Gypsy", "LTR_Copia", "TIR", "LINE", "Other_TE", "non_overlapping"))
      Feature_custom_labels <- c (
        "mRNA" = "Gene",
        "promoter" = "Promoter Region",
        "LTR_Gypsy" = "Gypsy LTR",
        "LTR_Copia" = "Copia LTR",
        "Other_TE" = "Other TE",
        "TIR" = "TIR",
        "LINE" = "LINE",
        "non_overlapping" = "Intergenic Region"
      )
      Out_data <- Out_data [!is.na (Out_data$Context), ]
      levels(Out_data$Context) <- gsub("_", " ", levels(Out_data$Context))
      Out_plot <- ggplot (data = Out_data, aes (x = Context, y = Freq, fill = Feature)) +
        geom_bar (position = "fill", stat = "identity") + 
        scale_fill_viridis (discrete = TRUE, option = "F", labels = Feature_custom_labels) +
        # scale_x_discrete (labels = c ("Random", "Spring\nvs\nSummer", "Spring\nvs\nAutumn", "Summer\nvs\nAutumn")) + 
        labs (title = NULL, x = "DMR Comparisons", y = "Proportion") + 
        coord_flip () + 
        theme_bw () + 
        theme (axis.ticks.x = element_blank (),
               #legend.position = "none",
               plot.title = element_text (hjust = 0.5)) + 
        ggtitle (paste ("DMRs Genic Feature Overlaps"))
      ggsave (plot = Out_plot, 
              filename = file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots", paste0 ("Parents_Offspring_DMR_Gene_TE_overlaps_CG_CHG", ".png")), 
              height = 12, width = 30, units = "cm")
      ggsave (plot = Out_plot, 
              filename = file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plots", paste0 ("Parents_Offspring_DMR_Gene_TE_overlaps_CG_CHG", ".pdf")), 
              height = 12, width = 30, units = "cm")
    }
  }
  
  ## Parent Offspring Heatmaps
  {
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
      All_DMRs <- lapply (Contexts, function (context) {
        if (file.exists (file.path (Working_dir, "OakProgeny_vs_par", paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", "AmbientParents_vs_AmbientOffspring", ".qs")))) {
          qread (file.path (Working_dir, "OakProgeny_vs_par", paste0 ("analysis/objects/Pairwise_Comparisons/DMRs_", context, "_", "AmbientParents_vs_AmbientOffspring", ".qs")))
        }
      })
      names (All_DMRs) <- Contexts
    }
    
    ## Recalculate level of methylation for individual
    {
      # Parent_Offsrping_Inds <- c ("P63872", "P63875", "P86734", "Aa6387J", "Aa8673J","Aa6387S", "Aa8673S")
      # for (individual in Parent_Offsrping_Inds) {
      #   CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports/", paste0 ("CX_report_", individual, ".qs" )))
      #   Heatmap_data <- lapply (Contexts, function (context) {
      #     print (paste ("working on", context))
      #     DMRcaller::analyseReadsInsideRegionsForCondition (All_DMRs [[context]], CX_report, context = context)
      #   })
      #   names (Heatmap_data) <- Contexts
      #   dir.create ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Heatmap_data")
      #   qsave (Heatmap_data, file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Heatmap_data", 
      #                                   paste0 ("Heatmap_Parent_Off_", individual, ".qs")))
      # }
    }
    
    ## Load data and plot heatmap
    {
      Parent_Offsrping_Inds <- c ("P63872", "P63875", "P86734", "Aa6387J", "Aa8673J","Aa6387S", "Aa8673S")
      All_heatmap_data <- lapply (Parent_Offsrping_Inds, function (individual) {
        qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Heatmap_data", 
                          paste0 ("Heatmap_Parent_Off_", individual, ".qs")))
      })
      names (All_heatmap_data) <- Parent_Offsrping_Inds
      
      Contexts <- c ("CG", "CHG", "CHH")
      for (context in Contexts) {
        if (exists ("Heatmap_matrix")) {
          rm ("Heatmap_matrix")
        }
        for (individual in Parent_Offsrping_Inds) {
          if (!exists ("Heatmap_matrix")) {
            Heatmap_matrix_data <- cbind (as.character (seqnames (All_heatmap_data [[individual]] [[context]])), 
                                          start = start(All_heatmap_data [[individual]] [[context]]),                               # Start positions
                                          end = end(All_heatmap_data [[individual]] [[context]]),                                   # End positions
                                          mcols (All_heatmap_data [[individual]] [[context]] [,paste0 ("proportion", context)]))
            Heatmap_matrix_data$ID <- paste (Heatmap_matrix_data$`as.character(seqnames(All_heatmap_data[[individual]][[context]]))`, Heatmap_matrix_data$start, Heatmap_matrix_data$end, sep = "_")
            Heatmap_matrix_data <- Heatmap_matrix_data [, c ("ID", paste0 ("proportion", context))]
            Group <- find_group (individual, All_analysis_groups)
            colnames (Heatmap_matrix_data) <- c ("ID", paste (individual, Group, sep = "_"))
            Heatmap_matrix <- as.data.frame (Heatmap_matrix_data)
          } else {
            Heatmap_matrix_data <- cbind (as.character (seqnames (All_heatmap_data [[individual]] [[context]])), 
                                          start = start(All_heatmap_data [[individual]] [[context]]),                               # Start positions
                                          end = end(All_heatmap_data [[individual]] [[context]]),                                   # End positions
                                          mcols (All_heatmap_data [[individual]] [[context]] [,paste0 ("proportion", context)]))
            Heatmap_matrix_data$ID <- paste (Heatmap_matrix_data$`as.character(seqnames(All_heatmap_data[[individual]][[context]]))`, Heatmap_matrix_data$start, Heatmap_matrix_data$end, sep = "_")
            Heatmap_matrix_data <- Heatmap_matrix_data [, c ("ID", paste0 ("proportion", context))]
            Group <- find_group (individual, All_analysis_groups)
            colnames (Heatmap_matrix_data) <- c ("ID", paste (individual, Group, sep = "_"))
            Heatmap_matrix <- left_join (Heatmap_matrix, as.data.frame (Heatmap_matrix_data), by = "ID")
          }
        }
        
        Heatmap_matrix <- column_to_rownames (Heatmap_matrix, var = "ID")
        ### This step may not be "technically" proper way of handeling, but not sure what else to do
        Heatmap_matrix [is.na (Heatmap_matrix)] <- 0
        if (max (nrow (Heatmap_matrix)) > 10000) {
          Heatmap_matrix <- Heatmap_matrix [sample (nrow (Heatmap_matrix), 10000, replace = FALSE), ]
        }
        
        colnames (Heatmap_matrix) <- gsub ("Ambient", "", colnames (Heatmap_matrix))
        colnames (Heatmap_matrix) <- gsub ("Parents", "", colnames (Heatmap_matrix))
        
        names (Heatmap_matrix)
        colnames (Heatmap_matrix) <- c ("Progeny1A", "Progeny1B", "Progeny2", 
                                        "Summer1", "Summer2", 
                                        "Autumn1", "Autumn2")
        
        png (file.path (Out_dir, paste0 ("ALL_DMR_heatmap_Par_vs_off", context, ".png")), width = 400, height = 400)
        pheatmap (Heatmap_matrix, 
                  main = context,
                  color = colorRampPalette (c ("white", "red"))(100),
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  legend = FALSE
        )
        dev.off()
        
        pdf (file.path (Out_dir, paste0 ("ALL_DMR_heatmap_Par_vs_off", context, ".pdf")))
        pheatmap (Heatmap_matrix, 
                  main = context,
                  color = colorRampPalette (c ("white", "red"))(100),
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  legend = FALSE
        )
        dev.off()
      }
      
    }
  }
  
  ## Ontology Plots
  {
    ### Load GFF file as GRanges object and filter
    {
      GFF_data_genes <- rtracklayer::import (Gene_annotation_file, format = "gff")
      GFF_data_mRNA <- GFF_data_genes [GFF_data_genes$type == "mRNA"]
      GFF_data_mRNA$score [is.na (GFF_data_mRNA$score)] <- 0 # for testing writing to bed file
      rtracklayer::export (GFF_data_mRNA, format = "BED", "/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Par_offspring_deeptools/allmRNA.bed")
      GFF_data_genes <- GFF_data_genes [GFF_data_genes$type == "mRNA"]
      
      ## Load TEs and remove partial matches and unknown TEs
      GFF_data_TEs <- rtracklayer::import (GFF_filepath_TEs, format = "gff")
      GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$type == "match"]
      GFF_data_TEs <- GFF_data_TEs [!is.na (GFF_data_TEs$Order)]
      GFF_data_TEs <- GFF_data_TEs [GFF_data_TEs$Order != "NA"]
      GFF_data_TEs <- GFF_data_TEs [seqnames (GFF_data_TEs) %in% seqnames (DMR_calling_chromosomes)]
      GFF_data_TEs <- unique (GFF_data_TEs)
      # Filter for TEs >500bp in length and where order 
      GFF_data_TEs_500bp <- GFF_data_TEs [width (GFF_data_TEs) > 500] # length (GFF_data_TEs) = 129220
      non_na_idx <- !(GFF_data_TEs_500bp$Superfamily == "NA")
      GFF_data_TEs_500bp$Order [non_na_idx] <- paste (GFF_data_TEs_500bp$Order [non_na_idx], GFF_data_TEs_500bp$Superfamily [non_na_idx], sep = "_")
      Cutoff_threshold <- 0.1 * length (GFF_data_TEs_500bp)
      Cutoff_TE_elements <- table (GFF_data_TEs_500bp$Order) < Cutoff_threshold
      mcols (GFF_data_TEs_500bp)$Order [mcols (GFF_data_TEs_500bp)$Order %in% names (Cutoff_TE_elements [Cutoff_TE_elements])] <- "Other_TE"
      mcols (GFF_data_TEs_500bp)$type <- mcols (GFF_data_TEs_500bp)$Order
      
      GFF_data_promoters <- GenomicRanges::flank (GFF_data_genes [GFF_data_genes$type == "mRNA"], width = promoter_width, start = TRUE)
      GFF_data_promoters$type <- "promoter"
      
      GFF_combined <- c (GFF_data_genes, GFF_data_TEs_500bp, GFF_data_promoters)
    }
    
    ### Overlap with Genes and Promotors
    {
      Pairwise_comparison_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakProgeny_vs_par/analysis/objects/Pairwise_Comparisons/"
      Contexts <- c ("CG", "CHG", "CHH")
      DMR_methylation <- c ("gain", "loss")
      Promoter_or_gene <- c ("promoter", "gene")
      min_sig_genes <- 10
      
      for (context in Contexts) {
        for (methylation in DMR_methylation) {
          for (element_of_interest in Promoter_or_gene) {
            print (context)
            print (methylation)
            print (element_of_interest)
            
            genes_of_interest <- paste0 ("/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/Progeny_", context, "_overlapping_", element_of_interest, "_", methylation, ".RData")
            All_context_DMRs <- qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_AmbientParents_vs_AmbientOffspring.qs")))
            All_context_DMRs <- All_context_DMRs [All_context_DMRs$regionType == methylation]
            if (element_of_interest == "gene") {
              GFF_file <- GFF_data_genes
            } else {
              GFF_file <- GFF_data_promoters
            }
            Overlaps_genes <- GFF_file [overlapsAny (GFF_file, All_context_DMRs)]
            save (Overlaps_genes, file = genes_of_interest)
            
            enrichment_RData <- paste0 (GO_out_folder, "/Enrichment_results_", basename (genes_of_interest))
            enrichment <- GO_enrichment (
              genes_of_interest = genes_of_interest,
              GO_database = GO_database, 
              gff_file = gff_file,
              out_folder = GO_out_folder,
              useInfo = "all",
              min_sig_genes = 10,
              fisher_threshold = 0.05)
            
            plot_title <- paste (context, element_of_interest, methylation)
            GO_plots(enrichment_RData = enrichment_RData,
                     show_ontology = TRUE,
                     min_sig_genes = 10,
                     max_y_lim = 1000, 
                     fisher_threshold = 0.02,
                     plot_title = plot_title,
                     plot_width = 43.5,
                     out_folder = "/rds/projects/l/lunadiee-epi-virtualmchine/Jack/GO_oak_seasons_generations/Test_dataset/output_Joe/Ontologies_oak_seasons/")
          }
        }
      }
    }
  }
}

# Supplimentary Figures (some are generated by the main figure sections)
{
  ## Venn Diagrams
  {
    ### Venn diagrams for progeny and seasons
    {
      #### Load data
      {
        ###### Random variables
        {
          library ("VennDiagram")
          
          Pairwise_comparison_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons"
          Groups <- c ("SpringElevatedParents", "SummerElevatedParents", "AutumnElevatedParents", "SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")
          
          context <- "CHH"
        }
        
        ###### Load DMRs from genetically identical parents across all 3 seasons
        {
          # Filepath where genetic controlled DMRs are stored
          Genetic_control_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakGeneticControl_eCO2/analysis/objects/Pairwise_Comparisons"
          Genetic_control_combinations <- c ("SpringAmbientParents_vs_SummerAmbientParents", "SpringAmbientParents_vs_AutumnAmbientParents", "SummerAmbientParents_vs_AutumnAmbientParents")
          All_genetic_control_context_DMRs <- lapply (Genetic_control_combinations, function (combination) {
            if (file.exists (file.path (Genetic_control_filepath, paste0 ("DMRs_", context, "_", combination, ".qs")))) {
              qread (file.path (Genetic_control_filepath, paste0 ("DMRs_", context, "_", combination, ".qs")))
            }
          })
          names (All_genetic_control_context_DMRs) <- Genetic_control_combinations
        }
        
        ###### Plot venn for seasons vs progeny
        {
          Contexts <- c ("CG", "CHG", "CHH")
          for (context in Contexts) {
            Pairwise_combinations <- combn (Groups, 2, paste, collapse = "_vs_")
            Reversed_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
            
            All_context_DMRs <- lapply (Reversed_combinations, function (DMR_pair) {
              if (file.exists (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
                qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
              }
            })
            names (All_context_DMRs) <- Reversed_combinations
            
            gr1 <- All_context_DMRs [["AmbientOffspring_vs_SpringAmbientParents"]]
            gr2 <- All_context_DMRs [["AmbientOffspring_vs_SummerAmbientParents"]]
            gr3 <- All_context_DMRs [["AmbientOffspring_vs_AutumnAmbientParents"]]
            
            # Create lists for Venn diagram
            list1 <- get_ranges(gr1)
            list2 <- get_ranges(gr2)
            list3 <- get_ranges(gr3)
            
            Venn_data <- list(
              "Set1" = unique(list1),
              "Set2" = unique(list2),
              "Set3" = unique(list3)
            )
            
            output_path <- file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all")
            venn_filename_base <- paste0 ("Venn_Seasons_Offspring_", context)
            
            venn_diagram <- venn.diagram (
              x = Venn_data,
              category.names = c ("Offspring\nvs\nSpring", "Offspring\nvs\nSummer", "Offspring\nvs\nAutumn"),
              fill = viridis (3, option = "turbo"),
              alpha = 0.2,
              print.mode = c ("raw", "percent"),
              main = context,
              main.cex = 1.5, 
              main.fontface = "bold", 
              filename = file.path (output_path, paste0 (venn_filename_base, ".png")),
              width = 15,
              height = 15,
              units = "cm",
              margin = 0.05,
              output = TRUE,
              cex = 1,
              cat.cex = 1,
              cat.dist = c (0.11, 0.11, 0.11)
            )
            
            pdf_file <- file.path(output_path, paste0(venn_filename_base, ".pdf"))
            venn_diagram <- venn.diagram(
              x = Venn_data,
              category.names = c ("Offspring\nvs\nSpring", "Offspring\nvs\nSummer", "Offspring\nvs\nAutumn"),
              fill = viridis (3, option = "turbo"),
              alpha = 0.2,
              print.mode = c ("raw", "percent"),
              main = context,
              filename = NULL,
              output = FALSE,
              width = 10,
              height = 10,
              units = "cm",
              margin = 0.1,
              output = TRUE,
              cex = 1.5,
              cat.cex = 1.5,
              cat.dist = c (0.13, 0.13, 0.13)
            )
            pdf (pdf_file)
            grid.draw(grobTree(venn_diagram))
            dev.off()
          }
          
        }
      }
      
      ####
      {
        
      }
    }
    
    ### Venn Diagrams for genetically identical individuals
    {
      Pairwise_comparison_filepath <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/Pairwise_Comparisons/"
      Pairwise_comparison_filepath_genetic <- "/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakGeneticControl_AmbientParents_eCO2/analysis/objects/Pairwise_Comparisons/"
      Groups <- c ("SpringElevatedParents", "SummerElevatedParents", "AutumnElevatedParents", "ElevatedOffspring", "SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")
      
      Contexts <- c ("CG", "CHG", "CHH")
      for (context in Contexts) {
        DMRs <- lapply (paste0 (Pairwise_comparison_filepath_genetic, "/DMRs_", unlist (as.list (seqnames (DMR_calling_chromosomes))), "_", context, "_Aa6387_vs_Aa8673.qs"), qread)
        All_DMRs_genetic <- do.call (c, DMRs)
        
        Pairwise_combinations <- combn (Groups, 2, paste, collapse = "_vs_")
        Reversed_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
        
        All_context_DMRs <- lapply (Reversed_combinations, function (DMR_pair) {
          if (file.exists (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
            qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
          }
        })
        names (All_context_DMRs) <- Reversed_combinations
        
        gr1 <- All_context_DMRs [["SummerAmbientParents_vs_SpringAmbientParents"]]
        gr2 <- All_context_DMRs [["AutumnAmbientParents_vs_SpringAmbientParents"]]
        gr3 <- All_DMRs_genetic
        
        list1 <- get_ranges(gr1)
        list2 <- get_ranges(gr2)
        list3 <- get_ranges(gr3)
        
        Venn_data <- list(
          "Set1" = unique(list1),
          "Set2" = unique(list2),
          "Set3" = unique(list3)
        )
        
        output_path <- file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all")
        venn_filename_base <- paste0 ("Venn_SpringSummer_SpringAutumn_Individual12_", context)
        
        venn_diagram <- venn.diagram (
          x = Venn_data,
          category.names = c ("Spring\nvs\nSummer", "Sping\nvs\nAutumn", "Individual1\nvs\nIndividual2"),
          fill = viridis (3, option = "turbo"),
          alpha = 0.2,
          print.mode = c ("raw", "percent"),
          main = context,
          main.cex = 1.5, 
          main.fontface = "bold", 
          filename = file.path (output_path, paste0 (venn_filename_base, ".png")),
          width = 15,
          height = 15,
          units = "cm",
          margin = 0.05,
          output = TRUE,
          cex = 1,
          cat.cex = 1,
          cat.dist = c (0.11, 0.11, 0.11)
        )
        
        pdf_file <- file.path(output_path, paste0(venn_filename_base, ".pdf"))
        venn_diagram <- venn.diagram(
          x = Venn_data,
          category.names = c ("Spring\nvs\nSummer", "Sping\nvs\nAutumn", "Individual1\nvs\nIndividual2"),
          fill = viridis (3, option = "turbo"),
          alpha = 0.2,
          print.mode = c ("raw", "percent"),
          main = context,
          filename = NULL,
          output = FALSE,
          width = 10,
          height = 10,
          units = "cm",
          margin = 0.1,
          output = TRUE,
          cex = 1.5,
          cat.cex = 1.5,
          cat.dist = c (0.13, 0.13, 0.13)
        )
        pdf (pdf_file)
        grid.draw(grobTree(venn_diagram))
        dev.off()
      }
    }
  }
  
  ## Heatmaps
  {
    
  }
  
  ## Barplots
  {
    
  }
  
  ## PCAs
  {
    ### PCA of TEs by class of TE
    {
      
    }
  }
  
  ## Ontologies
  {
    
  }
  
  ## Histograms/Violin plots
  {
    
  }
  
  ## Deeptools
  {
    
  }
  
  ## Overlapping DMRs
  {
    context <- "CG"
    
    All_context_DMRs <- lapply (Reversed_combinations, function (DMR_pair) {
      if (file.exists (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))) {
        qread (file.path (Pairwise_comparison_filepath, paste0 ("DMRs_", context, "_", DMR_pair, ".qs")))
      } else if (context == "Random") {
        function_generate_random_regions (sample_number, Chromosome_lengths, region_width)
      }
    })
    names (All_context_DMRs) <- Reversed_combinations
    
    Spring_summer_DMRs <- All_context_DMRs [["SummerAmbientParents_vs_SpringAmbientParents"]]
    Spring_autumn_DMRs <- All_context_DMRs [["AutumnAmbientParents_vs_SpringAmbientParents"]]
    Summer_autumn_DMRs <- All_context_DMRs [["AutumnAmbientParents_vs_SummerAmbientParents"]]
    
    Par_vs_Off <- "AmbientParents_vs_AmbientOffspring"
    Par_vs_Off_DMRsCG <- qread (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakProgeny_vs_par/analysis/objects/Pairwise_Comparisons",
                                           paste0 ("DMRs_CG", "_", Par_vs_Off, ".qs")))
    
    a <- Spring_summer_DMRs [overlapsAny (Spring_summer_DMRs, Spring_autumn_DMRs)]
    a <- a [overlapsAny (a, Summer_autumn_DMRs)]
    a <- a [overlapsAny (a, Par_vs_Off_DMRsCG)]
  }
  
  ## Save DMRs as csv files for supplimentary tables
  {
    ### Load data of DMRs
    {
      Pairwise_combinations <- combn (names (All_analysis_groups), 2, paste, collapse = "_vs_")
      Pairwise_combinations <- sapply (strsplit (Pairwise_combinations, "_vs_"), function(x) paste (rev (x), collapse = "_vs_"))
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
    }
    
    All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
    
    Relevant_DMR_comparisons <- c ("SummerAmbientParents_vs_SpringAmbientParents", 
                                   "AutumnAmbientParents_vs_SpringAmbientParents", 
                                   "AutumnAmbientParents_vs_SummerAmbientParents",
                                   "AmbientOffspring_vs_SpringAmbientParents", 
                                   "AmbientOffspring_vs_SummerAmbientParents", 
                                   "AmbientOffspring_vs_AutumnAmbientParents")
    
    for (context in Contexts) {
      for (comparison in Relevant_DMR_comparisons) {
        Out_csv_data <- All_DMRs [[context]][[comparison]]
        print (paste (context, length (Out_csv_data), "DMRs"))
        Out_filepath <- file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/PDF_figures/Supp_figures/ST/DMRs/",
                                   paste0 (context, "_", comparison, "_", "DMRs.csv"))
        write.csv (Out_csv_data, file = Out_filepath, row.names = FALSE)
      }
    }
  }
}
