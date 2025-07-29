## This script investigates the hypothesis that there is an overrepresentation of CG methylation
## at loci that are hypomethyalted in the CHH context to compensate

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
        ggsave(filename = paste0(out_folder, category, "_GO_graph_", out_file_ID, show_ontology, ".png"), 
               plot = GO_graphs[[category]], 
               width = plot_width, 
               height = 1 * nrow (filtered_data) + 2, units = "cm")
        ggsave(filename = paste0(out_folder, category, "_GO_graph_", out_file_ID, show_ontology, ".pdf"), 
               plot = GO_graphs[[category]], 
               width = plot_width, 
               height = 1 * nrow (filtered_data) + 2, units = "cm")
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

# For the heatmap of CG contexts, find the most overmethylated (comparing Spring and other seasons)
# and make heatmap of those loci in the CG context
{
  ## Load all heatmap data
  {
    All_individuals <- c (SpringAmbientParents, SummerAmbientParents, AutumnAmbientParents, AmbientOffspring)
    All_heatmap_data <- lapply (All_individuals, function (individual) {
      qread (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps", paste0 ("Heatmap_recomputed_ambient_par_and_ambient_off_", individual, ".qs")))
    })
    names (All_heatmap_data) <- All_individuals
  }
  
  ## Get DMRs from Spring vs Autumn CG where DMRs absolute difference > 0.4
  {
    context <- "CG"
    Spring_Autumn_0.4_diff_DMRs <- All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_0.4_diff_DMRs <- All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_Autumn_0.4_diff_DMRs <- unique (c (Spring_Autumn_0.4_diff_DMRs, Spring_Summer_0.4_diff_DMRs))
  }
  
  ## Compute heatmap
  {
    context <- "CHH"
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
  }
  
  ## Subset from heatmap for the Spring_Summer_Autumn_0.4_diff_DMRs and plot
  {
    Matrix_names <- paste (seqnames (Spring_Summer_Autumn_0.4_diff_DMRs), 
                           start (Spring_Summer_Autumn_0.4_diff_DMRs),
                           end (Spring_Summer_Autumn_0.4_diff_DMRs),
                           sep = "_")
    
    Matrix_CG_0.4 <- Heatmap_matrix [Matrix_names,]
    # Matrix_CHH_0.4 <- Heatmap_matrix_CHH [Matrix_names,]
    
    # Heatmap_matrix <- Matrix_CHH_0.4
    if (max (nrow (Heatmap_matrix)) > 10000) {
      Heatmap_matrix <- Heatmap_matrix [sample (nrow (Heatmap_matrix), 10000, replace = FALSE), ]
    }
    
    colnames (Heatmap_matrix) <- gsub ("Ambient", "", colnames (Heatmap_matrix))
    colnames (Heatmap_matrix) <- gsub ("Parents", "", colnames (Heatmap_matrix))
    
    colnames (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3",
                                    "Spring6","Summer1", "Summer2", "Autumn1", 
                                    "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
    
    png (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                    paste0 ("DMR_heatmap_recomputed_ambient_par_and_ambient_offspring_0.4_", context, ".png")), width = 400, height = 400)
    a <- pheatmap (as.matrix (Heatmap_matrix), 
                   main = context,
                   color = colorRampPalette (c ("white", "red"))(100),
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   fontsize_col = 14,
                   legend = TRUE
    )
    dev.off()
    pdf (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                    paste0 ("DMR_heatmap_recomputed_ambient_par_and_ambient_offspring_0.4_", context, ".pdf")))
    print (pheatmap (as.matrix (Heatmap_matrix), 
                     main = context,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     fontsize_col = 14,
                     legend = TRUE
    ))
    dev.off()
  }
}

# Compute methylation of CG regions in demethylated CHH contexts
{
  ## Get DMRs from Spring vs Autumn and Spring vs Summer CHH where DMRs absolute difference > 0.4
  {
    context <- "CHH"
    Spring_Autumn_0.4_diff_DMRs <- All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_0.4_diff_DMRs <- All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_Autumn_0.4_diff_DMRs <- unique (c (Spring_Autumn_0.4_diff_DMRs, Spring_Summer_0.4_diff_DMRs))
  }
  
  ## Compute methylation in CG context of CHH DMRs > 0.4 
  {
    All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
    for (individual in All_individuals) {
      CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports/", paste0 ("CX_report_", individual, ".qs" )))
      Heatmap_data <- lapply (Contexts, function (context) {
        print (paste ("working on", context))
        DMRcaller::analyseReadsInsideRegionsForCondition (Spring_Summer_Autumn_0.4_diff_DMRs, CX_report, context = context)
      })
      names (Heatmap_data) <- Contexts
      dir.create (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis"))
      print (paste ("saving Heatmap_data to", file.path (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis", paste0 ("Heatmap_", individual, ".qs")))))
      qsave (Heatmap_data, file.path (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis", paste0 ("Heatmap_recomputed_Spring_summer_spring_autumn_CG_replace_CHH_hypothesis_", individual, ".qs"))))
    }
  }
  
  ## Plot heatmap of the CG methylation from CHH DMRs where the methylation difference is > 0.4
  {
    context <- "CG"
    All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
    All_heatmap_data <- lapply (All_individuals, function (individual) {
      qread (file.path (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis", paste0 ("Heatmap_recomputed_Spring_summer_spring_autumn_CG_replace_CHH_hypothesis_", individual, ".qs"))))
    })
    names (All_heatmap_data) <- All_individuals
    
    
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
    
    ### Make plot
    {
      png (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("CG_Overrepresentation_heatmap_recomputed_ambient_par_and_ambient_offspring.png")), width = 400, height = 400)
      a <- pheatmap (as.matrix (Heatmap_matrix), 
                     main = context,
                     color = colorRampPalette (c ("white", "red"))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     fontsize_col = 14,
                     legend = TRUE
      )
      dev.off()
      pdf (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("CG_Overrepresentation_heatmap_recomputed_ambient_par_and_ambient_offspring.pdf")))
      print (pheatmap (as.matrix (Heatmap_matrix), 
                       main = context,
                       color = colorRampPalette (c ("white", "red"))(100),
                       show_rownames = FALSE,
                       show_colnames = TRUE,
                       fontsize_col = 14,
                       legend = TRUE
      ))
      dev.off()
    }
    
    # Plot CG methylation of 386 most differentially methylated DMRs in the CHH context after clustering on all CG DMRs 
    {
      Most_diff_DMRs <- rownames (Heatmap_matrix)
      
      All_individuals <- c (SpringAmbientParents, SummerAmbientParents, AutumnAmbientParents, AmbientOffspring)
      All_heatmap_data <- lapply (All_individuals, function (individual) {
        qread (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps", paste0 ("Heatmap_recomputed_ambient_par_and_ambient_off_", individual, ".qs")))
      })
      names (All_heatmap_data) <- All_individuals
      
      context <- "CHH"
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
      
      colnames (Heatmap_matrix) <- gsub ("Ambient", "", colnames (Heatmap_matrix))
      colnames (Heatmap_matrix) <- gsub ("Parents", "", colnames (Heatmap_matrix))
      
      colnames (Heatmap_matrix) <- c ("Spring5", "Spring4", "Spring1", "Spring2", "Spring3",
                                      "Spring6","Summer1", "Summer2", "Autumn1", 
                                      "Autumn2", "Autumn3", "Progeny1A", "Progeny1B", "Progeny2")
      
      annotation_row <- data.frame ("CHH_Diff" = factor (rownames (Heatmap_matrix) %in% Most_diff_DMRs,
                                                         levels = c (TRUE, FALSE)))
      rownames (annotation_row) <- rownames (Heatmap_matrix)
      
      png (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("Overrepresented_CG_heatmap_recomputed_ambient_par_and_ambient_offspring_", context, ".png")), width = 400, height = 400)
      a <- pheatmap (as.matrix (Heatmap_matrix), 
                     main = context,
                     color = colorRampPalette (c ("white", "red"))(100),
                     annotation_row = annotation_row,
                     annotation_colors = list (CHH_Diff = c("TRUE" = "black", "FALSE" = "white")),
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     fontsize_col = 14,
                     legend = TRUE
      )
      dev.off()
      
      pdf (file.path ("/rds/projects/l/lunadiee-epi-virtualmchine/Epigenomic_landscape_of_Oak/Scripts/Plot_all", 
                      paste0 ("Overrepresented_CG_heatmap_recomputed_ambient_par_and_ambient_offspring_", context, ".pdf")))
      print (pheatmap (as.matrix (Heatmap_matrix), 
                       main = context,
                       color = colorRampPalette (c ("white", "red"))(100),
                       annotation_row = annotation_row,
                       annotation_colors = list (CHH_Diff = c("TRUE" = "black", "FALSE" = "white")),
                       show_rownames = FALSE,
                       show_colnames = TRUE,
                       fontsize_col = 14,
                       legend = TRUE
      ))
      dev.off()
    }
    
    ## Compute Methylation in CG context of CHH DMRs with low stdev (WIP)
    {
      Spring_Autumn_lowstdev_diff_DMRs <- All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents [
        abs (mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
               mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
      
      Spring_Summer_lowstdev_diff_DMRs <- All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents [
        abs (mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
               mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    }
  }
}  

# Compute the methyaltion of CHH in the 2855 DMRs from the most overmethylated CG DMRs
{
  ## Get DMRs from Spring vs Autumn CG where DMRs absolute difference > 0.4
  {
    context <- "CG"
    Spring_Autumn_0.4_diff_DMRs <- All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$AutumnAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_0.4_diff_DMRs <- All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents [
      abs (mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion1"] - 
             mcols (All_DMRs [[context]]$SummerAmbientParents_vs_SpringAmbientParents) [, "proportion2"]) > 0.4]
    
    Spring_Summer_Autumn_0.4_diff_DMRs <- unique (c (Spring_Autumn_0.4_diff_DMRs, Spring_Summer_0.4_diff_DMRs))
  }
  
  ## Compute methylation in CHH context of CG DMRs > 0.4 
  {
    All_individuals <- Reduce (c, All_analysis_groups [c ("SpringAmbientParents", "SummerAmbientParents", "AutumnAmbientParents", "AmbientOffspring")])
    for (individual in All_individuals) {
      CX_report <- qread (file.path (Working_dir, Specondition, "analysis/objects/Individual_CX_reports/", paste0 ("CX_report_", individual, ".qs" )))
      Contexts <- "CHH"
      Heatmap_data <- lapply (Contexts, function (context) {
        print (paste ("working on", context))
        DMRcaller::analyseReadsInsideRegionsForCondition (Spring_Summer_Autumn_0.4_diff_DMRs, CX_report, context = context)
      })
      names (Heatmap_data) <- Contexts
      dir.create (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis"))
      print (paste ("saving Heatmap_data to", file.path (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis", paste0 ("Heatmap_", individual, ".qs")))))
      qsave (Heatmap_data, file.path (file.path (Working_dir, Specondition, "analysis/objects/Heatmaps/CG_replace_CHH_hypothesis", paste0 ("Heatmap_recomputed_Spring_summer_spring_autumn_CG_replace_CHH_hypothesis_CHH_methylation_in_CG_DMRs", individual, ".qs"))))
    }
  }
}

