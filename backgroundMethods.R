
medianWithoutNA<-function(x) {
   median(x[which(!is.na(x))])
}

getBackgroundList <- function(backgroundCDS){
    # fill a list where each named entry gives a vector with the gene sums
    #    across all empty droplets in that sample
    backgroundCounts = list()
    sample.ids = levels(as.factor(colData(backgroundCDS)$sample))

    # get the data for each subset
    for (sample.id in sample.ids){
        backgroundCounts[[sample.id]] = Matrix::rowSums(
            exprs(backgroundCDS[,colData(backgroundCDS)$sample == sample.id])
            ) # Jonathan has a step here to get only UMI < threshold, but I did that in the inital step
    }

    return(backgroundCounts)
}

# Requires: inputCDS is a cds holding cells whose expression values are to be corrected using the Packer method
# Modifies: nothing
# Effects: Calculates and returns a dense matrix. Rows are genes, columns are cells.
#          Counts are scaled by size factor so that the total UMIs per cell are equal. 
#          Counts are log2 transformed, and any non-expressed genes are dropped.
#          The matrix is then returned
getScaledExprs <- function(inputCDS){
    # Get exprs after scaling by size factor
    sfScaledExprs =  Matrix::t(Matrix::t(exprs(inputCDS)) / colData(inputCDS)$Size_Factor) # Double transpose to make matrix syntax work
    # Jonthan put this into log scale before finding the PCA. I'll try this as well. I'm presuming that the logic
    #   is to not let the PC choice/magnitude note get dominated by hyper expressed genes. Probably
    #   helpful here, since things like Kap are a screaming hot signal in proximal tube, etc.
    sfScaledExprs@x = log2(sfScaledExprs@x + 1) # pseudocount to not explode empty counts
    # Check for any overflow issues and remove non-expressed genes
    dim(sfScaledExprs)
    sfScaled_rowsums = Matrix::rowSums(sfScaledExprs)
    expressedRowsOnly = sfScaledExprs[is.finite(sfScaled_rowsums) & sfScaled_rowsums != 0,]
    dim(expressedRowsOnly)

    return(expressedRowsOnly)
}

# Requires: expressedRowsOnly is a dense matrix of genes X cells (rows x columns), should
#      be the output of getScaledExprs. PCs_to_use is an integer of how many principle compoonents to consider
getCellPCA <- function(expressedRowsOnly, PCs_to_use=50){
    # Get the PCs I want to keep
    useScaling = TRUE
    irlba_res <- monocle3:::sparse_prcomp_irlba(t(expressedRowsOnly), n=min(PCs_to_use, min(dim(expressedRowsOnly)) - 1),
                                    center=useScaling, scale=useScaling)

    # Return
    return(irlba_res)
}

# Requires: preprocessedExprs is the ouptut of getScaledExprs, a log2 transformed count matrix from the cellCDS
#         where non-expressed genes have been dropped and counts have been scaled by each cell's size factor
#         input_irlba_res is a matrix of cells in PC space
#         backgroundLIstInput is a named list, where each sample names a CDS with only that sample's background empty droplets
#         inputCDS is the CDS holding cell data
 getBG_loadings <- function(preprocessedExprs, input_irlba_res, backgroundListInput, inputCDS){
    print("Preprocessing background data")
    preprocessedExprs.t = t(preprocessedExprs)
    preprocessedExprs.t.center = Matrix::colMeans(preprocessedExprs.t)
    # Get the standard deviation of each gene so they can be standardized shortly
    # NOTE: This is modified from Jonathan's original code. Using this framing for fast computation on large matrices
    preprocessedExprs.t.varScaling = (Matrix::colMeans(preprocessedExprs.t^2) - # E[x^2] 
                             ((Matrix::colMeans(preprocessedExprs.t))**2)   # E[x] ^ 2
                                 )
    preprocessedExprs.t.varScaling = unname(sqrt(preprocessedExprs.t.varScaling))

    print("Getting background count totals")
    norm.background.counts = list()
    for (sample.id in names(backgroundListInput)){
        # Note: Jonathan had a line here to subset to the genes used for ordering. In his case, that was a subset of
        #    genes that were over-dispersed
        tmp = backgroundListInput[[sample.id]][names(preprocessedExprs.t.center)]
        tmp = tmp / (sum(backgroundListInput[[sample.id]]) / (colData(inputCDS)$n.umi / colData(inputCDS)$Size_Factor)[1] )
        norm.background.counts[[sample.id]] = (tmp - preprocessedExprs.t.center) / preprocessedExprs.t.varScaling
    }

    print("Transforming background to PC space")
    # Get the background in pca space
    backgroundPCA_list = list()
    for (sample.id in names(backgroundListInput)){
        backgroundPCA_list[[sample.id]] = as.vector(t(norm.background.counts[[sample.id]] %*% 
                                                        input_irlba_res$rotation))
    }
    #browser() 

    # Helper function
    magFunc = function(v) sqrt(sum(v^2))
    bgProj = list()
    # Get projection magnitude for each cell onto each bg in PCA space
    topDimPCA <- input_irlba_res$x
    for (sample.id in names(backgroundPCA_list)) {
        bgProj[[sample.id]] = ((topDimPCA %*% as.matrix(backgroundPCA_list[[sample.id]])) /
                                magFunc(backgroundPCA_list[[sample.id]]) )
    }

    return(bgProj)
 }

getBG_normed_PCA <- function(inputCDS, inputBG_proj, cellsInPCA_space){
    # Add this data to the main cell CDS for book-keeping/to help setting up the model matrix
    for (sample.id in names(inputBG_proj)){
        thisProjColName = paste0("bgProj_", sample.id)
        colData(inputCDS)[,thisProjColName] = as.vector(inputBG_proj[[sample.id]])
    }

    # Now set up a regression model to see the extent to which PCA scores for cells can be explained en masse
    #    by those cells' projection magnitude upon the background vectors
    model_matrix_background = model.matrix(
        formula(paste("~", paste("bgProj_", names(inputBG_proj), sep="", collapse=" + "))),
        data=colData(inputCDS), drop.unused.levels=T )

    # Fit it
    backgroundFit = limma::lmFit(t(cellsInPCA_space), model_matrix_background)
    # Pull the coefficients
    backgroundBeta = backgroundFit$coefficients[, -1, drop=F]
    backgroundBeta[is.na(backgroundBeta)] = 0

    # Regress out the background effects
    bgNormPCA = t(t(cellsInPCA_space) - backgroundBeta %*% t(model_matrix_background[,-1]))

    return(bgNormPCA)
}

getUMAP <- function(inputPCA, yToUse=NULL, targWeight=.5, aVal=1.58, bVal=.9){
    # UMAP Plotting
    umapComponents = 2
    verbose = TRUE

    extra_arguments = list(
        n_neighbors = 15,
        # min_dist = 0.1,
        a=aVal,
        b=bVal,
        metric = "cosine")

    umap_args_here <- c(list(
            X = inputPCA,
            y=yToUse,
            target_weight=targWeight,
          #  log = F,
            n_component = as.integer(umapComponents), 
            verbose = verbose #,  return_all = T
            ),
                   
            extra_arguments[names(extra_arguments) %in% 
                c("python_home", "n_neighbors", "metric", "n_epochs", 
                  "negative_sample_rate", "alpha", "init", "min_dist", 
                  "spread", "set_op_mix_ratio", "local_connectivity", 
                  "bandwidth", "gamma",
                   "a", "b",
                   "random_state", 
                  "metric_kwds", "angular_rp_forest", "verbose")])

    print("umap_args_here:")
    # print(umap_args_here)
    print(str(umap_args_here))

    set.seed(17)
    system.time({
        umapRes <- do.call(uwot::umap, umap_args_here)
        })

    return(umapRes)
}

# Requires: inputCellCDS is a CDS of putative cells. Assumes pre-filtered to all be
#    above a basic UMI threshold, etc.
#
# inputBackgroundCDS is a CDS of empty droplets. This is assumed to have been filtered
#     before passing, e.g. already only be "cells" with <= 15 UMIs.
#     Also, coldata(inputBackgroundCDS) MUST CONTAIN A COLUMN "sample", a character vector, marking
#     the barcode-combo-of-origin for each given empty droplet
# Modifies: the inputCellCDS, adding columns for UMAP coordiantes 1/2 (after finding a 2d embedding)
#     using both no normalization (umap_X_no_norm) or after regressing out background in PC space using Jonathan's
#     method (umapX_bg_norm)
# Effects: Calculates background per Jonathan's method, returns a named list including
#      1) CDS with UMAP coordinates using the correction method
#      2) The cells' coordinates in PCA space (topDimPCA)
#      3) The backgrounds as "cells" in PCA space (bgNormPCA)
#      4) The magnitude of cells' projections onto the various backgrounds (backgroundProjections)
#      2-4 are for convenience, 1 is probably mostly what's useful
findBackground_Packer_Method <- function(inputCellCDS, inputBackgroundCDS, 
        PCs_to_use=50, # Can specifiy number of principle components to find, otherwise defaults to 50
        returnOnlyMatrix=TRUE # If false, return a list including a modified inputCellCDS, adjusted matrix of cells in pca space,
                            # background vectors in PC space, and projections of all cells onto background vectors
        ){
    set.seed(17)

    backgroundList <- getBackgroundList(inputBackgroundCDS)

    print("Getting PCA for cells")
    # Get the PCA info
    scaledExprs = getScaledExprs(inputCellCDS)
    irlba_res_inputCellCDS <- getCellPCA(scaledExprs, PCs_to_use)
    topDimPCA <- irlba_res_inputCellCDS$x # This is the data of cells in PC space, PCs x cells (rows x columns)

    print("Calculating background")
    # Get background counts
    backgroundProjections <- getBG_loadings(scaledExprs, irlba_res_inputCellCDS, 
                    backgroundList, inputCellCDS)
    # Fit a model to account for PCs of each cell using the projections as features. Get an updated
    #    matrix after regressing out the fit accounted for by this model
    print("Regressing out projections onto background(s)")
    bgNormPCA = getBG_normed_PCA(inputCellCDS, backgroundProjections, topDimPCA)

    if (returnOnlyMatrix){
            return(bgNormPCA)
        } else{
             return( list("CDS"=inputCellCDS, "topDimPCA" = topDimPCA,
         "bgNormPCA"=bgNormPCA, "backgroundProjections" = backgroundProjections ))
        }
}


# Requires: filtered CDS is a cds object
#           correctedCounts is a dgCMatrix produced by SoupX with updated UMI values
#           after background subtraction. Note that the default is a different type of matrix
#           which needs to be manually cast to a dgCMatrix
updateWithNewCounts <- function(originalCDS, correctedCounts){
    # Make the new CDS, updating the counts matrix (don't decompress!)
    updatedCDS = new_cell_data_set(correctedCounts,
         cell_metadata = as.data.frame(colData(originalCDS)),
         gene_metadata = as.data.frame(rowData(originalCDS)) )      
    return (updatedCDS)
}

# Requires:
#    inputCDS and inputBackgroundCDS are monocle CDS objects. inputBackgroundCDS should contain
#        "cells" below a small UMI cutoff such that all "cells" are presumably empty droplets
#    inputSoupXassumption is an number 0-100, which tells SoupX what is the assumed percentage
#         of UMIs that are - on average - from background
#    PCs_to_use is a number of principal components to consider. Must be <<< number of expressed genes.
#       Typically only use something like 50-100, and due to the sparse computation it should be
#        far fewer than number of genes input
# Modifies:
#     inputCDS
# Effects:
#       returns a new CDS object with the same number of rows/columns as before, but now with
#        the exprs matrix replaced with new values where background has been subtracted by 
#      soupX according to the input assumption and provided background data
correctCDSbySoupX_all_at_once <- function(inputCDS, inputBackgroundCDS, inputSoupXassumption, PCs_to_use){
    set.seed(17)
    # First, make these into a SoupChannel object
    cellMetaDataDF = data.frame("nUMIs" = as.numeric(colData(inputCDS)$n.umi))
    row.names(cellMetaDataDF) <- rownames(colData(inputCDS))
    scObj <- SoupChannel(tod=exprs(inputBackgroundCDS),
                        toc=exprs(inputCDS),
                        metaData=cellMetaDataDF)

    # Set the contamination fraction
    scObj = setContaminationFraction(scObj, (inputSoupXassumption * .01))
    print("scObj set up for SoupX Correction")

    # Pre-process, cluster, and use these clusters to improve the SoupX correction
    inputCDS <- preprocess_cds(inputCDS, method="PCA", num_dim=PCs_to_use)
    inputCDS <- reduce_dimension(inputCDS)
    inputCDS <- cluster_cells(inputCDS, k=7)
    colData(inputCDS)$cluster_label <- as.character(clusters(inputCDS))
    print("Cells clustered for use in SoupX calculations")
    
    # Try to set clusters, if available:
    clusterAddition <- tryCatch({
        scObj = setClusters(scObj, colData(inputCDS)$cluster_label)
        clusterAddition <- TRUE
        }, error=function(err){
            print("Error loading partitions for SoupX. Running 
                only on a cell-by-cell basis")
            clusterAddition <- FALSE
            })

    # Now correct the counts
    scOutput = adjustCounts(scObj)
    print("Counts adjusted")
    # Convert to the right kind of sparse Matrix to put back into a CDS
    formattedSCoutput = as(scOutput, "dgCMatrix")

    # Generate a new CDS
    inputCDS <- updateWithNewCounts(inputCDS, formattedSCoutput)
    # Return it
    return(inputCDS)
}

# An organizer function to stash commands to generate plots covering stats
#   like UMI distributions per sample, mitochondrial RNA content per sample,
#   or doublet stats by sample
getQC_information <- function(inputCDS, processingNote){
    # First, get a simple barplot showing UMI distributions
    # # Barplot of counts
    ggplot(as.data.frame(colData(inputCDS)), aes(x=sample, y=log_10_umi)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0("./plots/", processingNote, "QC_Output_UMIs_by_Sample.png"))

    # See Mitochondrial Read distributions
    ggplot(as.data.frame(colData(inputCDS)), aes(x=sample, y=perc_mitochondrial_umis)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0("./plots/", processingNote, "QC_Output_Mito_Perc_By_Sample.png"))

    # See doublet numbers (Scrublet Estimates)
    ggplot(as.data.frame(colData(inputCDS)), aes(x=sample, y=scrublet_score)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0("./plots/", processingNote, "QC_Output_Scrublet_Scores_By_Sample.png"))

    # Get counts of number of cells:
    myCountData = (as.data.frame(colData(inputCDS))  %>% 
            group_by(sample) %>% summarise(n=n()))
    print(myCountData, n=21)
}

# Same basic method as correctCDSbySoupX_all_at_once, but does the whole process
#   separately for cells and background combos for each sample, as listed in the
#    "sample" column of colData(inputCDS)
correctCDSbySoupX_individually <- function(inputCDS, 
                inputBackgroundCDS, inputSoupXassumption, PCs_to_use){
    # Loop over a list of all unique $sample entries in the background CDS. 
    sampleVec = as.character(levels(as.factor(colData(inputCDS)$sample)))
    correctedList = list()
    print("Starting To correct with SoupX INdividually")
    for (eachSample in sampleVec){
        # Get the subset
        inputSubset = inputCDS[,colData(inputCDS)$sample == eachSample]
        bgSubset = inputBackgroundCDS[,colData(inputBackgroundCDS)$sample == eachSample]
        # Run background correction on this subset
        print(paste0("Correcting ", eachSample))
        correctedSubset = correctCDSbySoupX_all_at_once(inputSubset, bgSubset,
                                inputSoupXassumption, PCs_to_use)
        correctedList = c(correctedList, eachSample=correctedSubset)
    }

    print("Merging individually corrected CDS list")
    inputCDS <- monocle3::combine_cds(correctedList)

    return(inputCDS)
}



logOfN_fact_stirling <- function(inputN){
    # Give back ln(N!), found using the sterling approx
    return( log(sqrt(2*inputN*pi)) + 
                inputN*(log(inputN / exp(1))) )
}



findBG_log_likelihood_en_masse <- function(inputCDS, backgroundCDS){
    print("Precomputing factorials for possible gene counts")
    # First, pre-compute the factorials for all observed gene counts
    cellUMIs = Matrix::colSums(exprs(inputCDS))
    maxGeneCellCount = max(exprs(inputCDS))
    umiFactLookupVec = rep(1, (1+maxGeneCellCount))
    #names(umiFactLookupVec) = seq(from=0, to=maxCellUMI)
    # Hard set logs of 0, 1, 2
    umiFactLookupVec[1] = log(1)
    umiFactLookupVec[2] = log(1)
    umiFactLookupVec[3] = log(2) # Starts at 1...
    # Loops over values [3,maxCellUMI], finds sterling approx
    for (eachInd in 4:(maxGeneCellCount+1)){
        # Watch offset. index "4" = the log(3!) entry
        umiFactLookupVec[eachInd] = logOfN_fact_stirling(eachInd - 1) 
    }

    print("Getting background frequencies")
    # Second, get the background frequencies in the soup. Also note zero-background and remove
    #    (would explode the ln(0) part of the likelihood)
    fData(backgroundCDS)$n.umi <- Matrix::rowSums(exprs(backgroundCDS))
    bgProportions = fData(backgroundCDS)$n.umi / (sum(fData(backgroundCDS)$n.umi))
    # Need to only calculate over values where the background has nonzero counts (or kills chiSq)
    indicesToUse = which(bgProportions != 0)
    # Only use these background values
    bgProportions = bgProportions[indicesToUse]
    logBgProp = log(bgProportions) # Will use this value repeatedly

    # Generate a holder vec for the log likelihoods
    logLikeVec = rep(1, ncol(inputCDS))

    print("Finding log liklihoods for all cells")
    # Loop through and calculate
    for (eachCellInd in 1:ncol(inputCDS)){
        log_of_NfactTerm = logOfN_fact_stirling(cellUMIs[eachCellInd])

        # Now, get the counts. Note removal of no-background genes
        cellCounts = as.vector(exprs(cellCDS)[indicesToUse,eachCellInd])
        termOneSum = sum(cellCounts * logBgProp)
        termTwoSum = -1 * sum(umiFactLookupVec[cellCounts + 1])
        #termTwoSum = -1 * sum(sapply(cellCounts, logOfN_fact_stirling))
        
        # Add together for the final log likelihood
        logLikeVec[eachCellInd] = log_of_NfactTerm + termOneSum + termTwoSum

         if (eachCellInd%%100 == 0){
            cat(eachCellInd, "\n")
        }
    }

    # Add to the cds as a new column in colData, and return
    colData(inputCDS)$bgLogLike = logLikeVec
    return(inputCDS)
}





