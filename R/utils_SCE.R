
#' Extract information from SingleCellExperiment object to create CNA object
#' 
#' @param sce A SingleCellExperiment object. 
#' @param sample_key A character string denoting the name of the sample-level 
#' identifier (e.g. DonorID). 
#' @param sample_vars A character string or vector specifying which columns of 
#' `colData(sce)` should be included. This should include any columns for batch
#' and covariates. 
#' @param graph A KNN graph in sparseMatrix format, or the name of such a graph
#' stored in `colPair(sce)`. 
#' 
#' @returns An `rnca` data object (list)
#' 
#' @export 
createObject.SingleCellExperiment <- function(sce, 
                                              sample_key, 
                                              sample_vars, 
                                              graph) {    
    if(any(c(missing(sce), missing(sample_key), missing(sample_vars), missing(graph)))){
        stop("Parameters missing with no default.")
    }
    if(!isa(graph, "Matrix") && !(is.character(graph) && length(graph) == 1)){
        stop("'graph' must be a sparseMatrix or character string")
    }
    if(isa(graph, "Matrix")){
        if(!all(dim(graph) == ncol(sce))){
            stop("'graph' must be a square sparseMatrix with dims equal to 'ncol(sce)'")
        }
    } else if(is.character(graph)){
        stopifnot(graph %in% SingleCellExperiment::colPairNames(sce))
        graph <- SingleCellExperiment::colPair(sce, graph, asSparse = TRUE)
    }
    sample_vars <- c(sample_key, sample_vars)
    samplem_df <- SummarizedExperiment::colData(sce) |> 
        dplyr::select(dplyr::all_of(sample_vars)) |> 
        dplyr::distinct() |> 
        as.data.frame()
    n <- length(unique(samplem_df[, sample_key, drop = TRUE]))
    if(nrow(samplem_df) != n){
        stop("'sample_vars' must collapse to a unique value for each value in 'sample_key'",
             "Please check that sample_vars are sample-level variables.")
    }
    rownames(samplem_df) <- NULL
    
    obs_df <- SummarizedExperiment::colData(sce)[, sample_key, drop = FALSE] |> 
        as.data.frame()
    obs_df$CellID <- rownames(obs_df)
    rownames(obs_df) <- NULL
    
    rcna_object <- list(
        samplem = samplem_df,
        obs = obs_df, 
        connectivities = graph,
        samplem_key = sample_key,
        obs_key = 'CellID',
        N = nrow(samplem_df)
    )
    return(rcna_object)
}

#' Perform CNA on and assign the results to the input SingleCellExperiment object. 
#' 
#' @param sce A SingleCellExperiment object. 
#' @param test_var Variable to test for association. 
#' @param sample_key String denoting the name of the sample-level identifier (e.g. DonorID). 
#' @param graph A KNN graph in sparseMatrix format or the name of such a graph
#' stored in `colPair(sce)`. 
#' @param batches A character string. Name of batch variable. Currently only one 
#' categorical variable allowed. 
#' @param covs A character string or vector. Name(s) of other (numerical) covariates 
#' to control for. 
#' @param verbose Logical controlling verbosity
#' @param ... Passed to [association()]
#' 
#' @return A SingleCellExperiment object with NAM embeddings stored in 
#' `reducedDim(sce, "CNA")`, results stored in `metadata(sce)$CNA`, and correlation
#' values stored in `colData(sce)` as `cna_ncorrs`, `cna_ncorrs_fdr05`, and 
#' `cna_ncorrs_fdr10`
#' 
#' 
#' @export 
association.SingleCellExperiment <- function(sce, 
                                             test_var, 
                                             sample_key, 
                                             # sample_vars, 
                                             graph, 
                                             batches = NULL, 
                                             covs = NULL, 
                                             verbose = TRUE, 
                                             ...) {
    
    ## (1) format data 
    covs_keep <- test_var
    if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
    if (!is.null(covs)) covs_keep <- c(covs_keep, covs)
    rcna_data <- createObject.SingleCellExperiment(sce, 
                                                   sample_key, 
                                                   sample_vars = covs_keep, 
                                                   graph)
    yvals <- rcna_data$samplem[[test_var]]
    if(is.character(yvals)) {
        warning("'test_var' points to a character vector. Converting it to numeric via factor", 
                immediate. = TRUE)
        yval <- as.numeric(as.factor(yvals))
    }
    
    ## (2) do association
    cna_res <- association(data = rcna_data, 
                           y = test_var, 
                           batches = batches, 
                           covs = covs,
                           verbose = verbose,
                           ...) 

    # cna_res$samplem_df = rcna_data$samplem
    
    ## (3) save results 
    SingleCellExperiment::reducedDim(sce, "CNA") <- LinearEmbeddingMatrix(
        sampleFactors = cna_res$NAM_embeddings, 
        featureLoadings = cna_res$NAM_loadings, 
        factorData = S4Vectors::DataFrame(svs = cna_res$NAM_svs),
        metadata = list()
    )
    # put the rest into metadata
    S4Vectors::metadata(sce)$CNA <- cna_res[!names(cna_res) %in% 
                                                c("NAM_embeddings", "NAM_loadings", "NAM_svs")]
    # corrs goes into colData for easy plotting
    sce$cna_ncorrs <- cna_res$ncorrs[colnames(sce), , drop=TRUE]
    
    ## NOTE: If threshold was NULL, then no cells passed the significance threshold 
    sce$cna_ncorrs_fdr05 <- rep(0, ncol(sce))
    if (!is.null(cna_res$fdr_5p_t)) {
        idx_passed <- which(abs(sce$cna_ncorrs) >= cna_res$fdr_5p_t)
        sce$cna_ncorrs_fdr05[idx_passed] <- sce$cna_ncorrs[idx_passed]
    }
    
    sce$cna_ncorrs_fdr10 <- rep(0, ncol(sce))
    if (!is.null(cna_res$fdr_10p_t)) {
        idx_passed <- which(abs(sce$cna_ncorrs) >= cna_res$fdr_10p_t)
        sce$cna_ncorrs_fdr10[idx_passed] <- sce$cna_ncorrs[idx_passed]
    }
    
    return(sce)
}