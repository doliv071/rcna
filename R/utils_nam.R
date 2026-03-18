#' Perform a single diffusion mapping step.
#' 
#' @param a A sparseMatrix M by M (cells x cells) adjacency matrix (connectivities) 
#' whose m,m'-th entry indicates the similarity between cells m and m' in the graph
#' @param s A sparseMatrix M by N (cells x samples) indicator matrix assigning 
#' the m-th cell to the n-th sample from which it was drawn.
#' 
#' @note This function ignores distances and only uses unweighted connectivities
#' 
#' @returns A M by N sparseMatrix.
#' 
#' @keywords internal
diffuseStep <- function(a, s) {
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(res) 
}

#' Calculate kurtosis accounting for batch
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param batch A factor containing batch information
#' 
#' @returns a vector of kurtosis values for each column (cell) in the NAM
#' 
#' @keywords internal
batchKurtosis <- function(NAM, batch) {
    stopifnot(is.factor(batch))
    stopifnot(isa(NAM, "Matrix"))
    
    bI <- Matrix::fac2sparse(batch)
    # mean batch NAM
    bNAM <- Matrix::t(bI %*% NAM) %*% Matrix::Diagonal(x = 1 / Matrix::rowSums(bI))
    res <- Kurtosis(bNAM, margin = 1)
    return(res)
}

#' Filter NAM based on batch level kurtosis
#'
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param batch A factor containing batch information or NULL. If NULL, no
#' batch QC is performed (this is for pass-through). \cr
#' Default: NULL
#' @param min.threshold A numeric scalar defining the minimum value which the 
#' threshold should take.\cr
#' Default: 6
#' @param verbose Logical. Controls verbosity.\cr
#' Default: FALSE
#' 
#' @note min.threshold is 6 here, to keep functionality the same as initially 
#' written, but there is no mention of this in the paper.
#'
#' @returns A list containing the filtered NAM (NAM) and an index vector indicating
#' which cells were kept (keep).
#' 
#' @keywords internal
qcNAM <- function(NAM, 
                  batch = NULL, 
                  min.threshold = 6,
                  verbose = FALSE) {

    if (is.null(batch) || length(unique(batch)) == 1) {
        if(verbose) message('Only one unique batch supplied to qc')
        keep <- rep(TRUE, ncol(NAM))
        res <- list(NAM = NAM, keep = keep)
    } else {
        medkurt <- Kurtosis(NAM, margin = 1) |> median()
        kurtoses <- batchKurtosis(NAM, batch)
        # TODO: add batch level kurtosis message
        # if(verbose) message("")
        
        if(verbose) message("Median batch kurtosis: ", round(median(kurtoses), 3))
        threshold <- max(min.threshold, 2*median(kurtoses))
        if(verbose) message("Throwing out neighborhoods with batch kurtosis >= ", threshold)
        keep <- which(kurtoses < threshold)
        if(verbose) message("keeping ", length(keep), " neighborhoods")
        # if(verbose) message("Median batch kurtosis after filtering: ", 
        #                     round(kurtoses, 3))
        res <- list(NAM = NAM[, keep, drop = FALSE], keep = keep)
    }
    return(res) 
}

#' Perform batch correction of a NAM with OLS regression
#'
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param cov.mat a matrix of covariates. \cr
#' Default: NULL
#' @param batch a factor containing the batch to correct. \cr
#' Default: NULL
#'
#' @returns A list containing the corrected NAM, the annihilator matrix (M), 
#' and OLS degrees of freedom consumed (r)
#' 
#' @keywords internal
olsNAM <- function(NAM, 
                   cov.mat = NULL, 
                   batch = NULL){
    if(!is.null(batch)){
        B <- Matrix::sparse.model.matrix(~0+batch)
    } 
    if(!is.null(cov.mat)){
        C <- Scale(as(cov.mat, "sparseMatrix"))
        
    }
    X <- cbind(B, C)
    N <- nrow(NAM)
    I <- Matrix::Diagonal(n = N)
    # catch if we got here without a batch or covariate matrix
    if(is.null(X)){
        warning("no batch or cov.mat supplied.", immediate. = TRUE)
        res <- list(NAM_ = NAM, 
                    M = I,
                    r = 0)
    } else {
        dfs <- ncol(X)
        H <- X %*% Matrix::solve(Matrix::crossprod(X, X), Matrix::t(X))
        M <- I - H
        NAM_ <- M %*% NAM 
        colnames(NAM_) <- colnames(NAM)
        rownames(NAM_) <- rownames(NAM)
        res <- list(NAM_ = NAM_, 
                    M = M,
                    r = dfs)
    }
    return(res)
}

#' Perform batch correction of a NAM using (partial) ridge regression
#'
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param cov.mat a matrix of covariates. \cr
#' Default: NULL
#' @param batch a factor containing the batch to correct. \cr
#' Default: NULL
#' @param ridges A vector of ridges to try or NULL. If NULL ridges are calculated
#' automatically centered on `svd(X)$d[1]^2`. \cr
#' Default: NULL
#' @param partial A character string or NULL. If NULL, ridge regression is performed
#' on all variables. If a character, specify which set of variables ridge regression
#' should be used on, the other set does not get penalized. Must be one of "batches" or
#' "covariates". \cr 
#' Default: NULL
#'
#' @returns A list containing the corrected NAM, the Hat matrix (M), and the 
#' effective degree of freedom (`sum(diag(H))` from the ridge hat matrix) (r)
#' 
#' @keywords internal
ridgeNAM <- function(NAM, 
                     cov.mat = NULL, 
                     batch = NULL, 
                     ridges = NULL,
                     partial = NULL){
    if(!is.null(partial)){
        partial <- match.arg(partial, c("batches", "covariates"))
    }
    if(!is.null(partial) && any(c(is.null(cov.mat), is.null(batch)))){
        warning("'partial' is not null, but only 'cov.mat' or 'batch' supplied.",
                " Using standard ridge regression.")
    }
    dfs <- 0
    if(!is.null(batch) && !is.null(cov.mat)){
        B <- Matrix::sparse.model.matrix(~0+batch)
        C <- as(cov.mat, "sparseMatrix")
        X <- Scale(cbind(B, C))
        
        if(is.null(partial)){ 
            L <- Matrix::Diagonal(n = ncol(X))
        } else if(partial == "batches"){ # only do ridge on batches, fit OLS to covariates
            L <- Matrix::Diagonal(x = c(rep(1, ncol(B)), rep(0, ncol(X) - ncol(B))))
            dfs <- ncol(C)
        } else { # only do ridge on covariates, fit OLS to batches
            L <- Matrix::Diagonal(x = c(rep(0, ncol(B)), rep(1, ncol(X) - ncol(B))))
            dfs <- ncol(B)
        }
    } else if(!is.null(batch)){
        X <- Matrix::sparse.model.matrix(~0+batch) |> Scale()
        L <- Matrix::Diagonal(n = ncol(X))
    } else if (!is.null(cov.mat)) {
        X <- Scale(as(cov.mat, "sparseMatrix"))
        L <- Matrix::Diagonal(n = ncol(X))
    }
    N <- nrow(NAM)
    I <- Matrix::Diagonal(n = N)
    if(is.null(ridges)){
        sv_max <- svd(X)$d[1]^2
        ridges <- sv_max * 10^seq(3, -3, length.out = 20)
        # ridges <- c(10^seq(-5, 5, length.out = 20), sqrt(.Machine$double.eps))
    } 
    # find the best ridge penalty
    gcv <- vapply(ridges, function(l) {
        H <- X %*% Matrix::solve(Matrix::crossprod(X, X) + l * L, Matrix::t(X))
        resids <- (I - H) %*% NAM
        num <- (1/N) * sum(resids^2)
        den <- ((1/N) * (1 - sum(Matrix::diag(H))))^2
        num/den
    }, numeric(1)) |> round(3)
    
    lambda <- ridges[which.min(gcv)]
    H <- X %*% Matrix::solve(Matrix::crossprod(X, X) + lambda * L, Matrix::t(X))
    M <- I - H
    NAM_ <- M %*% NAM
    colnames(NAM_) <- colnames(NAM)
    rownames(NAM_) <- rownames(NAM)
    kurtoses <- batchKurtosis(NAM_, batch) |> median()
    message("with ridge ", lambda, " median batch kurtosis = ", kurtoses)
    dfs <- dfs + sum(Matrix::diag(H))
    return(list(NAM_ = NAM_, 
                M = M,
                r = dfs))
}

#' remove batch effects (including continuous covariates)
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param covs.mat A matrix of covariates.\cr
#' Default NULL
#' @param batch A factor containing batch information.\cr
#' Default NULL
#' @param method A character string specifying what method to use for batch
#' and covariate correction. Options are: c("ridge", "ols"), Setting to "ridge"
#' will reproduce original `rcna` behavior.
#' @param partial.by A character string. One of NULL, "batches", or "covariates". 
#' If `method == "ridge"` then ridge penalties are only applied to those variables.
#' If NULL then all variables are penalized. Setting to "batches" will reproduce
#' original `rcna` behavior.\cr
#' Default: NULL
#' @param ridges A numeric vector of ridge penalties to try. If NULL, a default
#' set of ridges centered on `svd(X)$d[1]^2` is used.\cr
#' Default NULL
#' 
#' @note assumes that covariates are not categorical
#' 
#' @returns A list containing the residuals NAM (NAM_), the annihilator matrix 
#' used to remove unwanted effects (M), and the approximate degrees of freedom
#' consumed (r).
#' 
#' @keywords internal
residNAM <- function(NAM, 
                     covs.mat = NULL, 
                     batch = NULL, 
                     method = c("ridge", "ols"),
                     partial.by = c("batches", "covariates"),
                     ridges = NULL, 
                     verbose = FALSE) {
    method <- match.arg(method, choices = c("ridge", "ols"))
    
    NAM_ <- Scale(NAM, center = TRUE, scale = FALSE)
    if(!is.null(batch) || !is.null(covs.mat)){
        if(method == "ols"){
            NAMres <- olsNAM(NAM_, covs.mat, batch)
        } else {
            if(!is.null(partial.by) && (!is.null(batch) && !is.null(covs.mat))){
                partial.by <- match.arg(partial.by, choices = c("batches", "covariates"))
            } else {
                partial.by <- NULL
            }
            NAMres <- ridgeNAM(NAM_, covs.mat, batch, ridges, partial = partial.by)
        }
    } else {
        # pass through to generate proper results list. 
        I <- Matrix::Diagonal(n = nrow(NAM))
        NAMres <- list(NAM_ = NAM, 
                       M = I,
                       r = 0)
    }
    NAMres$NAM_ <- Scale(NAMres$NAM_, center = FALSE, scale = TRUE)
    return(NAMres)
}

#' Compute the singular-value decomposition of a NAM matrix.
#' 
#' @param NAM An M x N (Cells x Samples) neighborhood association matrix as 
#' calculated by [buildNAM()]
#' @param n.pcs maximum number of components to use in the SVD. Passed to the `k`
#' parameter of [RSpectra::svds()]. If NULL or greater than `0.5 * min(dim(NAM))`
#' then the full svd is calculated. \cr
#' Default: NULL
#' 
#' @note svs are actually eigenvalues, not SVs. I squared them to be consistent 
#' with python code. 
#' @note svdNAM scales NAM columns, so this is correlation, not covariance
#' 
#' @returns A list containing the U, D, and V' components of the SVD.
#' 
#' @keywords internal
svdNAM <- function(NAM, n.pcs = NULL) {
    if (is.null(n.pcs) || n.pcs > .5 * min(dim(NAM))) {
        svd_res <- svd(NAM)
    } else {
        svd_res <- RSpectra::svds(NAM, k = n.pcs)
    }
    
    U <- svd_res$u[, seq_len(n.pcs)]
    colnames(U) <- paste0('PC', seq_len(n.pcs))
    rownames(U) <- rownames(NAM)
    V <- svd_res$v[, seq_len(n.pcs)]
    colnames(V) <- paste0('PC', seq_len(n.pcs))
    rownames(V) <- colnames(NAM)
    
    res <- list(U = U, 
                svs = svd_res$d^2, 
                V = V)
    return(res)
}

#' Calculate Neighborhood Association Matrix (NAM) using diffusion mapping
#' 
#' @param data A rcna data object (list). See [nam()] or [association()] for details.
#' @param n.steps A numeric, the number of steps to take during the random walk. 
#' If NULL, then `kurt.delta` is checked at each step to decide whether or not 
#' to take the next step. \cr
#' Default: NULL
#' @param min.steps A numeric value controlling the minimum number of steps to take. 
#' Only used if n.steps is NULL. \cr
#' Default: 3
#' @param max.steps A numeric value controlling the maximum number of steps to take. 
#' Only used if n.steps is NULL. \cr
#' Default: 15
#' @param kurt.delta The minimum change in kurtosis between steps under which
#' the walk stops (if min.steps have already been made). Only used if n.steps is NULL.\cr
#' Default: 3
#'
#' @note The original code had a bug here where `prevmedkurt <- medkurt` was 
#' updated before calculating `prevmedkurt - medkurt` and so only a single diffusion
#' step was performed.
#' 
#' @returns A sparseMatrix containing the NAM (A cell by sample neighborhood
#' association matrix).
#' 
#' @keywords internal
buildNAM <- function(data, 
                     n.steps = NULL, 
                     min.steps = 3L,
                     max.steps = 15L,
                     kurt.delta = 3, 
                     verbose = FALSE) {
    
    f <- as.formula(paste0("~ 0 + ", data$samplem_key))    
    s <- Matrix::sparse.model.matrix(f, data$obs)
    colnames(s) <- gsub(data$samplem_key, '\\1', colnames(s))
    rownames(s) <- data$obs[[data$obs_key]]
    s <- s[, data$samplem[[data$samplem_key]]] ## Necessary? 
    
    prevmedkurt <- Inf
    ## CHECK: number of iterations matches 
    for (i in seq_len(max.steps)) {
        s <- diffuseStep(data$connectivities, s)
        medkurt <- median(Kurtosis(propTable(s, margin = 2), margin = 1))
        # if(verbose) message("Median kurtosis: ", round(medkurt, 3))
        # if n.steps is null then take at least min.steps steps
        if (is.null(n.steps)) {
            kd <- prevmedkurt - medkurt
            if(verbose) message("Median kurtosis = ", round(medkurt, 3), 
                                " Kurtosis change = ", round(kd, 3), " at step ", i)
            if (kd < kurt.delta & i > min.steps) {
                message("stopping after ", i, " steps")
                break 
            }
            prevmedkurt <- medkurt
        } else if (i == n.steps) {
            break
        }
    }  
    sfin <- propTable(s, margin = 2)
    medkurt <- median(Kurtosis(sfin, margin = 1))
    if(verbose) message("Final kurtosis: ", round(medkurt, 3))
    snorm <- Matrix::t(sfin)
    rownames(snorm) <- data$samplem[[data$samplem_key]]
    colnames(snorm) <- data$obs[[data$obs_key]]
    return(snorm)
}

#' Build and decompose a neighborhood abundance matrix 
#' 
#' @param data a list containing: 
#' samplem = (sample-level metadata), 
#' obs = (cell-level metadata), 
#' connectivities  = (sparse symmetric weighted adjacency matrix),
#' samplem_key = (character string indicating the column in samplem uniquely identifying samples),
#' obs_key = (character string indicating the column in obs uniquely identifying cells),
#' N = nrow(samplem_df).
#' @param y A character string specifying the column in `data$samplem` where 
#' the variable of interest is stored.
#' @param batches A character string denoting the column in `data$samplem` where
#' batch information is stored. Only a single batch variable is allowed. \cr
#' Default: NULL 
#' @param covs A character string or vector denoting the column(s) `data$samplem`
#' where continuous covariate information is stored.\cr
#' Default: NULL
#' @param n.steps Numeric scalar controlling the number of steps to take during 
#' the random walk.\cr
#' Default: NULL
#' @param min.steps Numeric scalar. Minimum number of steps to take on the random
#' walk.\cr
#' Default: 3L
#' @param max.steps Numeric scalar. Maximum number of steps to take on the random
#' walk.\cr
#' Default: 15L
#' @param kurtosis.delta Numeric scalar. Minimum kurtosis change below which the
#' random walk can stop (after travelling at least `min.steps`).\cr
#' Default: 3
#' @param min.batch.kurtosis Sets the minimum threshold below which cells are kept
#' during QC. This ensures that when kurtosis is not extreme between batches
#' then cells are not filtered. `threshold <- max(min.batch.kurtosis, 2*median(kurtoses))`.\cr
#' Default: 6
#' @param max.frac.pcs The number of PCs to calculate for SVD. The minimum is 10 
#' and the maximum is the number of samples minus 1. \cr
#' Default: 0.15
#' @param method A character string specifying what method to use for batch
#' and covariate correction. Options are: c("ridge", "ols"), Setting to "ridge"
#' will reproduce original behavior.
#' @param partial.by A character string. One of NULL, "batches", or "covariates". If 
#' `method == "ridge"` then ridge penalties are only applied to those variables.
#' If null then all variables are penalized. Setting to "batches" will reproduce
#' original behavior.\cr
#' Default: NULL
#' @param ridges A numeric vector of ridges to try during ridge regression. When
#' NULL ridges are calculated automatically. \cr
#' Default: NULL
#' @param suffix A character scalar for optionally adding a suffix to the output
#' components. \cr
#' Default: ""
#' @param verbose Logical scalar controlling verbosity. \cr
#' Default: FALSE
#' @param filter.samples STUB. Not used currently.
#' 
#' @return A named list. All keys are optionally suffixed by \code{suffix}.
#' \describe{
#'   \item{y}{Numeric vector of the (possibly coerced) response variable.}
#'   \item{raw.NAM.T}{sparseMatrix. Transposed raw NAM (cells × samples) before
#'     QC or residualisation.}
#'   \item{qc.NAM.T}{sparseMatrix. Transposed NAM after batch-kurtosis QC.}
#'   \item{resid.NAM.T}{sparseMatrix. Transposed NAM after residualisation.}
#'   \item{_M}{The annihilator matrix from [residNAM()].}
#'   \item{_r}{Numeric. Degrees of freedom consumed by the correction.}
#'   \item{NAM_sampleXpc}{Matrix (samples × PCs): left singular vectors.}
#'   \item{NAM_svs}{Numeric vector of squared singular values.}
#'   \item{NAM_varexp}{Numeric vector of variance explained per PC.}
#'   \item{NAM_nbhdXpc}{Matrix (cells × PCs): right singular vectors.}
#'   \item{keptcells}{Integer index vector of retained cell columns after QC.}
#'   \item{_batches}{Factor of batch labels, or \code{NULL}.}
#'   \item{covariates}{sparseMatrix of covariates, or \code{NA_real_} if none.}
#' }
#' 
#' @export 
nam <- function(data, y,
                batches = NULL, 
                covs = NULL, 
                # random walk parameters passed to buildNAM
                n.steps = NULL, 
                min.steps = 3L,
                max.steps = 15L,
                kurtosis.delta = 3,
                # minimum batch kurtosis passed to qcNAM
                min.batch.kurtosis = 6,
                max.frac.pcs = 0.15, 
                # passed to residNAM
                method = c("ridge", "ols"),
                partial.by = NULL,
                ridges = NULL,
                # general
                suffix = '',
                filter.samples = NULL,
                verbose = FALSE) {
    if(missing(data) || missing(y)){
        stop("missing parameters with no defaults.")
    }
    stopifnot(all(c("samplem", "obs", "connectivities", 
                    "samplem_key", "obs_key", "N") %in% names(data)))
    stopifnot(is.null(batches) || (is.character(batches) && length(batches) == 1))
    stopifnot(is.null(covs) || (is.character(covs)))
    # passed to buildNAM
    stopifnot(is.null(n.steps) || (is.numeric(n.steps) && length(n.steps) == 1))
    stopifnot(is.numeric(min.steps) && length(min.steps) == 1)
    stopifnot(is.numeric(max.steps) && length(max.steps) == 1)
    stopifnot(is.numeric(kurtosis.delta) && length(kurtosis.delta) == 1)
    # passed to qcNAM
    stopifnot(is.numeric(min.batch.kurtosis) && length(min.batch.kurtosis) == 1)
    stopifnot(is.numeric(max.frac.pcs) && length(max.frac.pcs) == 1)
    # passed to residNAM
    stopifnot(is.character(method) && length(method) == 1)
    method <- match.arg(method, c("ridge", "ols"))
    stopifnot(is.character(method) && length(method) == 1)
    stopifnot(is.null(partial.by) || (is.character(partial.by) && length(n.steps) == 1))
    partial.by <- match.arg(partial.by, c("batches", "covariates"))
    stopifnot(is.null(ridges) || is.numeric(ridges))
    # general
    stopifnot(is.character(suffix) && length(suffix) == 1)
    
    if(!is.null(n.steps) && n.steps > max.steps){
        warning("'n.steps' is larger than 'max.steps', updating 'max.steps' to allow", 
                " at least 'n.steps'.", call. = FALSE, immediate. = TRUE)
        max.steps <- n.steps
    }
    
    res <- list()
    ## TODO: Filtering should happen at top level, either during data object
    ##       build or as a helper to manipulate existing data object
    
    if (is.null(batches)) {
        batch_vec <- NULL
    } else {
        if(!batches %in% colnames(data$samplem)){
            stop("Could not find batches: ", batches, " in `colnames(data$samplem)`")
        }
        # batches must be categorical or factor
        # TODO: handle multiple batches
        batch_vec <- dplyr::pull(data$samplem, dplyr::one_of(batches)) |> 
            as.factor()
    }
    if (is.null(covs)) {
        cov_mat <- NULL
    } else {
        # covariates must be numeric or factors
        # TODO: actually check that they are numeric
        cov_mat <- data$samplem[, covs, drop = FALSE] |> 
            as("sparseMatrix")
        # # covariates must be numeric or factors
        # cov_df <- data$samplem[, covs, drop = FALSE] 
        # cform <- paste0("~ 0 + ", paste0(colnames(cov_df), collapse = " + ")) |> 
        #     as.formula()
        # cov_mat <- Matrix::sparse.model.matrix(cform, cov_df)
    }
    
    if (verbose) message('Constructing NAM')
    NAM <- buildNAM(data, 
                    n.steps = n.steps, 
                    min.steps = min.steps,
                    max.steps = max.steps, 
                    kurt.delta = kurtosis.delta, 
                    verbose = verbose)
    if (verbose) message('QC-ing NAM')
    res_qc_nam <- qcNAM(NAM, 
                        batch = batch_vec, 
                        min.threshold = min.batch.kurtosis, 
                        verbose = verbose) 
    # y is returned in order to avoid recalculating NAM if not desired.
    # TODO: handle multiple y?
    y <- data$samplem[, y, drop = TRUE]
    if(is.character(y)){
        y <- as.numeric(as.factor(y))
    } else {
        y <- as.numeric(y)
    }
    
    
    ## (3) Decompose NAM 
    ## TODO: check if double brackets appropriate for multiple covs and/or batches
    ## TODO: check with Y&L if covs should be numerical and batches categorical
    ## NOTE: don't really need a separate SVD function, since it can be done in one line
    if (verbose) message('Residualizing NAM')
    res_resid_nam <- residNAM(res_qc_nam$NAM, 
                              covs.mat = cov_mat, 
                              batch = batch_vec, 
                              method = method,
                              partial.by = partial.by,
                              ridges = ridges,
                              verbose = verbose)
    if (verbose) message('Decomposing NAM')
    n_pcs <- max(10, ceiling(max.frac.pcs * nrow(data$samplem)))
    n_pcs <- min(n_pcs, nrow(data$samplem) - 1) ## make sure you don't compute all SVs    
    res_svd_nam <- svdNAM(res_resid_nam$NAM_, n.pcs = n_pcs)
    
    res[['y']] <- y
    res[['raw.NAM.T']] <- Matrix::t(NAM)
    res[['qc.NAM.T']] <- Matrix::t(res_qc_nam$NAM)
    res[['resid.NAM.T']] <- Matrix::t(res_resid_nam$NAM_)
    res[['_M']] <- res_resid_nam$M
    # TODO rename to df (degrees of freedom)
    res[['_r']] <- res_resid_nam$r
    res[['NAM_sampleXpc']] <- res_svd_nam$U
    res[['NAM_svs']] <- res_svd_nam$svs
    res[['NAM_varexp']] <- res_svd_nam$svs / nrow(res_svd_nam$U) / nrow(res_svd_nam$V)
    res[['NAM_nbhdXpc']] <- res_svd_nam$V
    res[['keptcells']] <- res_qc_nam[[2]]
    res[['_batches']] <- batch_vec
    res[['covariates']] <- cov_mat
    
    names(res) <- paste0(names(res), suffix)
    
    return(res)
}
