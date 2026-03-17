#' Fit y to the data
#' 
#' @param q A vector of values to fit (a la y)
#' @param k Number of components of U to fit
#' @param U the left singular vector matrix of SVD
#' @returns a list containing fitted values (qhat) and regression coefficients (beta) 
#'
#' @keywords internal
regress <- function(q, k, U) {
    Xpc <- U[, 1:k]
    beta <- crossprod(Xpc, q)
    qhat <- Xpc %*% beta
    res <- list(qhat = qhat,
                beta = beta)
    return(res)
}

#' Calculate the F-test p-value
#' 
#' @param yhat y fit to the data
#' @param ycond y corrected for batch and covariates
#' @param k Number of components of U to test
#' @param n number of samples
#' @param r degrees of freedom
#' @param log.p A logical whether the log transformed p-value value should be 
#' calculated.\cr
#' Default: FALSE
#' 
#' @returns A list containing p-value and r^2 for the F-test
#' 
#' @keywords internal
calcStats <- function(yhat, ycond, k, n, r, log.p = FALSE) {
    ssefull <- crossprod(yhat - ycond)
    ssered <- crossprod(ycond)
    deltasse <-  ssered - ssefull
    f <- (deltasse / k) / (ssefull / n)
    if(!log.p){
        p <- pf(f, k, n-(1+r+k), lower.tail = FALSE)
    } else {
        # N.B. this is an approximation of 1-p when p gets too small to calculate
        #      in 64bit space? For small f/df1/df2 the drift is noticeable though
        p <- -pf(f, k, n-(1+r+k), log.p = TRUE)
    }
    r2 <- 1 - ssefull/ssered
    res <- list(p = p, r2 = r2)
    return(res)
}

#' Calculate the p-value for all ks and return the minimum
#' 
#' @param M the annihilator matrix (M returned by [residNAM])
#' @param y A vector containing the variable of interest
#' @param ks A vector of ks to test.
#' @param U a matrix whose columns contain the left singular vectors of x 
#' (U returned by [svdNAM])
#' @param r the degrees of freedom (r returned by [residNAM])
#' @param use.logp A logical whether the log transformed p-value value should be 
#' calculated.\cr
#' Default: FALSE
#' 
#' @returns A list containing the k, p, and r^2 for k which produces the minimum
#' p-value.
#' 
#' @keywords internal
minpStats <- function(M, y, ks, U, r, use.logp = FALSE) {
    n <- length(y)
    zcond <- M %*% y
    zcond <- scale(zcond, center = FALSE, scale = TRUE)
    qhats <- lapply(ks, function(k) regress(zcond, k, U)$qhat)
    .tmp <- lapply(seq_along(ks), function(i) calcStats(qhats[[i]], zcond, ks[i], n, r, log.p = use.logp))
    ps <- vapply(.tmp, \(x) x$p, numeric(1))
    r2s <- vapply(.tmp, \(x) x$r2, numeric(1))
    k_ <- which.min(ps)
    res <- list(k = ks[k_], 
                p = ps[k_], 
                r2 = r2s[k_])
    return(res)
}

#' Calculate global and (optionally) local association tests.
#'
#' @param NAMsvd The list output from the [nam()] function.
#' @param y A vector with contrast variable value to be tested for association.
#' @param batches_vec A factor or numeric vector of batches to adjust for.
#' @param ks A numeric scalar selecting the number of components of the SVD
#' to test for global association. If null 4 values for k are selected between
#' n/50 and n/5. \cr
#' Default: NULL
#' @param Nnull A numeric specifying the number of null permutations.\cr
#' Default: 1000
#' @param force_permute_all A logical controlling whether permutations of `y` 
#' should preserve batch information. \cr
#' Default: FALSE 
#' @param local_test A logical, whether or not to perform local test of neighborhood
#' correlations. \cr
#' Default: TRUE
#' @param use.logp A logical whether the log transformed p-value value should be 
#' calculated. Set to TRUE if your FDRs are incredibly small or 0. \cr
#' Default: FALSE
#' @param seed A numeric seed to set. Set if you want to repeat exact permutations
#' as a previous run. If NULL, a random seed is chosen. \cr
#' Default: NULL
#'
#' @returns A list. 
#' list(p = pfinal, 
#' minps_null = minps_null,
#' k = k,
#' ncorrs = ncorrs, 
#' fdrs = fdrs,
#' fdr_5p_t = fdr_5p_t, 
#' fdr_10p_t = fdr_10p_t,
#' yhat = yhat, 
#' ycond = ycond,
#' ks = ks, 
#' beta = beta,
#' r2 = r2, 
#' r2_perpc = r2_perpc,
#' nullr2_mean = mean(nullr2s), 
#' nullr2_std = sd(nullr2s))
#' 
#' @keywords internal
innerAssociation <- function(NAMsvd, y, batches_vec, 
                             ks = NULL, 
                             Nnull = 1000, 
                             force_permute_all = FALSE, 
                             allow_duplicate_perms = TRUE,
                             local_test = TRUE, 
                             use.logp = FALSE,
                             seed = NULL) {
    if (is.null(seed)) {
        seed <- sample(1e6, 1)
    }
    set.seed(seed)
    if (force_permute_all) {
        batches_vec <- rep(1L, length(y))
    } 
    
    # prep data
    U <- NAMsvd[["NAM_sampleXpc"]]
    sv <- NAMsvd[["NAM_svs"]]
    V <- NAMsvd[["NAM_nbhdXpc"]]
    M <- NAMsvd[["_M"]] 
    r <- NAMsvd[["_r"]]
    y <- scale(y)
    n <- length(y)
    
    if (is.null(ks)) {
        # unique handles situations where n < 16
        ks <- seq(n/50, n/5, length.out = 4) |> ceiling() |> unique()
    } else if(length(ks) == 1){
        if(ks < dim(U)[2]){
            ks <- seq(1, ks, 1)
        } else {
            warning("user defined ks is too large. setting it to maximum allowable: ", 
                    dim(U)[2] - 1)
            ks <- seq(1, dim(U)[2] - 1, 1)
        }
    } 
    
    # get non-null f-test p-value
    mp <- minpStats(M, y, ks, U, r, use.logp = use.logp)
    k <- mp$k 
    p <- mp$p 
    r2 <- mp$r2
    if (k == max(ks)) {
        warning("data supported use of ", k, 
                " NAM PCs, which is the maximum considered.", 
                " Consider allowing more PCs by using the 'ks' argument.")        
    }
    
    # compute coefficients and r2 with chosen model
    ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
    rg <- regress(ycond, k, U)
    yhat <- rg$qhat
    beta <- rg$beta
    r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))**2
    
    # get neighborhood scores with chosen model
    ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta / n)
    rownames(ncorrs) <- rownames(V)
    
    # compute final p-value using Nnull null f-test p-values
    y_null <- conditional_permutation(batches_vec, y, Nnull, 
                                      duplicates.ok = allow_duplicate_perms, 
                                      seed = seed)
    .tmp <- apply(y_null, 2, \(z) minpStats(M, z, ks, U, r, use.logp = use.logp))
    minps_null <- vapply(.tmp, \(x) x$p, numeric(1))
    nullr2s <- vapply(.tmp, \(x) x$r2, numeric(1))
    # add sqrt(.Machine$double.eps) for floating point maths and add 1 to avoid pfinal = 0
    pfinal <- (sum(minps_null <= p + sqrt(.Machine$double.eps)) + 1)/(Nnull + 1)
    if (sum(minps_null <= p + sqrt(.Machine$double.eps)) == 0) {
        warning("global association p-value attained minimal possible value.", 
                "Consider increasing Nnull", immediate. = TRUE)
    }
    
    # get neighborhood fdrs if requested
    fdrs <- NULL
    fdr_5p_t <- NULL 
    fdr_10p_t <- NULL
    
    if (local_test) {
        message('computing neighborhood-level FDRs')
        Nnull <- min(1000, ncol(y_null))
        y_null <- y_null[, 1:Nnull]
        ycond_null <- scale(M %*% y_null, center = FALSE, scale = TRUE)
        gamma_null <- crossprod(U[, 1:k], ycond_null)
        ncorrs_null <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_null / n)))
        # use both ncorrs and ncorrs_null to make sure we observe full range
        maxcorr <- max(c(abs(ncorrs), abs(ncorrs_null))) + sqrt(.Machine$double.eps)
        # using bins way outside the range generates NaN since both
        # NULL and data have 0 counts in the bins. However, starting at 0
        # avoids shenanigans in empirical_fdrs()
        fdr_thresholds <- seq(0, maxcorr, maxcorr/400) # maxcorr/4
        fdrs <- empirical_fdrs(ncorrs, ncorrs_null, fdr_thresholds) # was fdr_vals
        # TODO: filter fdrs where num_detected = 0 ?
        # starting from last until we observe the first non-0?
        
        # find maximal FDR<5% and FDR<10% sets
        if (min(fdrs$fdr) > 0.05) {        
            fdr_5p_t <- NULL
        } else {
            fdr_5p_t <- min(subset(fdrs, fdr < 0.05)$threshold)        
        }
        if (min(fdrs$fdr) > 0.10) {        
            fdr_10p_t <- NULL
        } else {
            fdr_10p_t <- min(subset(fdrs, fdr < 0.1)$threshold)
        }
    }
    
    res <- list(p = pfinal, 
                minps_null = minps_null,
                k = k,
                ncorrs = ncorrs, 
                fdrs = fdrs,
                fdr_5p_t = fdr_5p_t, 
                fdr_10p_t = fdr_10p_t,
                yhat = yhat, 
                ycond = ycond,
                ks = ks, 
                beta = beta,
                r2 = r2, 
                r2_perpc = r2_perpc,
                nullr2_mean = mean(nullr2s), 
                nullr2_std = sd(nullr2s),
                seed = seed)
    
    return(res)
}

#' Main function to perform CNA association analysis
#' 
#' @param data a list containing: 
#' samplem = (sample-level metadata), 
#' obs = (cell-level metadata), 
#' connectivities  = (sparse symmetric weighted adjacency matrix),
#' samplem_key = (character string indicating the column in samplem uniquely identifying samples),
#' obs_key = (character string indicating the column in obs uniquely identifying cells),
#' N = nrow(samplem_df).
#' @param y A character string specifying the column in `samplem` containing
#' the variable of interest.
#' @param batches A character string specifying the column in `samplem` 
#' containing the batch variable. Only a single batch variable is currently supported. \cr
#' Default: NULL
#' @param covs A character string or vector specifying the column(s) in `samplem` 
#' containing the covariate variables.\cr
#' Default: NULL
#' @param n.steps Number of steps to take during the random walk. If specified then
#' exactly this many steps is taken on the random walk. \cr
#' Default: NULL
#' @param suffix A character string to be appended for unknown reasons.\cr
#' Default: ''
#' @param return.nam Logical controlling whether or not to return the NAM. \cr
#' Default: FALSE
#' @param filter.samples STUB. Currently ignored.\cr
#' @param Ks One of NULL, a numeric scalar, or a numeric vector. If NULL, then 
#' `seq(n/50, n/5, length.out = 4)` are checked. If a numeric scalar then 1:Ks are
#' checked. If a numeric vector, then all supplied values are checked. \cr
#' Default: NULL
#' @param N.nulls .\cr
#' Default: 1000
#' @param force.permute.all .\cr
#' Default: FALSE
#' @param allow.duplicate.perms .\cr
#' Default: TRUE
#' @param local.test .\cr
#' Default: TRUE
#' @param use.logp A logical whether the log transformed p-value value should be 
#' calculated. Set to TRUE if your FDRs are incredibly small or 0. \cr
#' Default: FALSE
#' @param seed .\cr
#' Default: NULL
#' @param verbose Logical controlling the verbosity of the function.\cr
#' Default: FALSE
#' @param ... Additional parameters passed to [nam()]
#' 
#' @return A list containing p, nullminps, k, ncorrs, fdrs, fdr_5p_t, fdr_10p_t, 
#' yhat, ycond, ks, beta, r2, r2_perpc, nullr2_mean, and nullr2_std. If 
#' `return.nam = TRUE`, then additionally NAM_embeddings, NAM_loadings, and NAM_svs
#' 
#' @export 
association <- function(data, nam.result, y, 
                        batches = NULL, 
                        covs = NULL, 
                        n.steps = NULL, 
                        suffix = '',
                        return.nam = FALSE, 
                        filter.samples = NULL,
                        # passed to inner association
                        Ks = NULL, 
                        N.nulls = 1000, 
                        force.permute.all = FALSE, 
                        allow.duplicate.perms = TRUE,
                        local.test = TRUE, 
                        use.logp = FALSE,
                        seed = NULL,
                        verbose = TRUE, 
                        ...) {
    if(missing(data) && missing(nam.result)){
        stop("One of 'data' or 'nam.result' must be specified.")
    } else if(!missing(data) && !missing(nam.result)){
        warning("both 'data' and 'nam.result' were supplied,", 
                " using nam.result for association test.", 
                immediate. = TRUE, call. = FALSE)
    }
    if(!missing(data)){
        stopifnot(batches %in% colnames(data$samplem))
        stopifnot(all(covs %in% colnames(data$samplem)))
    }
    # TODO: check NAs in batches and covariates.
    stopifnot(length(batches) == 1)
    if(!is.null(n.steps)){
        stopifnot(length(n.steps) == 1 && is.numeric(n.steps))
    }
    stopifnot(length(suffix) == 1 && is.character(suffix))
    stopifnot(length(return.nam) == 1 && is.logical(return.nam))
    stopifnot(length(N.nulls) == 1 && is.numeric(N.nulls))
    
    stopifnot(length(force.permute.all) == 1 && is.logical(force.permute.all))
    stopifnot(length(allow.duplicate.perms) == 1 && is.logical(allow.duplicate.perms))
    stopifnot(length(local.test) == 1 && is.logical(local.test))
    stopifnot(length(use.logp) == 1 && is.logical(use.logp))
    stopifnot(length(verbose) == 1 && is.logical(verbose))
    
    if(!is.null(Ks)){
        stopifnot(is.numeric(Ks))
    }
    
    ## TODO: add sample filtering 
    
    if (verbose) message('Build NAM PCs')
    if(missing(nam.result)){
        # formatting and error checking
        stopifnot(all(c("samplem", "obs", "connectivities", 
                        "samplem_key", "obs_key", "N") %in% names(data)))
        nam_res <- nam(data, y = y, 
                       batches = batches, 
                       covs = covs, 
                       filter.samples = filter.samples,
                       n.steps = n.steps, 
                       suffix = suffix, 
                       verbose = verbose,
                       ...)
    } else {
        nam_res <- nam.result
    }
    
    ## For association, batches needs to be a numeric vector
    if (is.null(batches)) {
        X <- rep(1, data$N) |> as.factor()
        # X <- Matrix::Matrix(rep(1, data$N), ncol = 1, sparse = TRUE)
    } else {
        X <- dplyr::pull(data$samplem, dplyr::one_of(batches)) |> 
            as.factor()
        
        # batches_df <- data$samplem[, batches, drop = FALSE] 
        # bform <- paste0("~ 0 + ", paste0(colnames(batches_df), collapse = " + ")) |> 
        #     as.formula()
        # X <- Matrix::sparse.model.matrix(bform, batches_df)
    }
    
    if (verbose) message('Perform association testing')
    ## TODO: add filter_samples to nam results
    #         y[nam_res[[paste0('_filter_samples', suffix)]]],
    #         batches[nam_res[paste0('_filter_samples', suffix)]] 
    res <- innerAssociation(NAMsvd = nam_res,
                            y = nam_res$y, 
                            batches_vec = X, 
                            ks = Ks, 
                            Nnull = N.nulls, 
                            force_permute_all = force.permute.all, 
                            allow_duplicate_perms = allow.duplicate.perms, 
                            local_test = local.test, 
                            use.logp = use.logp, 
                            seed = seed)
    if (return.nam) {
        res[[paste0('NAM_embeddings', suffix)]] <- nam_res[[paste0("NAM_nbhdXpc", suffix)]]
        res[[paste0('NAM_loadings', suffix)]] <- nam_res[[paste0("NAM_sampleXpc", suffix)]]
        res[[paste0('NAM_svs', suffix)]] <- nam_res[[paste0("NAM_svs", suffix)]]
    }
    # TODO: add info about kept cells
    #     vars(res)['kept'] = du['keptcells'+suffix]
    
    return(res)
}

