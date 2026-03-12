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
#' calculated
#' 
#' @returns A list containing p-value and r^2 for the F-test
#' 
#' @keywords internal
calcStats <- function(yhat, ycond, k, n, r, log.p = TRUE) {
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
#' 
#' @returns A list containing the k, p, and r^2 for k which produces the minimum
#' p-value.
#' 
#' @keywords internal
minpStats <- function(M, y, ks, U, r) {
    n <- length(y)
    zcond <- M %*% y
    zcond <- scale(zcond, center = FALSE, scale = TRUE)
    qhats <- lapply(ks, function(k) regress(zcond, k, U)$qhat)
    .tmp <- lapply(seq_along(ks), function(i) calcStats(qhats[[i]], zcond, ks[i], n, r))
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
                             seed = NULL) {
    if (is.null(seed)) {
        seed <- sample(1e6, 1)
        set.seed(seed)
    }
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
    mp <- minpStats(M, y, ks, U, r)
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
                                      duplicates.ok = allow_duplicate_perms)
    .tmp <- apply(y_null, 2, \(z) minpStats(M, z, ks, U, r))
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
#' @param y A vector with contrast variable value to be tested for association.
#' @param batches A character string to denote batch variables.\cr
#' Default: NULL
#' @param covs A character string or vector to denote covariate variables.\cr
#' Default: NULL
#' @param nsteps TBD.\cr
#' Default: NULL
#' @param suffix A character string to be appended for unknown reasons. N.B. 
#' There is a bug that will cause an error if this is combined with 
#' `return_nam = TRUE`.\cr
#' Default: ''
#' @param return_nam Logical controlling whether or not to return the NAM. See
#' above for known bug.\cr
#' Default: FALSE
#' @param verbose Logical controlling the verbosity of the function.\cr
#' Default: FALSE
#' @param force_recompute STUB. Ignored
#' 
#' @return A list containing p, nullminps, k, ncorrs, fdrs, fdr_5p_t, fdr_10p_t, 
#' yhat, ycond, ks, beta, r2, r2_perpc, nullr2_mean, and nullr2_std. If 
#' `return_nam = TRUE`, then additionally NAM_embeddings, NAM_loadings, and NAM_svs
#' 
#' @export 
association <- function(data, 
                        y, # TODO: y should be taken from data$samplem directly 
                        batches = NULL, 
                        covs = NULL, 
                        nsteps = NULL, 
                        suffix = '',
                        force_recompute = FALSE, 
                        return_nam = FALSE, 
                        verbose = TRUE, 
                        ...) {
    
    # formatting and error checking
    
    
    ## CHECK: need _df_to_array in R? 
    #     covs = _df_to_array(data, covs)
    #     batches = _df_to_array(data, batches)
    #     y = _df_to_array(data, y)
    
    ## TODO: check lengths
    #     if y.shape != (data.N,):
    #         raise ValueError(
    #             'y should be an array of length data.N; instead its shape is: '+str(y.shape))
    
    ## TODO: add sample filtering 
    #     if covs is not None:
    #         filter_samples = ~(np.isnan(y) | np.any(np.isnan(covs), axis=1))
    #     else:
    #         filter_samples = ~np.isnan(y)
    
    ## Here, data has all the du things 
    #     du = data.uns
    if (verbose) message('Build NAM PCs')
    nam_res <- nam(data, 
                   batches = batches, 
                   covs = covs, 
                   filter.samples = filter_samples,
                   nsteps = nsteps, 
                   suffix = suffix,
                   force.recompute = force_recompute)
    
    ## For association, batches needs to be a numeric vector
    if (is.null(batches)) {
        batches_vec <- rep(1, data$N)
    } else {
        batches_vec <- dplyr::select(data$samplem, dplyr::one_of(batches)) |>
            as.matrix() |>
            as.integer()
    }
    
    if (verbose) message('Perform association testing')
    ## TODO: add filter_samples to nam results
    #         y[nam_res[[paste0('_filter_samples', suffix)]]],
    #         batches[nam_res[paste0('_filter_samples', suffix)]] 
    res <- innerAssociation(NAMsvd = nam_res,
                            y = y, 
                            batches_vec = batches_vec, 
                            ...)
    
    if (return_nam) {
        #         res[[paste0('NAM_embeddings', suffix)]] <- nam_res$NAM_nbhdXpc
        #         res[[paste0('NAM_loadings', suffix)]] <- nam_res$NAM_sampleXpc
        #         res[[paste0('NAM_svs', suffix)]] <- nam_res$NAM_svs
        res[['NAM_embeddings']] <- nam_res$NAM_nbhdXpc
        res[['NAM_loadings']] <- nam_res$NAM_sampleXpc
        res[['NAM_svs']] <- nam_res$NAM_svs
    }
    # TODO: add info about kept cells
    #     vars(res)['kept'] = du['keptcells'+suffix]
    
    return(res)
}

