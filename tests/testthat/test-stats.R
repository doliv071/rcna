test_that("stats working as expected", {
    y <- rnorm(16)
    b <- sample(letters[1:4], 16, replace = TRUE) |> factor()
    g <- sample(LETTERS[1:2], 16, replace = TRUE) |> factor()
    nnull <- 10
    foo <- conditional_permutation(Y = y, B = b, num = nnull)
    expect_equal(dim(foo), c(length(y), nnull))
    
    NAMres <- nam(data = make_rcna_data(y.sd = 3), y = "y", batches = NULL, covs = NULL)
    
    ### This is pulle from innerNAM so we can work with pretend real data? should we?
    U <- NAMres[["NAM_sampleXpc"]]
    sv <- NAMres[["NAM_svs"]]
    V <- NAMres[["NAM_nbhdXpc"]]
    M <- NAMres[["_M"]] 
    r <- NAMres[["_r"]]
    y <- scale(NAMres[["y"]])
    n <- length(y)
    batches_vec = NAMres[["_batches"]]
    group_vec = NAMres[["_donors"]]
    
    ks <- seq(n/50, n/5, by = 1) |> ceiling() |> unique()
    k <- minpStats(M, y, ks, U, r)$k
    
    ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
    beta <- regress(ycond, k, U)$beta
    ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta / n)
    rownames(ncorrs) <- rownames(V)
    
    # compute final p-value using Nnull null f-test p-values
    y_null <- conditional_permutation(Y = y, 
                                      B = batches_vec, 
                                      G = group_vec,
                                      num = nnull, 
                                      duplicates.ok = TRUE)

    Nnull <- min(1000, ncol(y_null))
    y_null <- y_null[, 1:Nnull]
    ycond_null <- scale(M %*% y_null, center = FALSE, scale = TRUE)
    gamma_null <- crossprod(U[, 1:k], ycond_null)  
    ncorrs_null <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_null / n)))
    
    maxcorr <- max(c(abs(ncorrs), abs(ncorrs_null))) + sqrt(.Machine$double.eps)
    fdr_thresholds <- seq(0, maxcorr, maxcorr/400) 
    
    nthresh <- length(fdr_thresholds) - 1
    tails <- tail_counts(breaks = fdr_thresholds, z = ncorrs_null)[1:nthresh, ] |> t()
    expect_equal(dim(tails), c(nnull, nthresh))
    
    ranks <- tail_counts(breaks = fdr_thresholds, z = ncorrs)[1:nthresh, ] |> t()
    expect_equal(dim(ranks), c(1, nthresh))
    
    fdp <- sweep(tails, 2, ranks, '/')
    
    valid_thresh <- setdiff(1:ncol(fdp), unique(which(is.na(fdp), arr.ind = TRUE)[,2]))
    ret_thresh <- fdr_thresholds[valid_thresh]
    fdp <- fdp[,valid_thresh]
    
    fdp[fdp > 1] <- 1
    fdr <- Matrix::colMeans(fdp)
    
    fdr_monotone <- cummin(fdr)
    fdrs <- empirical_fdrs(z = ncorrs, znull = ncorrs_null, thresholds = fdr_thresholds)
    
})

