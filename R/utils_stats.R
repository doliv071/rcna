#' Permute Y conditioned on B
#' 
#' @param B a numeric batch vector
#' @param Y the variable of interest, scaled
#' @param num Number of permutations to perform
#' 
#' @returns A matrix of permutations of Y conditioned on B
#' 
#' @keywords internal
conditional_permutation <- function(B, Y, num, duplicates.ok = TRUE) {
    # TODO: for certain combinations of length(Y) and B it is better to compute
    #       all permutations rather than sampling and filtering duplicates.
    y_perm <- lapply(seq_len(num), function(i) {
        res <- split(seq_len(length(Y)), B) |> 
            lapply(\(idx) {
                data.frame(idx, val = sample(Y[idx]))
            }) |> do.call(what = rbind)
        res[order(res$idx), "val", drop = TRUE]
    })
    y_perm <- matrix(unlist(y_perm), nrow = length(y_perm), byrow = TRUE)
    dup_y <- duplicated(y_perm)
    if(sum(dup_y) > 0.2*num){
        warning("Found ", sum(dup_y), " duplicated permutations.", immediate. = TRUE)
    }
    if(!duplicates.ok){
        y_perm <- y_perm[!dup_y,]
    }
    return(t(y_perm))
}

#' calculate the number of correlations within each threshold
#'
#' @param breaks A vector of thresholds 
#' @param z A matrix of correlations 
#' 
#' @returns a vector of length `ncol(z)` containing the the number of neighborhoods 
#' exceeding each break (threshold)
#' 
#' @keywords internal
tail_counts <- function(breaks, z) {
    tailsize <- nrow(z)
    breaks <- breaks^2
    z <- z^2
    res <- apply(z, 2, function(zi) {
        binsums <- findInterval(zi, breaks) |> 
            tabulate(nbins = length(breaks) - 1) |> 
            cumsum()
        tailsize - binsums
    })
    return(res)
}

#' compute the empirical FDRs via FDP procedure
#'
#' @param z A correlation matrix with a single column containing the observed 
#' correlations
#' @param znull A correlation matrix with n-permutations columns containing the 
#' permuted null correlations
#' @param threshold a vector of bins for the distribution
#' 
#' @returns A data.frame containing thresholds, FDRs, and number of neighborhoods in each
#' threshold bin.
#' 
#' @keywords internal
empirical_fdrs <- function(z, znull, thresholds) {    
    n <- length(thresholds) - 1
    tails <- t(tail_counts(thresholds, znull)[1:n, ])
    ranks <- t(tail_counts(thresholds, z)[1:n, ])
    # compute FDPs
    fdp <- sweep(tails, 2, ranks, '/')
    # make sure we didn't over-shoot upper bins too much
    valid_thresh <- setdiff(1:ncol(fdp), unique(which(is.na(fdp), arr.ind = TRUE)[,2]))
    ret_thresh <- thresholds[valid_thresh]
    fdp <- fdp[,valid_thresh]
    fdp[fdp > 1] <- 1
    fdr <- Matrix::colMeans(fdp)
    # enforce monotonicity
    fdr_monotone <- cummin(fdr)
    num_detected <- vapply(ret_thresh, function(.t) sum(abs(z) > .t), numeric(1))
    fdrs <- data.frame(threshold = ret_thresh,
                       fdr = fdr_monotone, 
                       num_detected = num_detected) 
    return(fdrs)
}
