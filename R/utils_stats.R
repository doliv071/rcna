#' Permute Y conditioned on B
#' 
#' @param B A factor (or object coercible via \code{split()}) of length n
#' indicating batch membership. Permutations are performed within each level.
#' @param Y A numeric vector of length n. The variable to permute.
#' @param num Positive integer. Number of permutations to generate.
#' @param seed An integer random seed for reproducibility.
#' @param duplicates.ok Logical. If FALSE, duplicate permutation rows are
#' removed from the output, which may result in fewer than `num` columns.\cr
#' Default: TRUE.
#' 
#' @returns A numeric matrix of dimensions n × m, where `m <= num`. 
#' Each column is one permutation of `Y` within batch levels.
#'   
#' @keywords internal
conditional_permutation <- function(B, Y, num, seed = NULL, duplicates.ok = TRUE) {
    stopifnot(is.numeric(Y))
    stopifnot(length(B) == length(Y))
    stopifnot(is.numeric(num) && length(num) == 1 && num >= 1)
    stopifnot(is.logical(duplicates.ok) && length(duplicates.ok) == 1)
    if(is.null(seed)){
        seed <- sample(1e6, 1)
    }
    set.seed(seed)
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
#' @returns A numeric matrix of dimensions (\code{length(breaks) - 1} x \code{ncol(z)}), 
#' where entry [i, j] is the number of rows in column j of \code{z} whose squared 
#' value exceeds \code{breaks[i]^2}.
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
#' @param znull A correlation matrix with p-permutations columns containing the 
#' permuted null correlations
#' @param threshold a vector of bins for the distribution
#' 
#' @note The raw per-permutation FDP estimates are averaged across permutations
#' to obtain the FDR, which is then monotonised via `cummin()` to ensure
#' that the FDR curve is non-increasing as the threshold increases.
#' 
#' @returns A data.frame containing thresholds, FDRs, and number of neighborhoods 
#' in each threshold bin.
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
