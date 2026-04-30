make_rcna_data <- function(n_samples = 20, n_cells = 200, seed = 42, y.sd = 1) {
    set.seed(seed)
    sample_ids <- paste0("S", seq_len(n_samples))
    cell_ids   <- paste0("C", seq_len(n_cells))
    
    b <- rep(c("A", "B"), each = floor(n_samples / 2))
    d <- rep(paste0("D", 1:2), each = floor(n_samples / 2))
    
    if(n_samples %% 2 == 1){
        b <- factor(c(b, "B"))
        d <- factor(c(d, "D1"))
    }
    
    samplem <- data.frame(
        SampleID  = sample_ids,
        y         = rnorm(n_samples, sd = y.sd),
        batch     = b,
        donor     = d,
        covariate = rnorm(n_samples),
        stringsAsFactors = FALSE
    )
    
    obs <- data.frame(
        CellID   = cell_ids,
        SampleID = sample(sample_ids, n_cells, replace = TRUE),
        stringsAsFactors = FALSE
    )
    
    # sparse symmetric adjacency (connectivities)
    adj <- Matrix::rsparsematrix(n_cells, n_cells, density = 0.05, symmetric = TRUE)
    adj@x <- abs(adj@x)  # non-negative weights
    diag(adj) <- 0
    
    res <- list(
        samplem        = samplem,
        obs            = obs,
        connectivities = adj,
        samplem_key    = "SampleID",
        obs_key        = "CellID",
        N              = n_samples
    )
    return(res)
}