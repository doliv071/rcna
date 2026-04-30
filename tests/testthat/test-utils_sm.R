test_that("sparseMatrix function: Kurtosis is working as expected", {

    x <- Matrix::rsparsematrix(100, 200, density = 1, symmetric = FALSE,
                               rand.x = function(n) signif(rnorm(n, sd = 3), 2))
    K <- Kurtosis(x)
    K_col <- Kurtosis(x, margin = 2)
    K_row <- Kurtosis(x, margin = 1)
    
    k <- moments::kurtosis(as.vector(x))
    k_col <- moments::kurtosis(as.matrix(x))
    k_row <- moments::kurtosis(as.matrix(Matrix::t(x)))
    
    expect_true(all.equal(k, K))
    expect_true(all.equal(k_col, K_col))
    expect_true(all.equal(k_row, K_row))
    expect_error(Kurtosis(x, margin = 3))
    
})

test_that("sparseMatrix function: propTable is working as expected", {
    
    x <- Matrix::rsparsematrix(10, 20, density = 1, symmetric = FALSE,
                               rand.x = function(n) signif(runif(n, min = 0, max = 10), 2))
    pt2 <- proportions(as.matrix(x), 2) 
    pT2 <- propTable(x, margin = 2) |> as.matrix()
    
    pt1 <- proportions(as.matrix(x), 1) 
    pT1 <- propTable(x, margin = 1) |> as.matrix()
    all.equal(pt1, pT1)
    all.equal(pt2, pT2)
    
})

test_that("sparseMatrix function: propTable is working as expected", {
    
    x <- Matrix::rsparsematrix(10, 20, density = 1, symmetric = FALSE,
                               rand.x = function(n) signif(runif(n, min = 0, max = 10), 2))
    scs <- scale(as.matrix(x), center = TRUE, scale = TRUE)
    attributes(scs) <- list(dim = attributes(scs)$dim)
    sc <- scale(as.matrix(x), center = TRUE, scale = FALSE)
    attributes(sc) <- list(dim = attributes(sc)$dim)
    ss <- scale(as.matrix(x), center = FALSE, scale = TRUE)
    attributes(ss) <- list(dim = attributes(ss)$dim)
    s <- scale(as.matrix(x), center = FALSE, scale = FALSE)
    attributes(s) <- list(dim = attributes(s)$dim)
    
    Scs <- Scale(x, center = TRUE, scale = TRUE) |> as.matrix()
    Sc <- Scale(x, center = TRUE, scale = FALSE) |> as.matrix()
    Ss <- Scale(x, center = FALSE, scale = TRUE) |> as.matrix()
    S <- Scale(x, center = FALSE, scale = FALSE) |> as.matrix()
    
    all.equal(scs, Scs)
    all.equal(sc, Sc)
    all.equal(ss, Ss)
    all.equal(s, S)
    
})








