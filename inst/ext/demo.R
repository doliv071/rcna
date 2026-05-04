# 
# # # check that we can replicate the workflow from the anndata
#
library(anndataR)
library(SingleCellExperiment)
library(ggplot2)
sce_obj <- read_h5ad("~/Git/cna-rcna-test/test.h5ad", as = "SingleCellExperiment")

reducedDim(sce_obj, "X_umap") |> 
    cbind(colData(sce_obj)[, c("case_coef", "case_coef_fdr")]) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = case_coef)) +
    ggplot2::scale_color_gradient2(low = "darkblue", 
                                   high = "darkred") + 
    theme_classic()

reducedDim(sce_obj, "X_umap") |> 
    cbind(colData(sce_obj)[, c("male_coef", "male_coef_fdr")]) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = male_coef)) +
    ggplot2::scale_color_gradient2(low = "darkblue", 
                                   high = "darkred") + 
    theme_classic()

reducedDim(sce_obj, "X_umap") |> as.data.frame() |> 
    cbind(leiden = as.factor(sce_obj$leiden)) |> 
    ggplot() + 
    geom_point(aes(x = V1, y = -V2, color = leiden)) +
    theme_classic()

# Not easy to do gradient fill for violins
reducedDim(sce_obj, "X_umap") |> as.data.frame() |> 
    cbind(leiden = as.factor(sce_obj$leiden), 
          male_coef = sce_obj$male_coef) |> 
    ggplot(aes(x = leiden, y = male_coef, color = leiden)) + 
    geom_violin() +
    theme_classic()

# 
# # # get the association object including nam and connectivities and save to an rds 
#
library(reticulate)
library(HSSRscripts)
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test.py")


connectivities <- d$obsp["connectivities"]
connectivities <- as(connectivities, "CsparseMatrix")

obs <- d$obs |> rownames_to_column("cell_id") |> 
    dplyr::mutate(id = as.factor(id), case = as.factor(case), 
                  male = as.factor(male), batch = as.factor(batch))

samplem <- obs |> dplyr::distinct(id, .keep_all = T) |> 
    dplyr::mutate(cell_id = NULL, case_coef = NULL, case_coef_fdr = NULL)

rownames(samplem) <- paste0("sample_", samplem$id)

demo <- list(
    connectivities = connectivities,
    samplem = samplem,
    samplem_key = "id",
    obs = obs,
    obs_key = "cell_id",
    N = nrow(samplem)
)

saveRDS(demo, "inst/ext/demo_data.rds")

p <- builtins$vars(p)
saveRDS(p, "inst/ext/demo_result_nocov_nobatch.rds")




# 
# # # now let's see if rcna recapitulates p from connectivities
#
devtools::document(".")

rcna_demo <- readRDS("inst/ext/demo_data.rds")

a <- nam(data = rcna_demo, y = "case")
b <- association(rcna_demo, nam.result = a, y = "case")

p <- readRDS("inst/ext/demo_result_nocov_nobatch.rds")

# only attributes differ
all.equal(as(as.matrix(p$M), "CsparseMatrix"), 
          as(as(a$`_M`, "CsparseMatrix"), "dsCMatrix"))

p$namresid |> as.matrix() |> as("sparseMatrix") |> dim()
Matrix::t(a$resid.NAM.T) |> dim()

# 50 x 10000 sparse Matrix of class "dgCMatrix"
all.equal(as(as.matrix(p$namresid), "sparseMatrix"), 
          # 50 x 10000 sparse Matrix of class "dgCMatrix"
          Matrix::t(a$resid.NAM.T))


#################################
# because things seem to diverge quickly, going through step by step :/
##############
library(reticulate)
library(HSSRscripts)
devtools::document(".")
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test.py")
rcna_demo <- readRDS("inst/ext/demo_data.rds")

#
### checking diffusion step...
#
diffuse_step1 <- d1 |> as.matrix() |> as("CsparseMatrix")


f <- as.formula(paste0("~ 0 + ", rcna_demo$samplem_key))
s <- Matrix::sparse.model.matrix(f, rcna_demo$obs)
colnames(s) <- gsub(rcna_demo$samplem_key, '\\1', colnames(s))
rownames(s) <- rcna_demo$obs[[rcna_demo$obs_key]]
diffuseStep <- function(a, s) {
    stopifnot((isa(a, "Matrix") || is.matrix(a)) && nrow(a) == ncol(a))
    stopifnot(isa(s, "Matrix") || is.matrix(s))
    stopifnot(nrow(s) == ncol(a))
    
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(res) 
}

rdiffuse_s1 <- s <- diffuseStep(rcna_demo$connectivities, s)

all.equal(rdiffuse_s1, diffuse_step1)

rdiffuse_s2 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step2 <- d2 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s2, diffuse_step2)

s <- diffuseStep(rcna_demo$connectivities, s)
rdiffuse_s4 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step4 <- d4 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s4, diffuse_step4)

s <- diffuseStep(rcna_demo$connectivities, s)
s <- diffuseStep(rcna_demo$connectivities, s)
s <- diffuseStep(rcna_demo$connectivities, s)
rdiffuse_s8 <- s <- diffuseStep(rcna_demo$connectivities, s)

diffuse_step8 <- d8 |> as.matrix() |> as("CsparseMatrix")

all.equal(rdiffuse_s8, diffuse_step8)

#### diffusion steps match exactly.

#
### checking nam step...both took 4 steps
#
pysnorm <- as.matrix(snorm) |> as("CsparseMatrix")

rsnorm <- buildNAM(rcna_demo)

all.equal(pysnorm, rsnorm)
# [1] TRUE

#### nam steps match exactly.

#
### checking qcnam step...
#

# NO BATCH
pyqc_nam <- qc_nam |> as.matrix() |> as("CsparseMatrix")
rqcnam <- buildNAM(rcna_demo) |> qcNAM()

all.equal(pyqc_nam, rqcnam$NAM)
# TRUE

# qc works for both batch and batch+covs

#
### checking svdnam step...
#

svd_nam[[1]]
rnam <- buildNAM(rcna_demo)
resid_rnam <- residNAM(rnam)
rsvd_nam <- svdNAM(resid_rnam$NAM)
all.equal(abs(as.matrix(svd_nam[[1]])), abs(as.matrix(rsvd_nam[[1]])))
# TRUE (after much painful bug finding)

#### svd steps match


#########################
# NAM processing matches exactly
###

#
### checking association step...
#
rnam <- nam(rcna_demo, y = "case")


r1 <- sweep(rnam$resid.NAM.T, 2, scale(rnam$y), "*") |>  Matrix::colMeans()
r2 <- Matrix::colMeans(rnam$resid.NAM.T) * scale(rnam$y)[,1]
r3 <- Matrix::colMeans(scale(rnam$y) * Matrix::t(rnam$resid.NAM.T)) |> matrix(ncol = 1)
all.equal(r1, r2)
V <- rnam$NAM_nbhdXpc
sv <- rnam$NAM_svs
ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta / n)

rass <- association(rcna_demo, y = "case")

# global p is identical
identical(rass$p, p$p)
# TRUE

identical(rass$k, p$k)
# TRUE

all.equal(as.numeric(rass$yhat), as.numeric(p$yresid_hat))
# TRUE

all.equal(as.numeric(rass$ycond), as.numeric(p$yresid))
# TRUE

all.equal(as.numeric(rass$ncorrs), as.numeric(p$ncorrs))
# "Mean relative difference: 0.001598857"

all.equal(as.numeric(rass$r2_perpc), as.numeric(p$r2_perpc))
# TRUE

all.equal(p$fdr_5p_t, rass$fdr_5p_t)

p$fdrs$threshold

data.frame(py.thresh = p$fdrs$threshold, r.thresh = rass$fdrs$threshold[100:399])

plot(p$fdrs$threshold, rass$fdrs$threshold[100:399])
abline(0,1)

plot(p$fdrs$fdr, rass$fdrs$fdr[100:399])
abline(0,1)

plot(p$fdrs$num_detected, rass$fdrs$num_detected[100:399])
abline(0,1)

# these are always random
all.equal(rass$minps_null, as.numeric(p$nullminps))
plot(rass$minps_null, as.numeric(p$nullminps))


test.iters <- parallel::mclapply(1:100, \(i) {
    suppressWarnings(
        association(nam.result = rnam, y = "case", verbose = FALSE)$fdrs$fdr
    )
}, mc.cores = 12)

test.mat <- do.call(cbind, test.iters)

all.equal(test.mat[,1], test.mat[,2])

cor.p <- apply(test.mat[100:399,], 2, \(x){ 
    data.frame(cor = cor(p$fdrs$fdr, x), 
               mae = mean(abs(p$fdrs$fdr - x )))
}) |> do.call(what = rbind)

plot(cor.p)


##################################################
# One more time with batch and covariates 
######
library(reticulate)
library(HSSRscripts)
devtools::document(".")
use_virtualenv("~/Git/cna-rcna-test/cna_env/", required = TRUE)
builtins <- import_builtins()
# py_run_string("import vars")
source_python("~/Git/cna-rcna-test/test.py")
rcna_demo <- readRDS("inst/ext/demo_data.rds")


rcna_demo$samplem$male <- as.numeric(rcna_demo$samplem$male)
rass2 <- association(rcna_demo, y = "case", batches = "batch")
rnam2 <- nam(data = rcna_demo, y = "case", batches = "batch")

qc.NAM <- Matrix::t(rnam2$qc.NAM.T)
qc_nam <- qc_nam |> as.matrix() |> as("CsparseMatrix")
all.equal(qc.NAM, qc_nam)

prn <- p2$namresid |> as.matrix() |> as("Matrix")
rrn <- Matrix::t(rnam2$resid.NAM.T)
all.equal(prn, rrn)
plot(prn[,1], rrn[,1])


plot(rass2$fdrs$fdr[100:399], p2$fdrs$fdr)
abline(0,1)

plot(rass2$fdrs$threshold[100:399], p2$fdrs$threshold)
abline(0,1)








