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
saveRDS(p, "inst/ext/demo_result.rds")




# 
# # # now let's see if rcna recapitulates p from connectivities
#
devtools::document(".")

rcna_demo <- readRDS("inst/ext/demo_data.rds")


data = rcna_demo
data$obs$id <- as.factor(data$obs$id)

res <- list()
batch_vec <- NULL
group_vec <- NULL
cov_mat <- NULL
n.steps = NULL
min.steps = 3
max.steps = 15L
kurt.delta = 3
verbose = FALSE

# c("y", "raw.NAM.T", "qc.NAM.T", "resid.NAM.T", "_M", "_r", "NAM_sampleXpc", 
#   "NAM_svs", "NAM_varexp", "NAM_nbhdXpc", "keptcells")

NAM <- buildNAM(rcna_demo)

nam(data = rcna_demo, y = "case")
association(rcna_demo, y = "case")


