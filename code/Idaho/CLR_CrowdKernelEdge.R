################################################################################
## Estimate the crowding kernel with distances from the plant outline (Idaho)
##
## Written by Hung-wei Lin Nov 2026
##
## Sensitivity version of CLR_CrowdKernel.R: the same neighbour table, kernels, range and
## growth-form check, with every distance measured from the outline of the focal plant
## (0 for neighbours inside it) in place of its centre. With distances from the centre a large plant
## has no neighbours inside its own area, where the kernel is strongest, so its crowding comes out
## low. Files and figures end in '_edge'.
## Then compares the two ways of measuring distance: same plants and regressions, so the AIC of the
## main kernel (reach R_ref, each version's chosen penalty) can be compared directly
## ('CLR_CrowdKernel_origin.csv').

ring_origin <- "edge"
source(here("code", "Idaho", "CLR_CrowdKernel.R"))

# Centre or outline: comparison of the two kernels -----------------------------
fit_o <- lapply(c(centre = "_centre", edge = "_edge"), function(k) readRDS(file.path(der_dir, sprintf("CLR_CrowdKernel_fit%s.rds", k))))
origin_tab <- data.frame(origin = c("centre", "outline"),
                         log10_lambda = sapply(fit_o, function(f) log10(f$lambda[f$best])),
                         aic = sapply(fit_o, function(f) f$aic[f$best]), dof = sapply(fit_o, function(f) f$dof[f$best]),
                         R_main = sapply(fit_o, `[[`, "R_main"), R_forb = sapply(fit_o, `[[`, "R_forb"))
origin_tab$delta_aic <- origin_tab$aic - min(origin_tab$aic)
print(origin_tab, digits = 5, row.names = FALSE)
write.csv(origin_tab, file.path(tab_dir, "CLR_CrowdKernel_origin.csv"), row.names = FALSE)
rm(fit_o, origin_tab)
