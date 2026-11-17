################################################################################
## Version of the Idaho pipeline (read by both runscripts and the ARC task script)
##
## Written by Hung-wei Lin Nov 2026
##
## Four versions: climate predictor x origin of the crowding distances. Steps shared by every version
## (checks, tracking, Part A data, both crowding kernels) write to the 'Idaho' folders; the rest
## (crowding indices, Part A input, models) to the folders of the version, e.g. 'Idaho(1Env&centreK)'.

  # "1Env": one climate variable; "PC": principal components of the climate variables
clim_version <- "1Env"
  # "centre": crowding distances from the plant centre (Teller et al. 2016); "edge": from its outline
kernel_origin <- "centre"

ver_tag <- sprintf("Idaho(%s&%sK)", clim_version, kernel_origin)
