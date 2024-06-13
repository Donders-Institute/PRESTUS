library(data.table)
library(Matrix)
data_path <- '/project/3015999.02/andche_sandbox/TUS_sims/tusim/data/images_Sjoerd/'
data <- fread(paste0(data_path, "results_aggregated.csv"))

data[,rel_peak_pressure:=(.SD[medium=='water',max_pressure]-max_pressure)/.SD[medium=='water', max_pressure], by = .(subject_id, target, layer)]
eucl_norm<- function(ref_xyz, targ_xyz){
  ref_xyz <- Matrix(rep(as.matrix(ref_xyz), each = nrow(targ_xyz)), ncol = 3)
  targ_xyz <- Matrix(targ_xyz)
  sqrt(rowSums((ref_xyz-targ_xyz)^2))

}
data[,shift_in_focal_pos_mm:=eucl_norm(.SD[medium=='water',.(max_pos_x, max_pos_y, max_pos_z)], cbind(max_pos_x, max_pos_y, max_pos_z))/0.5, by = .(subject_id, target,layer)]

fwrite(data, paste0(data_path, "results_aggregated.csv"))
