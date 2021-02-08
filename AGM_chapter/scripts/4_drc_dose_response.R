# dose resposne effects

library(tidyverse)
library(MPStats)

cat("Fitting score by compound dose\n")

load("intermediate_data/well_scores.Rdata")
load("intermediate_data/compound_moa.Rdata")

score_by_dose_fits <- MPStats::fit_drc_score_by_dose(well_scores=well_scores)
save(score_by_dose_fits, file="intermediate_data/score_by_dose_fits.Rdata")
score_by_dose_fits %>% readr::write_tsv(
  path=paste0("product/score_by_dose_fits_", MPStats::date_code(), ".tsv"))

plot_all <- MPStats::plot_drc_score_by_dose(
  well_scores=well_scores,
  fits=score_by_dose_fits,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot_all,
  filename=paste0("product/score_by_dose_", MPStats::date_code(), ".pdf"),
  height=20,
  width=20)
ggplot2::ggsave(
  plot=plot_all,
  filename=paste0("product/score_by_dose_", MPStats::date_code(), ".png"),
  height=20,
  width=20)

compound_moa %>%
  dplyr::mutate(moa2 = ifelse(compound == "epothilone B", "cytoskeleton", moa2)) %>%
  plyr::d_ply("moa2", function(compounds){
    moa_type = compounds$moa2[1]
    cat("Plotting '", moa_type, "' compounds ...\n", sep="")
    plot <- MPStats::plot_drc_score_by_dose(
      well_scores = well_scores %>%
        dplyr::semi_join(compounds, by="compound"),
      fits = score_by_dose_fits %>%
        dplyr::semi_join(compounds, by="compound"),
      subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
      
    
    # use the number of rows x columns to figure out plot size
    facet_dims <- compounds %>% nrow() %>% ggplot2::wrap_dims()
    height <- facet_dims[1] * 2
    width <- facet_dims[2] * 2
    
    ggplot2::ggsave(
      plot=plot,
      filename=paste0("product/score_by_dose_", moa_type, "_", MPStats::date_code(), ".pdf"),
      height=height,
      width=width)
    ggplot2::ggsave(
      plot=plot,
      filename=paste0("product/score_by_dose_", moa_type, "_", MPStats::date_code(), ".png"),
      height=height,
      width=width)
  })


