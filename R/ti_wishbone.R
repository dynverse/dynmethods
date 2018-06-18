#' @importFrom viridis scale_colour_viridis
plot_wishbone <- function(prediction) {
  g <- ggplot() +
    geom_point(aes(comp_1, comp_2, color = time), prediction$dimred %>% mutate_at(c("comp_1", "comp_2"), dynutils::scale_minmax)) +
    viridis::scale_colour_viridis() +
    labs(colour = "Trajectory") +
    theme(legend.position = c(.92, .12))
}
