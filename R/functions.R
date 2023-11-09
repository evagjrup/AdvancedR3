#' Descriptive stats of each metabolite
#'
#' @param data Lipidomics df
#'
#' @return A data.frame/tibble
descriptive_stats <- function(data) {
    data %>%
        dplyr::group_by(metabolite) %>%
        dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, digits = 1)))
}

#' Plot distributions
#'
#' @param data lipidomics data
#'
#' @return metabolite distribution plot
plot_distributions <- function(data) {
    metabolite_distribution_plot <- ggplot2::ggplot(data,
                                                    ggplot2::aes(x = value)) +
        ggplot2::geom_histogram() +
        ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
    return(metabolite_distribution_plot)
}
