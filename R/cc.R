#' Cross correlation function that considers the relationship between two variables.
#'
#' @import tibble dplyr ggplot2 patchwork
#' @param x A numeric, time series, or xts variable.
#' @param y A numeric, time series, or xts variable.
#' @param max_lag A number that represents the maximum lag order for the ACF and PACF.
#' @param main_title Optional plot title.
#' @param output Return data or not.
#' @return The respective ACF and PACF functions.
#' @examples
#' cc(rnorm(100), rnorm(100))


cc  <- function(x, y, max_lag = NULL, main_title = NULL, output = NULL) {

  # colours
  tsm_pal <- list(
    default = tribble(~ colour, ~ hex,
                      'blue', '#4682B4',
                      'purple', '#B44683',
                      'gold', '#B4A446',
                      'darkblue','#31697E',
                      'darkpurple','#7E315C'
    ) |> deframe()
  )

  # number of observations
  num <- length(x)

  # conditions
  if (num > 49 & is.null(max_lag))
    max_lag = ceiling(10 + sqrt(num))
  if (num < 50 & is.null(max_lag))
    max_lag = floor(5 * log10(num))
  if (max_lag > (num - 1))
    stop("Number of lags exceeds number of observations")

  # calculations
  tmp <- ccf(x, y, max_lag, plot = FALSE)

  # confidence intervals - Enders 2014
  cf <- tibble::tibble(
    LAG = tmp$lag[,1,1],
    CCF = tmp$acf[,1,1],
  ) %>%
    dplyr::filter(LAG >= 0) %>%
    mutate(
      U = 2 * ((num - 0:max_lag)^(-1/2)),
      L = -U
    )

  # plot
  p1 <- cf %>%
    ggplot2::ggplot() +
    geom_bar(aes(x = LAG, y = CCF),
             stat = "identity",
             fill = tsm_pal$default[[2]]) +
    geom_line(aes(x = LAG, y = U),
              stat = "identity",
              colour = tsm_pal$default[[3]],
              linetype = "dashed",
              linewidth = 1) +
    geom_line(aes(x = LAG, y = L),
              stat = "identity",
              colour = tsm_pal$default[[3]],
              linetype = "dashed",
              linewidth = 1) +
    theme_light() +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0.6, 0.6, 0.6, 0.6, "cm"),
      legend.justification = c(0,1),
      legend.position = 'none',
      legend.title = element_blank()
    ) +
    expand_limits(y = c(-0.91, 0.91)) +
    labs(y = "CCF", x = NULL)

  cf_plot <- p1 +
    plot_layout(nrow = 1, byrow = FALSE) +
    plot_annotation(title = main_title)

  if (is.null(output)) {
    return(cf_plot)
  } else {
    return(list(cf_plot, cf))
  }

}
