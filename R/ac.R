#' Autocorrelation and partial autocorrelation function.
#'
#' @import tibble dplyr ggplot2 patchwork
#' @param data A numeric, time series, or xts variable.
#' @param max_lag A number that represents the maximum lag order for the ACF and PACF.
#' @param main_title Optional plot title.
#' @param output Return data or not.
#' @return The respective ACF and PACF functions.
#' @examples
#' ac(rnorm(100))


ac  <- function(data, max_lag = NULL, main_title = NULL, output = NULL) {

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

  # confidence intervals
  num <- length(data)
  U <- 2 / sqrt(num)
  L <- -U

  # conditions
  if (num > 49 & is.null(max_lag))
    max_lag = ceiling(10 + sqrt(num))
  if (num < 50 & is.null(max_lag))
    max_lag = floor(5 * log10(num))
  if (max_lag > (num - 1))
    stop("Number of lags exceeds number of observations")

  # calculations
  cf <- tibble::tibble(
    ACF = acf(data, max_lag, plot = FALSE)$acf[-1],
    PACF = as.numeric(pacf(data, max_lag, plot = FALSE)$acf),
    LAG = 1:max_lag
  ) |>
    dplyr::select(LAG, ACF, PACF)

  # plot
  p1 <- cf |>
    ggplot2::ggplot() +
    geom_bar(aes(x = LAG, y = ACF),
             stat = "identity",
             fill = tsm_pal$default[[2]]) +
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
    labs(y = "ACF", x = NULL) +
    geom_hline(yintercept = c(U, L), linetype = "dashed",
               color = tsm_pal$default[[3]], linewidth = 1)


  p2 <- cf |>
    ggplot2::ggplot() +
    geom_bar(aes(x = LAG, y = PACF),
             stat = "identity",
             fill = tsm_pal$default[[2]]) +
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
    labs(y = "PACF", x = NULL) +
    geom_hline(yintercept = c(U, L), linetype = "dashed",
               color = tsm_pal$default[[3]], linewidth = 1)

  cf_plot <- p1 + p2 +
    plot_layout(nrow = 1, byrow = FALSE) +
    plot_annotation(title = main_title)

  if (is.null(output)) {
    return(cf_plot)
  } else {
    return(list(cf_plot, cf))
  }

}

#ac(rnorm(100))
