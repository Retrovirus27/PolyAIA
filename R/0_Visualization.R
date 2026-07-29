#' Custom ggprism Theme - Deoxygenated Style
#'
#' A customized ggplot2 theme based on ggprism with grey facet strips,
#' bold text, and enhanced borders for publication-quality figures.
#'
#' @details
#' Theme specifications:
#' - Base size: 11pt
#' - Base line size: 0.25
#' - Strip background: grey95 with black border
#' - Strip text: bold, black, size 11
#' - Panel spacing: 0.5 lines
#' - Panel border: black, 0.5 linewidth
#' - Legend title: size 11
#' - Legend text: size 12
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom scuttle pooledSizeFactors
#'
#' @examples
#' \dontrun{
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   facet_wrap(~cyl) +
#'   th
#' }
#'
#'
#' @export th
th <- ggprism::theme_prism(
  base_size = 5,
  base_line_size = 0.25
) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "black", linewidth = 0.5),
    strip.text = ggplot2::element_text(colour = "black", size = 7),
    panel.spacing = ggplot2::unit(0.5, "lines"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )
th$legend.title <- ggplot2::element_text(size = 8)
th$legend.text  <- ggplot2::element_text(size = 7)

