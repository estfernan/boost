# plot_standard_standard_st
# plot_binary_st

##
## R package boost by Esteban Fernández, Xi Jiang, Suhana Bedi, and Qiwei Li
## Copyright (C) 2021
##
## This file is part of the R package boost.
##
## The R package boost is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package boost is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Spatial Transcriptomics Plot
##'
##' Plot the absolute, dichotomised, or relative gene expression levels in a
##' spatial transcriptomics (ST) dataset.
##'
##' @param expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   absolute, dichotomised, or relative gene expression levels. Each entry is
##'   an integer value that denotes the gene count at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param gene.name An optional character string that specifies the gene name.
##'   The default is `NULL` which appends nothing to the plot title.
##' @param binary A logical value that indicates if the gene expression levels
##'   are dichotomised. The default is `FALSE` to expect absolute or relative
##'   expression levels.
##' @param log.expr A logical value that indicates if the gene expression levels
##'   should be plotted at log-level. The default is `TRUE` to log-transform
##'   the gene expression levels if they are not dichotomised.
##'
##' @return A `ggplot` object.
##'
##' @examples
##' ## Need to implement this example.
##'
##' @export
##' @keywords plot
##'
st.plot <- function(
  expr, spots, gene.name = NULL,
  binary = FALSE,
  log.expr = TRUE
)
{
  if (log.expr & !binary)
  {
    expr <- log(expr + 1)
  }

  ifelse(
    binary,
    binary.st.plot(expr, spots, gene.name = gene.name),
    abs.st.plot(expr, spots, gene.name = gene.name)
  )
}

## /////////////////////////////////////////////////////////////////////////////
## Internal functions ----
## /////////////////////////////////////////////////////////////////////////////

##' Spatial Transcriptomics Plot
##'
##' Plot a gene in a spatial transcriptomics (ST) dataset.
##'
##' @param count An \eqn{n}-by-\eqn{p} numeric matrix \eqn{Y} that denotes the
##'   gene expression count table. Each entry is the read count for
##'   gene \eqn{j} collected at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param gene.name An optional character string that specifies the gene
##'   to plot. The default is `NULL` to plot the gene with the largest
##'   count value.
##'
##' @return A `ggplot` object.
##'
##' @noRd
##' @keywords plot
##'
##' @importFrom ggplot2 ggplot aes geom_point labs scale_color_distiller
##'   coord_fixed theme_light theme element_blank element_text
##'
abs.st.plot <- function(expr, spots, gene.name = NULL)
{
  append.name   <- !is.null(gene.name)

  x <- spots[, 1]
  y <- spots[, 2]
  g <- expr

  plot(
    p <- ggplot(mapping = aes(x = x, y = y, color = g)) +
      geom_point(shape = 19, size = 6.5) +
      labs(
        title    = paste0(
          "Expression Levels",
          ifelse(
            append.name,
            paste0(" for '", gene.name, "'"),
            ""
          )
        ),
        subtitle = paste0(
          "Each point represents an entry in the ",
          "geospatial profile of the ST dataset."
        )
      ) +
      scale_color_distiller(palette = "RdYlBu") +
      coord_fixed() +
      theme_light() +
      theme(
        aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title   = element_text(face = 2)
      )
  )

  return(p)
}

##' Plot Dichotomised Expression Levels
##'
##' Plot the dichotomised expression levels for a given gene.
##'
##' @param expr A numeric vector \eqn{p} of length \eqn{n} that denotes the
##'   dichotomised gene expression levels. Each entry is one if the gene is
##'   highly expressed at spot \eqn{i} and zero otherwise.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param gene.name An optional character string that specifies the gene name.
##'   The default is `NULL` which appends nothing to the plot title.
##'
##' @return A `ggplot` object.
##'
##' @noRd
##' @keywords plot
##'
##' @importFrom ggplot2 ggplot aes geom_point labs scale_color_manual
##'   coord_fixed theme_light theme element_blank element_text
##'
binary.st.plot <- function(expr, spots, gene.name = NULL)
{
  append.name   <- !is.null(gene.name)
  invalid.input <- max(expr) > 1

  if (invalid.input)
  {
    stop("value passed to 'expr' should be the dichotomised expression levels")
  }

  x <- spots[, 1]
  y <- spots[, 2]
  g <- as.factor(expr)

  plot(
    p <- ggplot(mapping = aes(x = x, y = y, color = g)) +
      geom_point(shape = 19, size = 6.5) +
      labs(
        title    = paste0(
          "Dichotomized Expression Levels",
          ifelse(
            append.name,
            paste0(" for '", gene.name, "'"),
            ""
          )
        ),
        subtitle = paste0(
          "Each point represents an entry in the ",
          "geospatial profile of the ST dataset."
        )
      ) +
      scale_color_manual(
        values = c("0" = "#4575B4", "1" = "#D73027"),
        labels = c("Low-Expression Group", "High-Expression Group")
      ) +
      coord_fixed() +
      theme_light() +
      theme(
        aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title   = element_text(face = 2)
      )
  )

  return(p)
}

##' Spatial Transcriptomics Plot
##'
##' Plot a gene in a spatial transcriptomics (ST) dataset.
##'
##' @param count An \eqn{n}-by-\eqn{p} numeric matrix \eqn{Y} that denotes the
##'   gene expression count table. Each entry is the read count for
##'   gene \eqn{j} collected at spot \eqn{i}.
##' @param spots An \eqn{n}-by-\eqn{2} numeric matrix \eqn{T} to represent the
##'   geospatial profile, where each row indicates the spot location in
##'   the grid.
##' @param gene.name An optional character string that specifies the gene
##'   to plot. The default is `NULL` to plot the gene with the largest
##'   count value.
##'
##' @return A `ggplot` object.
##'
##' @noRd
##' @keywords plot
##'
##' @importFrom ggplot2 ggplot aes geom_point labs scale_color_distiller
##'   coord_fixed theme_light theme element_blank element_text
##'
st.plot_ <- function(count, spots, gene.name = NULL)
{
  genes <- colnames(count)
  total <- colSums(count)

  default.gene  <- is.null(gene.name)
  missing.names <- is.null(genes)

  if (missing.names)
  {
    stop("provide column names for the gene expression count table")
  }

  if (default.gene)
  {
    gene.name <- genes[which.max(total)]
  }
  else
  {
    missing.gene <- !(gene.name %in% genes)

    if (missing.gene)
    {
      stop("value passed to 'gene.name' is not in the expression count table")
    }
  }

  x <- spots[, 1]
  y <- spots[, 2]
  g <- count[, gene.name]

  plot(
    p <- ggplot(mapping = aes(x = x, y = y, color = g)) +
      geom_point(shape = 19, size = 6.5) +
      labs(
        title    = paste0("Expression Levels for '", gene.name, "'"),
        subtitle = paste0(
          "Each point represents an entry in the ",
          "geospatial profile of the ST dataset."
        )
      ) +
      scale_color_distiller(palette = "RdYlBu") +
      coord_fixed() +
      theme_light() +
      theme(
        aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title   = element_text(face = 2)
      )
  )

  return(p)
}
