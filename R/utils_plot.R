#' Scatterplot with color 
#' 
#' @param coords Coordinate table 
#' @param y Vector for colors 
#' @param label (optional) name for y vector.\cr
#' Default: NULL
#' @param color.low,color.high Low and high color for correlation values.\cr
#' Default: muted("blue") and muted("red"), respectively
#' 
#' @returns A ggplot object 
#' 
#' @export 
dimplot.generic <- function(coords, y, 
                            label = NULL, 
                            color.low = muted("blue"),
                            color.high = muted("red")) {
    plt_df <- cbind(coords, y) |> 
        as.data.frame() |> 
        # shuffle WITHOUT replacement
        dplyr::slice_sample(prop = 1, replace = FALSE) 
    plt <- ggplot2::ggplot(ggplot2::aes_string(colnames(coords)[1], 
                                               colnames(coords)[2], 
                                               color = 'y')) + 
        ggplot2::geom_point(size = 0.5) + 
        ggplot2::scale_color_gradient2(low = color.low, 
                                       high = color.high) + 
        ggplot2::theme_classic(base_size = 16)    
    if (!is.null(label)) 
        plt <- plt + ggplot2::labs(color = label)
    return(plt)
}

#' Scatter plot of correlation values based on FDR threshold 
#' 
#' @param coords A data.frame of coordinates (x,y). 
#' @param res A result output from [association] function
#' @param fdr_thresh A numeric FDR threshold (\alpha), values larger than this
#' are colored gray.
#' 
#' @returns A ggplot object 
#' 
#' @export 
dimplot.ncorr <- function(coords, res, 
                          fdr_thresh = 0.05) {
    .thresh <- subset(res$fdrs, fdr < fdr_thresh)$threshold[1]    
    plt_df <- cbind(coords, correlation = res$ncorrs) |> 
        as.data.frame() |> 
        dplyr::mutate(passed = abs(correlation) > .thresh) |> 
        dplyr::slice_sample(prop = 1, replace = FALSE)
    
    plt <- ggplot2::ggplot(plt_df, ggplot2::aes_string(colnames(coords)[1], 
                                                       colnames(coords)[2])) + 
        ggplot2::geom_point(data = subset(plt_df, !passed), 
                            size = .5, color = 'grey', alpha = .05) + 
        ggplot2::geom_point(data = subset(plt_df, passed), 
                            ggplot2::aes(color = correlation),
                            size = .5,  alpha = .5) + 
        ggplot2::scale_color_gradient2(low = muted("blue"), high = muted("red")) + 
        ggplot2::theme_classic(base_size = 16)
    return(plt)
}

