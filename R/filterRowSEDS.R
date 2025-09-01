
#' @export
filterRowSEDS <- function(se, threshold = 1, nsample = 1, percent = NULL){

  if(!is.null(nsample) & !is.null(percent)){

    stop(sprintf("Choose one filtering criteria."))

  } else if(is.null(nsample) & is.null(percent)){

    stop(sprintf("Choose one filtering criteria."))

  } else if(!is.null(nsample)){  # filtering on the number of samples

    if(is.null(SummarizedExperiment::assay(se, "counts"))){
      stop(sprintf("No assay counts in the SE object."))
    }
    # check counts assay
    if (!("counts" %in% SummarizedExperiment::assayNames(se))) {
      stop("No assay 'counts' in the SummarizedExperiment object.")
    }
    # extracting counts assay
    mat <- SummarizedExperiment::assay(se, "counts")
    count_over_threshold <- rowSums(mat > as.numeric(threshold))
    keep <- count_over_threshold >= as.numeric(nsample)

    # filtering SummarizedExperiment
    se_filtered <- se[keep, ]

    return(se_filtered)

  } else {    # filtering on a percentage of samples

    # check counts assay
    if (!("counts" %in% SummarizedExperiment::assayNames(se))) {
      stop("No assay 'counts' in the SummarizedExperiment object.")
    }
    # estimate the number of samples
    nsample <- ((as.numeric(percent)/100)*ncol(se))

    mat <- SummarizedExperiment::assay(se, "counts")
    count_over_threshold <- rowSums(mat > as.numeric(threshold))
    keep <- count_over_threshold >= as.numeric(nsample)

    # filtering SummarizedExperiment
    se_filtered <- se[keep, ]

    return(se_filtered)

  }
}
