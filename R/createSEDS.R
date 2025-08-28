#'@export
createSEDS <- function(rnaseq, metadata, coldata_column = NULL, rowData = NULL, rowdata_column = NULL){

  # column set to id rownames
  if(!is.null(coldata_column)){

    colData <- S4Vectors::DataFrame(metadata)
    rownames(colData) <- colData[,as.character(coldata_column)]

  } else{
    # matching id used are rownames
    colData <- S4Vectors::DataFrame(metadata)

  }

  #filtering id not matching and casting to matrix
  rnaseq <- as.matrix(rnaseq[,which(colnames(rnaseq) %in% rownames(colData))])

  # if existing rowData and matching id column
  if(!is.null(rowdata_column) & !is.null(rowData)){

    rowData <- S4Vectors::DataFrame(rowData)
    rownames(rowData) <- rowData[,as.character(rowdata_column)]

  }
  if(!is.null(rowData)){

    # creation SE object with rowdata
    rse <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=rnaseq),
                                                      rowData=rowData,
                                                      colData=colData)
  }else{

    # creation SE object without rowdata
    rse <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=rnaseq),
                                                      colData=colData)
  }

  return(rse)
}
