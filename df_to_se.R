# authenticate with GCP using the token
# from a previous bigrquery::bq_auth()
do_gcs_auth <- function(){
  scope <-"https://www.googleapis.com/auth/cloud-platform"
  token <- gargle::token_fetch(scopes = scope)
  googleCloudStorageR::gcs_auth(token = token)  
}

# =========================================================
# taken from https://github.com/frederikziebell/frezielib

# is every observation containing only 
# the information from col1 and col2
# unique in the data.frame df?
check_unique <- function(df, col1, col2){
  !anyDuplicated(paste0(df[[col1]],df[[col2]]))
}

# is there a 1:n mapping between the values
# in col1 and the ones in col2
check_1_to_n <- function(df, col1, col2){
  df_unique <- dplyr::distinct(df[,c(col1, col2)])
  !anyDuplicated(df_unique[[col1]])
}

#' Convert data.frame to SummarizedExperiment
#' @param df The data.frame to convert
#' @param observation_id Name of the column that uniquely identifies observations (e.g. a sample ID)
#' @param feature_id Name of the column that uniquely identifies features (e.g. a gene ID)
#' @param value Character vector of one or more columns that contain measured values. Each element of
#' \code{value} will corresponds to a separate assay
#' @param observation_anno Character vector of columns that annotate observations (and will be put to the colData).
#' @param feature_anno Character vector of columns that annotate feature (and will be put to the rowData).
#' @param detect_colData_cols Whether to auto detect columns that can be added to the colData.
#' @param detect_rowData_cols Whether to auto detect columns that can be added to the rowData.
#' @param ambiguous Where to put ambiguous columns that can be put to both the colData and the rowData (requires
#' that \code{detect_colData_cols} and \code{detect_rowData_cols} are \code{TRUE}).
#' @param message_constant Whether to write a message listing which columns are constant.
df_to_se <- function(
    df,
    observation_id,
    feature_id,
    value,
    observation_anno = NULL,
    feature_anno = NULL,
    detect_colData_cols = TRUE,
    detect_rowData_cols = TRUE,
    ambiguous = c("both","none","rowData","colData"),
    message_constant = FALSE
) {
  
  ambiguous <- match.arg(ambiguous)
  
  # consistency checks
  if(!is.data.frame(df)){
    stop("df is not a data.frame.")
  }
  if(!observation_id %in% colnames(df)){
    stop("Column ", observation_col, " is not a column in the data frame.")
  }
  if(!feature_id %in% colnames(df)){
    stop("Column ", feature_id, " is not a column in the data frame.")
  }
  if(!all(value %in% colnames(df))){
    stop("One of the columns ", paste0(value, collapse=", "), " is not a column in the data frame.")
  }
  if(!is.null(observation_anno)){
    if(!all(observation_anno %in% colnames(df))){
      stop("One of the columns ", paste0(observation_anno, collapse=", "), " is not a column in the data frame.")
    }
  }
  if(!is.null(feature_anno)){
    if(!all(feature_anno %in% colnames(df))){
      stop("One of the columns ", paste0(feature_anno, collapse=", "), " is not a column in the data frame.")
    }
  }
  if(anyNA(df[[observation_id]])){
    stop("Observation ID column ",observation_id, " contains NAs.")
  }
  if(anyNA(df[[feature_id]])){
    stop("Feature ID column ",feature_id, " contains NAs.")
  }
  
  # check if there is exactly one row per (observation_id, feature_id) combination
  if(check_unique(df, feature_id, observation_id) == FALSE){
    stop(
      "The combination of ", feature_id," and ", observation_id, " is not unique.",
      "\nThis means there are at least two rows with the same entries in ",feature_id," and ",observation_id,"."
    )
  }
  
  # check if requested observation annotation has a 1:n mapping with observation_id
  if(!is.null(observation_anno)){
    for(col in observation_anno){
      if(check_1_to_n(df, observation_id, col)==FALSE){
        stop("Column ",col, " cannot be added to colData because there is no 1:n mapping between the values in column ",observation_id, " and the ones in column ", col,".")
      }
      
    }
    col_data_cols <- c(observation_id, observation_anno)
  } else {
    col_data_cols <- observation_id
  }
  
  # check if requested feature annotation has a 1:n mapping with feature_id
  if(!is.null(feature_anno)){
    for(col in feature_anno){
      if(check_1_to_n(df, feature_id, col)==FALSE){
        stop("Column ",col, " cannot be added to rowData because there is no 1:n mapping between the values in column ",observation_id, " and the ones in column ", col,".")
      }
    }
    row_data_cols <- c(feature_id, feature_anno)
  } else {
    row_data_cols <- feature_id
  }
  
  # remaining columns in df that are not observations, features, values
  # or already to be added to colData or rowData
  rest_cols <- setdiff(
    colnames(df), 
    c(feature_id, observation_id, value, observation_anno, feature_anno)
  )
  
  # notify about constant columns
  if(message_constant == TRUE) {
    constant_cols <- unlist(
      lapply(rest_cols, function(col){
        if(length(unique(df[[col]]))==1){
          col
        }
      })
    )
    if(length(constant_cols) > 0){
      message("Columns ", paste0(constant_cols, collapse=", "), " are constant.")
    }
  }
  
  # find out which rest_cols can be put to the colData
  rest_cols_col_data <- c()
  if(detect_colData_cols == TRUE){
    for(col in rest_cols){
      if(check_1_to_n(df, observation_id, col)){
        rest_cols_col_data <- c(rest_cols_col_data, col)
      }
    }
  }
  
  # find out which rest_cols can be put to the colData
  rest_cols_row_data <- c()
  if(detect_rowData_cols==TRUE){
    for(col in rest_cols){
      if(check_1_to_n(df, feature_id, col)){
        rest_cols_row_data <- c(rest_cols_row_data, col)
      }
    }
  }
  
  # put unambiguous rest_col to the respective rowData and colData
  col_data_cols <- c(col_data_cols, setdiff(rest_cols_col_data, rest_cols_row_data))
  row_data_cols <- c(row_data_cols, setdiff(rest_cols_row_data, rest_cols_col_data))
  # identify not distributed and ambiguous rest cols
  cols_ambiguous <- intersect(rest_cols_col_data, rest_cols_row_data)
  cols_not_distributed <- setdiff(rest_cols, c(rest_cols_col_data, rest_cols_row_data))
  
  # inform about columns that go neither to the colData
  # nor to the rowData
  if(length(cols_not_distributed) > 0 & detect_colData_cols == TRUE & detect_rowData_cols == TRUE){
    message("Columns ", paste0(cols_not_distributed, collapse = ", "), " could not be distributed to the colData or rowData. Are these additional assays?")
  }
  
  # deal with ambiguous columns
  if(length(cols_ambiguous) > 0) {
    if(ambiguous == "both"){
      message("Columns ", paste0(cols_ambiguous, collapse = ", "), " are ambiguous and were distributed to the colData and the rowData.")
      col_data_cols <- c(col_data_cols, cols_ambiguous)
      row_data_cols <- c(row_data_cols, cols_ambiguous)
    }
    if(ambiguous == "rowData"){
      row_data_cols <- c(row_data_cols, cols_ambiguous)
    }
    if(ambiguous == "colData"){
      col_data_cols <- c(col_data_cols, cols_ambiguous)
    }
  }
  
  # build data.frame with column info
  col_data <- dplyr::distinct(df[,col_data_cols,drop=F])
  # turn observation_id into rownames
  rownames(col_data) <- col_data[[observation_id]]
  col_data <- col_data[,!colnames(col_data) %in% observation_id, drop = F]
  
  # build data.frame with row info
  row_data <- dplyr::distinct(df[,row_data_cols,drop=F])
  # turn feature_id into rownames
  rownames(row_data) <- row_data[[feature_id]]
  row_data <- row_data[,!colnames(row_data) %in% feature_id, drop = F]
  
  # build assay matrices
  mat_list <- lapply(value, function(value_col){
    mat <- reshape2::acast(
      data = df, 
      formula = as.formula(paste0(feature_id, "~", observation_id)),
      value.var = value_col
    )
    # arrange assay matrix to the order 
    # given in rowData and colData
    mat[rownames(row_data),rownames(col_data)]
  })
  names(mat_list) <- value
  
  # create SE object
  SummarizedExperiment(
    assays = mat_list, 
    colData = col_data, 
    rowData = row_data
  )
}
# =========================================================
