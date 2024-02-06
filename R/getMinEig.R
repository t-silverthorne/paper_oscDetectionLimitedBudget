#' Helper function returns smallest eigenvalue
#' @export
getMinEig <- function(Mat,is_symmetric=F){
  return(eigen(Mat,only.values=T,symmetric = is_symmetric) %>% {.$values} %>%  min())
}