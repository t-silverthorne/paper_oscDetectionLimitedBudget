#' Helper function returns smallest eigenvalue
#' @export
getMinEig <- function(Mat,param){
  return(eigen(Mat,only.values=T) %>% {.$values} %>%  min())
}