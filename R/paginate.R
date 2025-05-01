paginate <- function(x, n_per_page){
  if(length(x)>1){
    id <- seq_along(x)
  }else if(length(x) %in% 1){
    if(is.wholenumber(x)){
      id <- seq(from = 1, to = x)
    }else{
      id <- 1L
    }
  }else if(length(x) %in% 0){
    id <- integer(length = 0L)
  }

  if(length(id)>0){
  n_pages <- ceiling(max(id)/n_per_page)

  pp <- rep(1:n_pages,
            each = n_per_page)
  pp[id]
  }else{
    numeric(0)
  }

}
