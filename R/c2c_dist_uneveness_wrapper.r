#' @author MaiChan Lau
#' @export
c2c_dist_uneveness_wrapper <- function(input_expr = NULL,
                                       alpha = NULL, data_type = 'linear'){
    
    d <-distance.function(input_expr)
    if(data_type == 'linear') {
        d_var <-c2c_dist_uneveness_linear(d, alpha=alpha)
    }else if(data_type == 'cyclical'){
        d_var <- c2c_dist_uneveness_cyclical(d, alpha=alpha)
    }
    return (d_var)
    
}