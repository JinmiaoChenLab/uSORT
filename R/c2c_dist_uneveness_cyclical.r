#' A c2c distance uneveness function
#' 
#' A c2c distance uneveness function allows you to compute the overall cell-to-cell (c2c)
#' distance variability or uneveness which is summed over small local windows surrounding 
#' individual cells.
#' 
#' @param d an expresssion matrix containing n-by-n pairwise distances computed for n cells
#' @param alpha window size which is given in the percentage of total no. of cells
#' @return a value of total uneveness 
#' @author MaiChan Lau
#' @export



c2c_dist_uneveness_cyclical<-function(d, alpha=0.3){
    
    ## Initialization
    num_cell<-nrow(d)
    total.var<-0
    window_size<-ceiling(num_cell*alpha)
    
    ## Compute uneveness for each cell
#     for(i in 1:num_cell){
#         window_around_cell<-d[-i,i,drop=F]
#         startID <- i - (as.integer(window_size/2))
#         endID <- i + (as.integer(window_size/2))
#         if(i<floor(window_size/2)+1) {
#             startID<-i
#             endID<-i+window_size
#         }
#         else{
#             startID<-floor(i-window_size/2)
#             endID<-ceiling(i+window_size/2)
#         }
#         
#         if(endID>(num_cell-1))endID<-(num_cell-1)
#         window_around_cell<-window_around_cell[startID:endID,,drop=F]
#         
#         ind.cell.var<-var(window_around_cell)
#         total.var<-total.var+ind.cell.var
#     }
    
    #for(i in (as.integer(window_size/2)):(num_cell-(as.integer(window_size/2)))){
    for(i in 1:num_cell){
        window_around_cell<-d[i,]
        window_around_cell_extended <- c(window_around_cell,window_around_cell,window_around_cell)
        
        startID <- i - (as.integer(window_size/2))
        endID <- i + (as.integer(window_size/2))

        if(startID<1){
            extended_startID <- num_cell-i
        }else{
            extended_startID <- num_cell+i}
        
        extended_endID <- num_cell + endID
        
        #if(startID<1)startID<-1
        #if(endID>num_cell) endID <- num_cell
        
        #cat(i,'\t',startID,'\t',endID,'\n')
        #cat(i,'\t',extended_startID,'\t',extended_endID,'\n')
        
        window_around_cell_extended<-window_around_cell_extended[extended_startID:extended_endID]
        window_around_cell_extended <- window_around_cell_extended[-i]
        ind.cell.var<-var(window_around_cell_extended)
        total.var<-total.var+ind.cell.var
        #cat(i, '\t', 'ind.cell.var = ', ind.cell.var, '\n')
    }
    
    return(total.var)
}
