test_Rwanderlust <- function() {
    set.seed(15)
    shuffled_iris <- iris[sample(150, 150, replace = FALSE), ]
    data <- shuffled_iris[ ,1:4]
    wishbone <- Rwanderlust(data = data, num_waypoints = 100, waypoints_seed = 2)
    checkEquals(length(wishbone), 3)
}
