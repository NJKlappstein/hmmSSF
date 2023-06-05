##' Get step lengths and turning angles in a case/control dataframe
##'
##' @param data Data frame with columns 'obs' (1: observation, 0: control),
##' 'x', and 'y' (and possibly others)
##'
##' @return Data frame passed as input with two extra columns: step and angle
##'
##' @export

get_step_and_angle <- function(data) {
    ## Observed locations
    obs <- subset(data, obs == 1)
    n_obs <- nrow(obs)

    ## Number of control steps per observation (must always be the same)
    n_zeros <- nrow(data)/n_obs - 1

    ## Start points of steps (include NAs at the top because no
    ## step length for first stratum)
    start_step <- matrix(
        rep(as.matrix(rbind(NA, obs[-n_obs, c("x", "y")])),
            each = n_zeros + 1),
        ncol = 2)

    ## Start points of angles (include NAs for the top two rows
    ## because no turning angle for first two strata)
    start_angle <- matrix(
        rep(as.matrix(rbind(NA, NA, obs[-c(n_obs - 1, n_obs), c("x", "y")])),
            each = n_zeros + 1),
        ncol = 2)

    ## Mid points of angles
    mid_angle <- matrix(
        rep(as.matrix(rbind(NA, NA, obs[-c(1, n_obs), c("x", "y")])),
            each = n_zeros + 1),
        ncol = 2)

    ## End points of steps and angles
    end <- as.matrix(data[,c("x", "y")])

    ## Step lengths
    steps <- sqrt(rowSums((end - start_step)^2))

    ## Vector for previous step
    previous_vector <- mid_angle - start_angle
    previous_bearing <- atan2(previous_vector[,2], previous_vector[,1])
    ## Vector for this step
    this_vector <- end - mid_angle
    this_bearing <- atan2(this_vector[,2], this_vector[,1])
    ## Turning angle between previous step and this step
    angles <- this_bearing - previous_bearing

    # Keep angles between -pi and pi
    angles[which(angles <= (-pi))] <- angles[which(angles <= (-pi))] + 2*pi
    angles[which(angles > pi)] <- angles[which(angles > pi)] - 2*pi

    data <- cbind(data, step = steps, angle = angles)

    first_locs <- which(data$ID[-1] != data$ID[-nrow(data)])+1
    data$step[first_locs] <- NA
    data$angle[first_locs] <- NA
    second_strat <- data$stratum[first_locs] + 1
    data$angle[which(data$stratum %in% second_strat)] <- NA

    return(data)
}
