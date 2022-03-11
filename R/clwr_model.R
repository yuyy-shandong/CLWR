#' Eliminating the Unmeasured Confounders and Estimating
#' bidirectional causal effect and lagged causal effect
#'
#' @param data an optional data frame containing the variables in the model.
#' @param x1 the variable name of exposure at time t1
#' @param x2 the variable name of exposure at time t2
#' @param x3 the variable name of exposure at time t3
#' @param y1 the variable name of outcome at time t1
#' @param y2 the variable name of outcome at time t2
#' @param y3 the variable name of outcome at time t3
#' @param centre the variable name of centre
#'
#'
#' @return  \code{est_xy} the casual effect estimation
#' and it's 95%CI of X at time t2 on Y at time t3
#' @return  \code{est_xx} the casual effect estimation
#' and it's 95%CI of X at time t2 on X at time t3
#' @return  \code{est_yx} the casual effect estimation
#' and it's 95%CI of Y at time t2 on X at time t3
#' @return  \code{est_yy} the casual effect estimation
#' and it's 95%CI of Y at time t2 on Y at time t3

clwr_model <- function(data = data, y1 = "y1",
    y2 = "y2", y3 = "y3", x1 = "x1",
    x2 = "x2", x3 = "x3", centre = "centre") {
    coff_data <- clwr_coff(data = data,
        y1 = "y1", y2 = "y2", y3 = "y3",
        x1 = "x1", x2 = "x2", x3 = "x3",
        centre = "centre")

    colnames(coff_data) <- c("cy3x1",
        "cx2x1", "cy2x1", "cx3y1", "cy2y1",
        "cx2y1", "big_se", "big_coff")

    data <- na.omit(data)
    centr_num <- unique(data$centre)

    nn_centre <- length(centr_num)
    we_matrix <- matrix(rep(0, nn_centre *
        nn_centre), nrow = nn_centre,
        ncol = nn_centre)
    we_matrix2 <- matrix(rep(0, nn_centre *
        nn_centre), nrow = nn_centre,
        ncol = nn_centre)
    cen_data <- matrix(rep(1, nn_centre *
        nn_centre), nrow = nn_centre,
        ncol = nn_centre)
    cen_data[upper.tri(cen_data)] <- 0
    cen_data[lower.tri(cen_data)] <- 0
    for (ppp in 1:nn_centre) {
        for (qqq in 1:nn_centre) {
            we_matrix[ppp, qqq] <- (1 / coff_data$big_se[ppp]) *
                (1 / coff_data$big_se[qqq]) *
                cen_data[ppp, qqq]
            we_matrix2[ppp, qqq] <- (1 / coff_data$big_se[ppp]) *
                (1 / coff_data$big_se[qqq]) *
                cen_data[ppp, qqq]

        }
    }
    xo_xo <- data.frame(1, coff_data$cx2x1,
        coff_data$cy2x1, -coff_data$cy2y1,
        -coff_data$cx2y1)
    xo_xo <- as.matrix(xo_xo)
    omega <- we_matrix
    beta1 <- solve(t(xo_xo) %*% omega %*%
        xo_xo) %*% t(xo_xo) %*% omega %*%
        coff_data$big_coff

    weigt_est_xy <- beta1[2]
    weigt_est_yy <- beta1[3]
    weigt_est_yx <- beta1[4]
    weigt_est_xx <- beta1[5]
    all_coef <- data.frame(coff_data$big_coff,
        coff_data$cx2x1, coff_data$cy2x1,
        -coff_data$cy2y1, -coff_data$cx2y1)
    weigt_est_xy_boot <- NULL
    weigt_est_yy_boot <- NULL
    weigt_est_yx_boot <- NULL
    weigt_est_xx_boot <- NULL

    for (boo in 1:1000) {
        hanghao <- sample(1:nn_centre,
            nn_centre, replace = T)
        we_matrix_boot <- matrix(rep(0,
            nn_centre * nn_centre), nrow = nn_centre,
            ncol = nn_centre)
        for (ppp in 1:nn_centre) {
            we_matrix_boot[ppp, ppp] <- we_matrix[hanghao[ppp],
                hanghao[ppp]]
        }
        data_boot <- all_coef[hanghao,
            ]

        xo_xo_boot <- data.frame(1, data_boot[,
            2], data_boot[, 3], data_boot[,
            4], data_boot[, 5])
        xo_xo_boot <- as.matrix(xo_xo_boot)
        omega_boot <- we_matrix_boot
        beta1_boot <- try(solve(t(xo_xo_boot) %*%
            omega_boot %*% xo_xo_boot) %*%
            t(xo_xo_boot) %*% omega_boot %*%
            data_boot[, 1], silent = T)
        if ("try-error" %in% class(beta1_boot)) {
            (next)()
        }
        weigt_est_xy_boot[boo] <- beta1_boot[2]
        weigt_est_yy_boot[boo] <- beta1_boot[3]
        weigt_est_yx_boot[boo] <- beta1_boot[4]
        weigt_est_xx_boot[boo] <- beta1_boot[5]
    }


    weigt_est_xy_ci <- quantile(weigt_est_xy_boot,
        c(0.025, 0.975), na.rm = T)
    weigt_est_yy_ci <- quantile(weigt_est_yy_boot,
        c(0.025, 0.975), na.rm = T)
    weigt_est_yx_ci <- quantile(weigt_est_yx_boot,
        c(0.025, 0.975), na.rm = T)
    weigt_est_xx_ci <- quantile(weigt_est_xx_boot,
        c(0.025, 0.975), na.rm = T)
    result <- list()
    est_xy <- c(weigt_est_xy, weigt_est_xy_ci)
    names(est_xy) <- c(paste0("The causal effect of ",
        x2, " on ", y3), "est_2.5%",
        "est_97.5%")
    result$est_xy <- est_xy

    est_xx <- c(weigt_est_xx, weigt_est_xx_ci)
    names(est_xx) <- c(paste0("The causal effect of ",
        x2, " on ", x3), "est_xx_2.5%",
        "est_xx_97.5%")
    result$est_xx <- est_xx

    est_yx <- c(weigt_est_yx, weigt_est_yx_ci)
    names(est_yx) <- c(paste0("The causal effect of ",
        y2, " on ", x3), "est_yx_2.5%",
        "est_yx_97.5%")
    result$est_yx <- est_yx

    est_yy <- c(weigt_est_yy, weigt_est_yy_ci)
    names(est_yy) <- c(paste0("The causal effect of ",
        y2, " on ", y3), "est_yy_2.5%",
        "est_yy_97.5%")
    result$est_yy <- est_yy

    return(result)

}
