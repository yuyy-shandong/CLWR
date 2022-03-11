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

clwr_coff <- function(data = data, y1 = "y1",
    y2 = "y2", y3 = "y3", x1 = "x1",
    x2 = "x2", x3 = "x3", centre = "centre") {

    data <- data[, c(y1, y2, y3, x1,
        x2, x3, centre)]
    colnames(data) <- c("y1", "y2", "y3",
        "x1", "x2", "x3", "centre")

    data <- na.omit(data)
    cy3x1 <- NULL
    cx2x1 <- NULL
    cy2x1 <- NULL
    cx3y1 <- NULL
    cy2y1 <- NULL
    cx2y1 <- NULL
    big_se <- NULL
    big_coff <- NULL

    centr_num <- unique(data$centre)

    len <- length(centr_num)
    for (cen in 1:len) {
        data_s <- data[data$centre ==
            centr_num[cen], ]
        ### 1
        model1 <- lm(y3 ~ x1, data = data_s)
        cy3x1[cen] <- summary(model1)$coefficients[2,
            1]
        ### 2
        model2 <- lm(x2 ~ x1, data = data_s)
        cx2x1[cen] <- summary(model2)$coefficients[2,
            1]
        ### 3
        model3 <- lm(y2 ~ x1, data = data_s)
        cy2x1[cen] <- summary(model3)$coefficients[2,
            1]
        # 2-1
        model11 <- lm(x3 ~ y1, data = data_s)
        cx3y1_y <- summary(model11)$coefficients[2,
            1]
        cx3y1[cen] <- cx3y1_y * var(data_s$y1) / var(data_s$x1)

        ### 2-2
        model21 <- lm(y2 ~ y1, data = data_s)
        cy2y1_y <- summary(model21)$coefficients[2,
            1]
        cy2y1[cen] <- cy2y1_y * var(data_s$y1) / var(data_s$x1)

        ### 2-3
        model31 <- lm(x2 ~ y1, data = data_s)
        cx2y1_y <- summary(model31)$coefficients[2,
            1]
        cx2y1[cen] <- cx2y1_y * var(data_s$y1) / var(data_s$x1)

        big_coff[cen] <- cy3x1[cen] -
            cx3y1[cen]

        cy3x1_boot <- NULL
        cx3y1_boot <- NULL
        big_coff_boot <- NULL

        for (boot_i in 1:1000) {
            sample_new <- sample(1:dim(data_s)[1],
                dim(data_s)[1], replace = T)
            data_boot <- data_s[sample_new,
                ]

            ### 1
            model1_boot <- lm(y3 ~ x1,
                data = data_boot)
            cy3x1_boot[boot_i] <- summary(model1_boot)$coefficients[2,
                1]
            model11_boot <- lm(x3 ~ y1,
                data = data_boot)
            cx3y1_y_boot <- summary(model11_boot)$coefficients[2,
                1]
            cx3y1_boot[boot_i] <- cx3y1_y_boot *
                var(data_boot$y1) / var(data_boot$x1)
            big_coff_boot[boot_i] <- cy3x1_boot[boot_i] -
                cx3y1_boot[boot_i]
        }

        big_se[cen] <- sd(big_coff_boot,
            na.rm = T)
    }

    result_data <- data.frame(cy3x1,
        cx2x1, cy2x1, cx3y1, cy2y1, cx2y1,
        big_se, big_coff)

    return(result_data)

}
