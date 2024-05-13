#' Accurate Posterior Predictive Plots for Ordbetareg Models
#'
#' The standard [brms::pp_check] plot available via `brms`
#' is not accurate for ordbetareg models because an
#' ordered beta regression has both continuous and discrete
#' components. This function implements a bar plot and a density
#' plot for the continuous and discrete elements separately,
#' and will return accurate posterior predictive plots
#' relative to the data.
#'
#' @param model A fitted [ordbetareg] model.
#' @param type Default is "both" for creating both a
#' discrete (bar) and continuous (density) plot. Can also be
#' "discrete" for only the bar plot for discrete values (0/1) or
#' "continuous" for continuous values (density plot).
#' @param ndraws Number of posterior draws to use to calculate estimates and show in plot.
#' Defaults to 10.
#' @param cores Number of cores to use to produce posterior predictive distribution. Defaults to NULL or 1 core.
#' @param group A factor variable of the same number of
#' rows as the data that is used to broduce grouped
#' (faceted) plots of the posterior distribution.
#' @param new_theme Any additional themes to be added to ggplot2 (default is NULL).
#' @param outcome_label A character value that will replace the name of the outcome in the plot
#' (default is the name of the response variable in the data frame).
#' @param animate Whether to animate each posterior draw for continuous
#' distributions (defaults to FALSE). Requires installation of the
#' `gganimate` and `transformr` R packages.
#' @param reverse_bounds Whether to plot data using the original bounds in the data
#' (i.e. not 0 and 1).
#' @param facet_scales The option passed on to the `facet_wrap` function in
#'  `ggplot2` for the type of scale for facetting if passing a variable for
#'  `group`. Defaults to `"fixed"` scales but can be set to `"free_y"` to allow
#'  probability density/bar count scales to vary or `"free"` to allow both x
#'  and y axes to vary (i.e., also outcome axis ticks).
#' @return If "both", prints both plots and returns a list of both plots as
#' `ggplot2` objects. Otherwise, prints and returnst the specific plot
#' as a `ggplot2` object.
#' @examples
#'
#' # need a fitted ordbetareg model
#'
#' data("ord_fit_mean")
#'
#' out_plots <- pp_check_ordbeta(ord_fit_mean)
#'
#' # view discrete bar plot
#'
#' out_plots$discrete
#'
#' # view continuous density plot
#'
#' out_plots$continuous
#'
#' # change title using ggplot2 ggtitle function
#'
#' out_plots$discrete + ggplot2::ggtitle("New title")
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr select
#' @importFrom utils packageVersion
pp_check_ordbeta <- function(model=NULL,
                             type="both",
                             ndraws=10,
                             cores=NULL,
                             group=NULL,
                             new_theme=NULL,
                             outcome_label=NULL,
                             animate=FALSE,
                             reverse_bounds=TRUE,
                             facet_scales="fixed") {

    #define globals

    pred <- draw <- true <- NULL

    # need to get posterior predictive distribution

    full_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL)

    outcome <- model$data[[1]]

    if(!is.null(outcome_label)) {

      outcome_label <- outcome_label

    } else {

      outcome_label <- names(model$data)[1]

    }

    this_star <- transformr::poly_star()

    if(reverse_bounds) {

        # revert to original scale

        l_bound <- model$lower_bound
        up_bound <- model$upper_bound


        full_dist <- apply(full_dist, 1, function(c) {

          c <- (c *(up_bound - l_bound)) + l_bound

          c

        }) %>% t

        outcome <- (outcome *(up_bound - l_bound)) + l_bound
        # deal with floating point errors

        tol <- 1e-5

        outcome[abs(outcome - up_bound) < tol] <- up_bound
        outcome[abs(outcome - l_bound) < tol] <- l_bound

    } else {

      l_bound <- 0
      up_bound <- 1

    }

    output <- list()

    if(is.null(group)) {

      group <- rep(1, length(outcome))

    }

    if(type %in% c("both","discrete")) {


          plot_data_bar <- lapply(unique(group), function(g) {

                apply(full_dist[,group==g], 1, function(c) {

                  tibble(num_ones=sum(c==up_bound),
                         num_zeroes=sum(c==l_bound),
                         num_cont=sum(c<up_bound & c>l_bound),
                         group=g)

                }) %>% bind_rows

            }) %>% bind_rows

          plot_data_bar$iter <- 1:nrow(plot_data_bar)

          plot_data_bar <- bind_rows(select(plot_data_bar, var="num_ones", group) %>% mutate(true=sum(outcome==up_bound,na.rm=T),
                                                                                             type=as.character(round(up_bound,3))),
                                     select(plot_data_bar, var="num_zeroes", group) %>% mutate(true=sum(outcome==l_bound,na.rm=T),
                                                                                               type=as.character(round(l_bound,3))),
                                     select(plot_data_bar, var="num_cont", group) %>% mutate(true=sum(outcome>l_bound & outcome<up_bound,na.rm=T),
                                                                                                   type=paste0('(',round(l_bound,3),
                                                                                                               ',',round(up_bound,3),
                                                                                                               ')'))) %>%
                          mutate(type=factor(type,
                                             levels=c(round(l_bound,3),
                                                      paste0('(',round(l_bound,3),
                                                             ',',round(up_bound,3),
                                                             ')'),
                                 round(up_bound,3))))

          true_data_bar <- lapply(unique(group), function(g) {

                        tibble(type=c(round(l_bound,3),
                                         round(up_bound,3),
                                         paste0('(',round(l_bound,3),
                                                                 ',',round(up_bound,3),
                                                                 ')')),
                                  true=c(sum(outcome[group==g]==l_bound,na.rm=T),
                                         sum(outcome[group==g]==up_bound,na.rm=T),
                                         sum(outcome[group==g]>l_bound & outcome[group==g]<up_bound,na.rm=T))) %>%
                                  mutate(type=factor(type,
                                                     levels=c(round(l_bound,3),
                                                              paste0('(',round(l_bound,3),
                                                                     ',',round(up_bound,3),
                                                                     ')'),
                                                              round(up_bound,3))),
                                         group=g)

                            }) %>% bind_rows



          bar_plot <- ggplot(plot_data_bar,
                             aes(x=type,y=var)) +
            geom_col(data=true_data_bar, aes(y=true),width=.5,fill="gray") +
            stat_summary(fun.data="median_hilow",geom = "pointrange") +
            labs(y="Observed and Estimated Counts",
                 x=paste0("No. of Discrete and Continuous Observations in ", outcome_label),
                 caption="Estimates are posterior medians with 5% and 95% high/low quantiles.\nCount of discrete values in data are represented by gray bars.") +
            ggtitle("Posterior Predictions for Types of Responses in Data",
                    subtitle=paste0("Outcome: ",outcome_label)) +
            theme_minimal()

          # add in facets if groups>1

          if(length(unique(plot_data_bar$group))>1) {

            bar_plot <- bar_plot + facet_wrap(~group,
                                              scales=facet_scales)
          }


          output$discrete <- bar_plot + new_theme

    }

    if(type %in% c("both","continuous")) {

      # need to change the posterior predict function to get the argument from ...

      cont_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL,ntrys=100)

      if(reverse_bounds) {

        # revert to original scale

        l_bound <- model$lower_bound
        up_bound <- model$upper_bound


        cont_dist <- apply(cont_dist, 1, function(c) {

          c <- (c *(up_bound - l_bound)) + l_bound

          c

        }) %>% t

      } else {

        l_bound <- 0
        up_bound <- 1

      }

      # if prediction was continuous, get continuous value from the outcome

      plot_data_dens <- lapply(unique(group), function(g) {

        counter <- 0

        apply(cont_dist[,group==g], 1, function(col) {

          counter <<- counter + 1

          tibble(pred=col[!(col %in% c(l_bound,up_bound))],
                          draw=counter,
                 group=g)

        }) %>% bind_rows

      }) %>% bind_rows

      true_data <- lapply(unique(group), function(g) {

        tibble(true=outcome[!(outcome %in% c(l_bound,up_bound)) & group==g]) %>%
          mutate(group=g)

      }) %>% bind_rows

      if(animate) {

        if (!requireNamespace("gganimate", quietly = TRUE)) {
          stop("Please install the `gganimate` package to use the animate option in this plot.", call. = FALSE)
        }

        if (!requireNamespace("transformr", quietly = TRUE)) {
          stop("Please install the `transformr` package to use the animate option in this plot.", call. = FALSE)
        }


        cont_plot <- plot_data_dens %>%
          ggplot(aes(x=pred)) +
          geom_density(alpha=0.7,bounds=c(l_bound,up_bound)) +
          theme_minimal() +
          labs(y="Probability Density",x=paste0("Continuous values of ",outcome_label),
               caption="Each black line is a draw from the posterior distribution.\nData distribution represented by gray line.") +
          gganimate::transition_time(draw) +
          gganimate::ease_aes('elastic-in-out') +
          ggtitle(paste0("Posterior Predictions for Data Bounded from ",
                         round(l_bound,3), " to ", round(up_bound,3)),
                  subtitle=paste0("Outcome: ",outcome_label))

      } else {

        cont_plot <- plot_data_dens %>%
          ggplot(aes(x=pred,group=draw)) +
          geom_density(alpha=0.7,bounds=c(l_bound,up_bound)) +
          theme_minimal() +
          labs(y="Probability Density",x=paste0("Continuous values of ",outcome_label),
               caption="Each black line is a draw from the posterior distribution.\nData distribution represented by gray line.") +
          ggtitle(paste0("Posterior Predictions for Data Bounded from ",
                         round(l_bound,3), " to ", round(up_bound,3)),
                  subtitle=paste0("Outcome: ",outcome_label))


      }

      # add in facets if groups>1

      if(length(unique(plot_data_bar$group))>1) {

        cont_plot <- cont_plot + facet_wrap(~group,
                                            scales=facet_scales)
      }

      cont_plot <- cont_plot + geom_density(data=true_data,
                                            mapping=aes(x=true),
                                            linewidth=2,colour="gray",alpha=0.7,
                                            bounds=c(l_bound,up_bound)) +
        new_theme

      output$continuous <- cont_plot


    }

    return(output)

}

