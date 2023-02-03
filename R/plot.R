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
#' @param theme Any additional theme arguments to be added to ggplot2 (default is NULL).
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
                             theme=NULL) {

    #define globals

    pred <- draw <- true <- NULL

    # need to get posterior predictive distribution

    full_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL)

    outcome <- model$data[[1]]

    output <- list()

    if(is.null(group)) {

      group <- rep(1, length(outcome))

    }

    if(type %in% c("both","discrete")) {


          plot_data_bar <- lapply(unique(group), function(g) {

                apply(full_dist[,group==g], 1, function(c) {

                  tibble(num_ones=sum(c==1),
                         num_zeroes=sum(c==0),
                         group=g)

                }) %>% bind_rows

            }) %>% bind_rows

          plot_data_bar <- bind_rows(select(plot_data_bar, var="num_ones", group) %>% mutate(true=sum(outcome==1,na.rm=T),type=1),
                                     select(plot_data_bar, var="num_zeroes", group) %>% mutate(true=sum(outcome==0,na.rm=T),type=0))
          true_data_bar <- tibble(type=c(0,1),
                                  true=c(sum(outcome==0,na.rm=T),sum(outcome==1,na.rm=T)))
          bar_plot <- ggplot(plot_data_bar,
                             aes(x=type,y=var)) +
            geom_col(data=true_data_bar, aes(y=true),width=.5,fill="gray") +
            stat_summary(fun.data="median_hilow",geom = "pointrange") +
            labs(y="Observed and Estimated Counts",
                 x=paste0("Lower/Higher Discrete Bounds in ", names(model$data)[1]),
                 caption="Estimates are posterior medians with 5% and 95% high/low quantiles.\nCount of discrete values in data are represented by gray bars.") +
            ggtitle("Posterior Predictions for Discrete Values",
                    subtitle=paste0("Outcome: ",names(model$data)[1])) +
            theme_minimal()

          # add in facets if groups>1

          if(length(unique(plot_data_bar$group))>1) {

            bar_plot <- bar_plot + facet_wrap(~group)
          }


          output$discrete <- bar_plot

    }

    if(type %in% c("both","continuous")) {

      # need to change the posterior predict function to get the argument from ...

      cont_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL,ntrys=100)

      # if prediction was continuous, get continuous value from the outcome

      plot_data_dens <- lapply(unique(group), function(g) {

        counter <- 0

        apply(cont_dist[,group==g], 1, function(col) {

          counter <<- counter + 1

          tibble(pred=col[!(outcome %in% c(0,1))],
                          true=outcome[!(outcome %in% c(0,1))],
                          draw=counter,
                 group=g)

        }) %>% bind_rows

      }) %>% bind_rows


      cont_plot <- plot_data_dens %>%
        ggplot(aes(x=pred,group=draw)) +
        geom_density(alpha=0.7) +
        theme_minimal() +
        labs(y="Probability Density",x=paste0("Continuous values of ",names(model$data)[1]),
             caption="Each black line is a draw from the posterior distribution.\nData distribution represented by gray line.") +
        ggtitle("Posterior Predictions for Continuous Data",
                subtitle=paste0("Outcome: ",names(model$data)[1]))

      # use different syntax depending on version of ggplot2

      if(packageVersion('ggplot2')=="3.4.0") {

        cont_plot <- cont_plot + geom_density(aes(x=true),linewidth=2,colour="gray",alpha=0.7)

      } else {

        cont_plot <- cont_plot + geom_density(aes(x=true),size=2,colour="gray",alpha=0.7)

      }

      output$continuous <- cont_plot


    }

    return(output)

}

