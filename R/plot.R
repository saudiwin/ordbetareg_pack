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
#' @param dv If you fit a model with multiple DVs/responses,
#' pass the name of the DV as a character value.
#' Note: this must be the same as the name of the column
#' in the data used to fit the model.
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
                             dv=NULL,
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

    if(!is.null(dv)) {

      full_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL,
                                     resp=dv)
    } else {

      full_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL)

    }

    if(!is.null(dv)) {

      outcome <- model$data[[dv]]

    } else {

      outcome <- model$data[[1]]

    }

    if(!is.null(outcome_label)) {

      outcome_label <- outcome_label

    } else {

      if(!is.null(dv)) {

        outcome_label <- dv

      } else {

        outcome_label <- names(model$data)[1]
      }

    }

    this_star <- transformr::poly_star()

    if(reverse_bounds) {

        # revert to original scale

        if(!is.null(dv)) {

          l_bound <- model$lower_bound[[dv]]
          up_bound <- model$upper_bound[[dv]]

        } else {


          l_bound <- model$lower_bound
          up_bound <- model$upper_bound


        }




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

          # define a simple function to calculate quantiles

          collapse_func <- function(x) {


            tibble(ymin=quantile(x, .05),
                   ymax=quantile(x, .95),
                   y=median(x))

          }

          bar_plot <- ggplot(plot_data_bar,
                             aes(x=type,y=var)) +
            geom_col(data=true_data_bar, aes(y=true),width=.5,fill="gray") +
            stat_summary(fun.data=collapse_func,geom = "pointrange") +
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
      if(!is.null(dv)) {

        cont_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL,ntrys=100,
                                       resp=dv)
      } else {

        cont_dist <- posterior_predict(model,ndraws=ndraws,cores=NULL,ntrys=100)
      }


      if(reverse_bounds) {

        # revert to original scale

        if(!is.null(dv)) {

          l_bound <- model$lower_bound[[dv]]
          up_bound <- model$upper_bound[[dv]]

        } else {


          l_bound <- model$lower_bound
          up_bound <- model$upper_bound


        }


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

#' Theme for Heiss plot
#' @noRd
.theme_ongo <- function(base_size = 11, base_family = "", prior = FALSE) {
  ret <- theme_bw(base_size, base_family) +
    theme(panel.background = element_rect(fill = "#ffffff", colour = NA),
          title = element_text(size = rel(1), family = base_family, face = "bold"),
          plot.subtitle = element_text(size = rel(0.8),
                                       family = base_family, face = "plain"),
          plot.caption = element_text(margin = margin(t = 10), size = rel(0.75), hjust = 0,
                                      family = base_family, face = "plain"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.15),
          panel.spacing = unit(1, "lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = rel(0.8),
                                    family = base_family, face = "bold"),
          axis.title.x = element_text(hjust = 0, margin = margin(t = 10)),
          axis.title.y = element_text(hjust = 1, margin = margin(r = 10)),
          axis.title.y.right = element_text(hjust = 0, margin = margin(l = 10)),
          legend.position = "bottom",
          legend.title = element_text(size = rel(0.7), vjust = 0.5,
                                      family = base_family, face = "plain"),
          legend.key.size = unit(0.7, "line"),
          legend.key = element_blank(),
          legend.spacing = unit(0.1, "lines"),
          legend.justification = "left",
          legend.margin = margin(t = -5, b = 0, l = 0, r = 0),
          strip.text = element_text(size = rel(0.9), hjust = 0,
                                    family = base_family, face = "bold"),
          strip.background = element_rect(fill = "white", colour = NA))

  if (prior) {
    ret <- ret +
      theme(panel.grid.major = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            panel.border = element_blank())
  } else {
    ret
  }
}


#' Heiss Plot for Predicted Proportions of Bounded Scale Components
#'
#' The Heiss plot, developed by the statistician Andrew Heiss, is a plot of the predicted proportions of components on a bounded scale that are grouped by the unique levels of a grouping variable or factor (such as a random effect) in the model. The plot excels at showing how the scale components--that is, the bottom, middle continuous, and top ends of the scale--vary with a discrete variable while also capturing posterior uncertainty. This plot was the winner of the 2023 ordbetareg Visualization Prize.
#'
#' For more details of the plot, see:
#'
#' Heiss, Andrew and Ye, Meng. "Enforcing Boundaries: China’s Overseas NGO Law and Operational Constraints for Global Civil Society." Working Paper, 2023. <https://stats.andrewheiss.com/compassionate-clam/notebook/manuscript.html>.
#'
#' @param object A fitted [ordbetareg()] model object.
#' @param grouping_fac A character string indicating the name of the discrete column in the data used for grouping predictions. Must be a valid column name that was passed to [ordbetareg()].
#' @param recode_group_labels Optional. A character vector of new labels for the grouping factor levels. Must match the number and order of unique levels/values in `grouping_fac`.
#' @param ndraws Optional. The number of posterior draws to use for predictions. If `NULL`, all available draws are used.
#' @param show_category_perc_labels Logical. Whether to display category percentage labels on the plot. Defaults to `TRUE`.
#' @param category_label_font_size The `ggplot2` font size for the labels on the
#' scale components (if `show_category_perc_labels` is `TRUE`). Defaults to 3.
#' @param strip_text_font A `ggplot2::element_text` object defining the font style for facet strip text. Defaults to `element_text(face = "plain", size = 9)`.
#' @param plot_title Title of the plot. Defaults to "Predicted Proportions of Bounded Scale Components".
#' @param plot_subtitle Subtitle of the plot. Defaults to a message indicating the grouping variable.
#' @param plot_caption Caption text for the plot. Defaults to a detailed description of the plot contents.
#' @param plot_caption_width Width (in characters) at which the caption is wrapped. Defaults to 60.
#' @param calc_func A function used to calculate the central tendency of predictions. Defaults to `mean`.
#' @param lb Lower bound for uncertainty intervals. Defaults to 0.05 (5th percentile).
#' @param upb Upper bound for uncertainty intervals. Defaults to 0.95 (95th percentile).
#' @param plot_font_size Base font size for the plot. Defaults to 11.
#' @param plot_font Base font family for the plot. Defaults to an empty string (uses system default).
#' @param y_axis_label Label for the y-axis. Defaults to "Predicted Proportions".
#' @param legend_name Legend title. Defaults to "Scale Components".
#' @param component_colors A character vector of colors for the plot components (bottom, continuous, top). Defaults to `c("#ef8737", "#bb292c", "#62205f")`.
#' @param component_labels A character vector of labels for the scale/outcome components (bottom, continuous, top). Defaults to `c("0", "(0-1)", "1")`.
#' @param ... Additional arguments passed to [posterior_epred_ordbeta())].
#'
#' @return A `ggplot2` object representing the predicted proportions of the components.
#'
#' @examples
#' # Load a fitted model object and create a plot for
#' # distinct values of the factor education
#' #
#' # data('ord_fit_mean')
#' #
#' # plot_heiss(ord_fit_mean,ndraws=100)
#' #
#' # See introductory package vignette for more information on function options
#'
#' @import ggplot2
#' @importFrom dplyr select mutate pull group_by summarize n left_join
#' @importFrom stringr str_wrap str_extract
#' @importFrom tidyr gather unnest
#' @importFrom scales label_percent label_number
#' @importFrom tibble as_tibble
#' @importFrom stats quantile
#'
#' @export
plot_heiss <- function(object,
                       grouping_fac=NULL,
                       recode_group_labels=NULL,
                       ndraws=NULL,
                       show_category_perc_labels=TRUE,
                       category_label_font_size=3,
                       strip_text_font=element_text(face="plain",size=9),
                       plot_title="Predicted Proportions of Bounded Scale Components",
                       plot_subtitle=paste0("By Unique Values of ",grouping_fac),
                       plot_caption="Plot shows predicted proportions of the components of a bounded scale, i.e. the predicted (expected) probability of the top value of the scale, the intermediate continuous values, and the bottom value of the scale. The predictions are subset for unique values of a grouping factor. The predictions are shown for multiple posterior draws to indicate uncertainty. Labels on components indicate posterior quantiles for the probability of that component for each level of the grouping variable.",
                       plot_caption_width=70,
                       calc_func=mean,lb=.05,upb=.95,
                       plot_font_size=11,
                       plot_font="",
                       y_axis_label="Predicted Proportions",
                       legend_name="Scale Components",
                       component_colors=c("#ef8737", "#bb292c", "#62205f"),
                       component_labels=c("0","(0-1)","1"),
                       ...) {

  if(is.null(grouping_fac)) stop("To use this function, you must pass the name of a discrete column in the data to the grouping_fac option.")

  if(!(grouping_fac %in% names(object$data))) stop("Please pass the name of one of the columns int the model data, which are: ", paste0(names(object$data), collapse=","))

  # factor to group by

  grouping_fac_data <- select(object$data, grouping_fac=grouping_fac) %>%
    mutate(rownum=1:n())

  if(length(unique(pull(grouping_fac_data, grouping_fac)))>100) stop("The column you passed has more than 100 values. That will overload the plotting function. Totally wack, bro.")

  if(! any(class(pull(grouping_fac_data, grouping_fac)) %in% c("factor","ordered"))) {

    grouping_fac_data <- mutate(grouping_fac_data,
                                grouping_fac = factor(grouping_fac))

  }

  if(!is.null(recode_group_labels) && length(recode_group_labels) != levels(pull(grouping_fac_data, grouping_fac))) stop("If you pass a character vector of new grouping factor labels to recode_group_labels, this character vector must be of the same length as the number of levels/unique values in the grouping factor.")

  if(is.null(recode_group_labels)) recode_group_labels <- levels(pull(grouping_fac_data, grouping_fac))

  if(any("ordered" %in% class(pull(grouping_fac_data, grouping_fac)))) {

    grouping_fac_data <-     grouping_fac_data %>%
      mutate(grouping_fac = ordered(grouping_fac,
                                    labels=stringr::str_wrap(recode_group_labels,
                                                             width=10)))

  } else {

    grouping_fac_data <-     grouping_fac_data %>%
      mutate(grouping_fac = factor(grouping_fac,
                                    labels=stringr::str_wrap(recode_group_labels,
                                                             width=10)))

  }

  # predict output first

  preds_local1 <- object %>%
    posterior_epred_ordbeta(component="bottom",
                            ndraws=ndraws,...) %>%
    as_tibble(.name_repair="universal_quiet") %>%
    mutate(.draw=1:n(),
           component=component_labels[1]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  preds_local2 <- object %>%
    posterior_epred_ordbeta(component="continuous",
                            ndraws=ndraws,...) %>%
    as_tibble(.name_repair="universal_quiet") %>%
    mutate(.draw=1:n(),
           component=component_labels[2]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  preds_local3 <- object %>%
    posterior_epred_ordbeta(component="top",
                            ndraws=ndraws,...) %>%
    as_tibble(.name_repair="universal_quiet") %>%
    mutate(.draw=1:n(),
           component=component_labels[3]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  # combine into one big tibble

  preds_local <- bind_rows(preds_local1,
                           preds_local2,
                           preds_local3) %>%
    left_join(grouping_fac_data,by="rownum")


  preds_local_plot <- preds_local %>%
    group_by(.draw, component, grouping_fac) %>%
    summarize(y=calc_func(pred),
              ymin=quantile(pred, upb),
              ymax=quantile(pred, lb)) %>%
    mutate(component=factor(component, levels=rev(component_labels)))

  preds_local_text <- preds_local_plot %>%
    group_by(component, grouping_fac) %>%
    summarize(median_prop = tibble(y_agg=median(y),
                                   ymin=quantile(y, lb),
                                   ymax=quantile(y, upb))) %>%
    unnest(median_prop) %>%
    group_by(grouping_fac) %>%
    mutate(y_plot = (y_agg / 2) + lag(cumsum(y_agg), default = 0)) %>%
    mutate(y_plot = 1 - y_plot) %>%
    mutate(prop_nice = label_percent(accuracy = 1)(y_agg)) %>%
    mutate(prop_ci_nice = paste0(label_number(accuracy = 1, scale = 100)(ymin),
                                 "–",
                                 label_percent(accuracy = 1)(ymax)))

  # resize plot caption

  if(!is.null(plot_caption) || length(plot_caption)>1) plot_caption <- stringr::str_wrap(plot_caption,plot_caption_width)

  outplot <- preds_local_plot %>%
    ggplot(aes(x = .draw, y = y)) +
    geom_area(aes(fill = component), position = position_stack()) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
    scale_fill_manual(name=legend_name,
                      values = component_colors) +
    facet_wrap(~grouping_fac, strip.position = "bottom", nrow = 1) +
    labs(x = NULL, y = y_axis_label, fill = NULL,
         caption=plot_caption) +
    .theme_ongo(base_size=plot_font_size,
                base_family=plot_font) +
    theme(
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(5, "pt"),
      strip.text = strip_text_font
    ) +
    ggtitle(plot_title,
            subtitle=plot_subtitle)

  if(show_category_perc_labels) outplot <- outplot + geom_text(
    data = preds_local_text,
    aes(x = max(preds_local_plot$.draw)-min(quantile(preds_local_plot$.draw,.35)), y = y_plot, label = prop_ci_nice),
    size = category_label_font_size, fontface = "bold", color = "white"
  )

  outplot

}

