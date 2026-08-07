# make the Heiss plot

library(ordbetareg)
library(dplyr)
library(purrr)
library(tidybayes)
library(marginaleffects)
library(tidyr)
library(scales)
library(ggplot2)

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

data("pew")

# prepare data

model_data <- select(pew,therm,
                     education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                     income="F_INCOME_FINAL")

# It takes a while to fit the models. Run the code
# below if you want to load a saved fitted model from the
# package, otherwise use the model-fitting code

data("ord_fit_mean")

plot_heiss <- function(object,
                       grouping_fac=NULL,
                       recode_group_labels=NULL,
                       ndraws=NULL,
                       show_category_perc_labels=TRUE,
                       strip_text_font=element_text(face="plain",size=9),
                       plot_title="Predicted Proportions of Bounded Scale Components",
                       plot_subtitle=paste0("By Unique Values of ",grouping_fac),
                       plot_caption="Plot shows predicted proportions of the components of a bounded scale, i.e. the predicted (expected) probability of the top value of the scale, the intermediate continuous values, and the bottom value of the scale. The predictions are subset for unique values of a grouping factor. The predictions are shown for multiple posterior draws to indicate uncertainty.",
                       plot_caption_width=60,
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
      mutate(grouping_fac = ordered(grouping_fac,
                                    labels=stringr::str_wrap(recode_group_labels,
                                                             width=10)))

  }

  # predict output first

  preds_local1 <- object %>%
    posterior_epred_ordbeta(component="bottom",
                            ndraws=ndraws,...) %>%
    as_tibble %>%
    mutate(.draw=1:n(),
           component=component_labels[1]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  preds_local2 <- object %>%
    posterior_epred_ordbeta(component="continuous",
                            ndraws=ndraws,...) %>%
    as_tibble %>%
    mutate(.draw=1:n(),
           component=component_labels[2]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  preds_local3 <- object %>%
    posterior_epred_ordbeta(component="top",
                            ndraws=ndraws,...) %>%
    as_tibble %>%
    mutate(.draw=1:n(),
           component=component_labels[3]) %>%
    gather(key="rownum",value="pred",-component,-.draw) %>%
    mutate(rownum=as.numeric(stringr::str_extract(rownum, '[0-9]+')))

  # combine into one big tibble

  preds_local <- bind_rows(preds_local1,
                           preds_local2,
                           preds_local3) %>%
    left_join(grouping_fac_data)


  preds_local_plot <- preds_local %>%
    group_by(.draw, component, grouping_fac) %>%
    summarize(y=calc_func(pred),
              ymin=quantile(pred, upb),
              ymax=quantile(pred, lb)) %>%
    mutate(component=factor(component, levels=component_labels))

  preds_local_text <- preds_local_plot %>%
    group_by(component, grouping_fac) %>%
    summarize(median_prop = median_qi(y)) %>%
    unnest(median_prop) %>%
    group_by(grouping_fac) %>%
    mutate(y_plot = (y / 2) + lag(cumsum(y), default = 0)) %>%
    mutate(y_plot = 1 - y_plot) %>%
    mutate(prop_nice = label_percent(accuracy = 1)(y)) %>%
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
    size = 3, fontface = "bold", color = "white"
  )

  outplot

}

plot_heiss(ord_fit_mean, grouping_fac="education")





