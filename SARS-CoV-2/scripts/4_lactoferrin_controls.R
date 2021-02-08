


library(pzfx)
library(ggplot2)

moi <- readr::read_csv(file="raw_data/lactoferrin_moi_raw_200515.csv") %>%
    dplyr::filter(MOI <= 10)
moi_treatment <- moi %>%
    dplyr::group_by(Condition, Cell_Line, MOI) %>%
    dplyr::summarize(
        mean = mean(Percent_Infection_Well_Level),
        sem =  mean(Percent_Infection_Well_Level) / dplyr::n())

## Bar plot version ##
p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position=c(.8, .8)) +
    ggplot2::geom_bar(
        data=moi_treatment,
        mapping=ggplot2::aes(
            x=MOI,
            y=mean,
            fill=Condition),
        size=1.3,
        position=ggplot2::position_dodge(width=1),
        stat="identity") +
    ggplot2::geom_errorbar(
        data=moi_treatment,
        mapping=ggplot2::aes(
            x=MOI,
            ymin=mean-sem,
            ymax=mean+sem,
            color=Condition),
        size=1.3,
        position=ggplot2::position_dodge(width=1),
        stat="identity") +
    ggplot2::geom_jitter(
        data=moi,
        mapping=ggplot2::aes(
            x=MOI,
            y=Percent_Infection_Well_Level,
            color=Condition,
            group=MOI),
        size=1.5,
        position=ggplot2::position_dodge(width=1),
        stat="identity") +
    ggplot2::scale_x_discrete("Multiplicity of infection") +
    ggplot2::scale_y_continuous("% cells infectected", limits=c(0,50)) +
    ggplot2::facet_wrap(~Cell_Line)

p %>% ggplot2::ggsave(
    filename="product/figures/lactoferrin/moi_bar_200515.pdf",
    width=7, height=4)     

### Line version ###
p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position=c(.8, .8)) +
    ggplot2::geom_line(
        data=moi_treatment,
        mapping=ggplot2::aes(
            x=log10(MOI+1),
            y=mean,
            color=Condition),
        size=1.3,
        stat="identity") +
    ggplot2::geom_jitter(
        data=moi,
        mapping=ggplot2::aes(
            x=log10(MOI+1),
            y=Percent_Infection_Well_Level,
            color=Condition,
            group=MOI),
        size=1.5,
        width=.02,
        stat="identity") +
    ggplot2::scale_x_continuous(
        "Multiplicity of infection",
        breaks=log10(1+c(0,1,5,10,50,100)),
        labels=c(0,1,5,10,50,100)) +
    ggplot2::scale_y_continuous("% cells infectected", limits=c(0,50)) +
    ggplot2::facet_wrap(~Cell_Line)

p %>% ggplot2::ggsave(
    filename="product/figures/lactoferrin/moi_line_200515.pdf",
    width=7, height=4)     


#######################################
transferrin <- pzfx::read_pzfx(path="raw_data/lactoferrin_apo_holo_transferrin_200515.pzfx") %>%
    dplyr::mutate(trial=1) %>%
    tidyr::pivot_longer(
        cols=-trial,
        names_to="treatment",
        values_to="percent_infection") %>%
    tidyr::separate(col=treatment, into=c("treatment", "value_type"), sep="_") %>%
    dplyr::mutate(value_type = value_type %>% stringr::str_to_lower()) %>%
    tidyr::pivot_wider(names_from="value_type", values_from="percent_infection")


###########################
rtPCR <- pzfx::read_pzfx(path="raw_data/lactoferrin_rtPCR_200515.pzfx") %>%
    dplyr::rename(hololactoferrin_dose = `[Hololactoferrin] (ug/mL)`) %>%
    tidyr::pivot_longer(
        cols=-hololactoferrin_dose,
        names_to="replicate") %>%
    dplyr::mutate(replicate = replicate %>% stringr::str_extract("[0-9]$") %>% as.numeric())

rtPCR_treatment <- rtPCR %>%
    dplyr::group_by(hololactoferrin_dose) %>%
    dplyr::summarize(mean = mean(value))


p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        data=m,
        mapping=ggplot2::aes(
            x=log10(MOI+1),
            y=mean,
            color=Condition),
        size=1.3,
        stat="identity") +
    ggplot2::geom_jitter(
        data=moi,
        mapping=ggplot2::aes(
            x=log10(MOI+1),
            y=Percent_Infection_Well_Level,
            color=Condition,
            group=MOI),
        size=1.5,
        width=.02,
        stat="identity") +
    ggplot2::scale_x_continuous(
        "Multiplicity of infection",
        breaks=log10(1+c(0,1,5,10,50,100)),
        labels=c(0,1,5,10,50,100)) +
    ggplot2::scale_y_continuous("% cells infectected", limits=c(0,50)) +
    ggplot2::facet_wrap(~Cell_Line)

p %>% ggplot2::ggsave(
    filename="product/figures/lactoferrin/moi_line_200515.pdf",
    width=7, height=4)     

