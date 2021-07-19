Phage abundance profiles
================
Sam Barnett
19 July, 2021

-   [Introduction](#introduction)
-   [Phage abundance](#phage-abundance)
-   [Streptomyces abundances](#streptomyces-abundances)
-   [phage:bacterial ratio](#phagebacterial-ratio)
-   [More info about the Streptomyces OTUs](#more-info-about-the-streptomyces-otus)

Introduction
------------

From the metagenomic-SIP study we identified a vOTU that likely infects a Streptomyces host. To better understand the relationship between this phage and its possible hosts, we used qPCR to track the population dynamics of the two. We used a primer-probe designed from the contig sequence to track vOTU count in microcosms over time following the carbon pulse. We then used a Streptomyces specific rpoB gene primer set to both track rpoB counts over time and to track counts of specific Streptomyces OTUs using amplicon sequencing.

### Initiation

``` r
library(dplyr)
library(ggplot2)
library(ggtext)
library(knitr)

library(phyloseq)
library(grid)
library(gridExtra)

library(readxl)

source("/Users/sambarnett/Documents/Misc_code/paul_tol_colors.R")


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
```

Phage abundance
---------------

We measured phage abundance using a probe-primer set designed from the vOTU metagenome sequence and running qPCR.

``` r
phage_plate1.df = read_xlsx("/Users/sambarnett/Documents/Buckley Lab/Phage_project/PCR_results/Plate1_Run1_Phage/Plate1_run1_results.xlsx")
phage_plate2.df = read_xlsx("/Users/sambarnett/Documents/Buckley Lab/Phage_project/PCR_results/Plate2_Run1_Phage/Plate2_run1_results.xlsx")
```

### Getting count values based on standards

Plate 1

``` r
## Standard curve
phage_plate1_std.df = phage_plate1.df %>%
  filter(grepl("Std", Content)) %>%
  mutate(control_conc = as.numeric(Sample),
         Cq = as.numeric(Cq))

phage_plate1.lm = lm(Cq~log10(control_conc), data=phage_plate1_std.df)
phage_plate1.lm
```

    ## 
    ## Call:
    ## lm(formula = Cq ~ log10(control_conc), data = phage_plate1_std.df)
    ## 
    ## Coefficients:
    ##         (Intercept)  log10(control_conc)  
    ##              20.913               -3.586

``` r
phage_plate1_slope = phage_plate1.lm$coefficients[2]
phage_plate1_intercept = phage_plate1.lm$coefficients[1]

ggplot(data=phage_plate1_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point() +
  geom_abline(slope = phage_plate1_slope, intercept = phage_plate1_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
## Get concentrations
phage_plate1_unknowns.df = phage_plate1.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-phage_plate1_intercept)/phage_plate1_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc,
         Replicate = as.numeric(Replicate))

ggplot(data=phage_plate1_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point(size=3, shape=1) +
  geom_point(data=phage_plate1_unknowns.df, aes(x=log10(mean_conc), y=mean_Cq), size=3, color="red") +
  geom_abline(slope = phage_plate1_slope, intercept = phage_plate1_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-3-2.png)

Plate 2

``` r
## Standard curve
phage_plate2_std.df = phage_plate2.df %>%
  filter(grepl("Std", Content)) %>%
  mutate(control_conc = as.numeric(Sample),
         Cq = as.numeric(Cq))

phage_plate2.lm = lm(Cq~log10(control_conc), data=phage_plate2_std.df)
phage_plate2.lm
```

    ## 
    ## Call:
    ## lm(formula = Cq ~ log10(control_conc), data = phage_plate2_std.df)
    ## 
    ## Coefficients:
    ##         (Intercept)  log10(control_conc)  
    ##              20.911               -3.579

``` r
phage_plate2_slope = phage_plate2.lm$coefficients[2]
phage_plate2_intercept = phage_plate2.lm$coefficients[1]

ggplot(data=phage_plate2_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point() +
  geom_abline(slope = phage_plate2_slope, intercept = phage_plate2_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
## Get concentrations
phage_plate2_unknowns.df = phage_plate2.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-phage_plate2_intercept)/phage_plate2_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc,
         Replicate = as.numeric(Replicate))

ggplot(data=phage_plate2_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point(size=3, shape=1) +
  geom_point(data=phage_plate2_unknowns.df, aes(x=log10(mean_conc), y=mean_Cq), size=3, color="red") +
  geom_abline(slope = phage_plate2_slope, intercept = phage_plate2_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-4-2.png)

### Examining vOTU abundance over time

First by concentration

``` r
phage_unknowns.df = rbind(mutate(phage_plate1_unknowns.df, plate = "plate1"),
                          mutate(phage_plate2_unknowns.df, plate = "plate2")) %>%
  filter(!(Sample %in% c("Kit", "WCR1", "WCR2", "WCR3")))

phage_unknowns.sum = phage_unknowns.df %>%
  group_by(Day) %>%
  summarize(sum_mean_log_conc = mean(mean_log_conc),
            sd_log_conc = sd(mean_log_conc),
            n_samples = n()) %>%
  as.data.frame() %>%
  mutate(SE_log_conc = sd_log_conc/sqrt(n_samples)) %>%
  rename(mean_log_conc = sum_mean_log_conc)

ggplot(data=phage_unknowns.df, aes(x=Day, y=mean_log_conc)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_log_conc-sd_log_conc, ymax = mean_log_conc+sd_log_conc)) +
  geom_line(data=phage_unknowns.sum) +
  labs(x="Time since substrate addition (Day)", y="Log10 phage concentration (pM)")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
ggplot(data=phage_unknowns.sum, aes(x=Day, y=mean_log_conc)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_log_conc-SE_log_conc, ymax = mean_log_conc+SE_log_conc)) +
  geom_line() +
  labs(x="Time since substrate addition (Day)", y="Log10 phage concentration (pM)")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-5-2.png)

``` r
ggplot(data=phage_unknowns.sum, aes(x=Day, y=10^mean_log_conc)) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean_log_conc-SE_log_conc, ymax = mean_log_conc+SE_log_conc)) +
  geom_line() +
  labs(x="Time since substrate addition (Day)", y="Phage concentration (pM)")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-5-3.png)

Next by estimated phage counts

``` r
phage_plate1_unknowns.df = phage_plate1.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-phage_plate1_intercept)/phage_plate1_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc,
         Replicate = as.numeric(Replicate)) %>%
  mutate(count = mean_conc*1E-9*100*4*6.02214076E14,
         plate = "Plate1")

phage_plate2_unknowns.df = phage_plate2.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-phage_plate2_intercept)/phage_plate2_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc,
         Replicate = as.numeric(Replicate)) %>%
  mutate(count = mean_conc*1E-9*100*4*6.02214076E14,
         plate = "Plate2")

phage_plate_unknowns.df = rbind(phage_plate1_unknowns.df, phage_plate2_unknowns.df) %>%
  filter(!(Sample %in% c("Kit", "WCR1", "WCR2", "WCR3")))

phage_plate_unknowns.sum = phage_plate_unknowns.df %>%
  filter(Sample != "Kit") %>%
  group_by(Day) %>%
  summarize(mean_count = mean(count),
            sd_count = sd(count),
            n_samples = n()) %>%
  as.data.frame() %>%
  mutate(SE_count = sd_count/sqrt(n_samples))

phage_count.plot = ggplot(data=phage_plate_unknowns.sum, aes(x=Day, y=mean_count/1E6)) +
  geom_errorbar(aes(ymin = mean_count/1E6-SE_count/1E6, ymax = mean_count/1E6+SE_count/1E6), width=0, size=0.5) +
  geom_line(size=1) +
  geom_point(size=2, fill="white", shape=21) +
  labs(x="Time since substrate addition (Day)", y="vOTU count per g soil\n(x1,000,000)") +
  theme_bw()
phage_count.plot
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-6-1.png)

Measure linear growth and decay

``` r
# Growth between day 1 and 6
phage_days_1_6.df = phage_plate_unknowns.df %>%
  mutate(ln_count = log(count)) %>%
  filter(Day <= 6)
phage_days_1_6.model = lm(ln_count ~ Day, data=phage_days_1_6.df)
summary(phage_days_1_6.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_count ~ Day, data = phage_days_1_6.df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3561 -1.0903 -0.7893  1.7274  2.5595 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  10.2466     0.7704   13.30 6.03e-09 ***
    ## Day           1.1352     0.1967    5.77 6.49e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.566 on 13 degrees of freedom
    ## Multiple R-squared:  0.7192, Adjusted R-squared:  0.6976 
    ## F-statistic: 33.29 on 1 and 13 DF,  p-value: 6.492e-05

``` r
# Decay between day 6 and 48
phage_days_6_48.df = phage_plate_unknowns.df %>%
  mutate(ln_count = log(count)) %>%
  filter(Day >= 6)
phage_days_6_48.model = lm(ln_count ~ Day, data=phage_days_6_48.df)
summary(phage_days_6_48.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_count ~ Day, data = phage_days_6_48.df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.84715 -0.43154  0.08868  0.32795  1.02407 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 16.368910   0.218717  74.841  < 2e-16 ***
    ## Day         -0.025567   0.007463  -3.426  0.00301 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5368 on 18 degrees of freedom
    ## Multiple R-squared:  0.3947, Adjusted R-squared:  0.3611 
    ## F-statistic: 11.74 on 1 and 18 DF,  p-value: 0.003013

Streptomyces abundances
-----------------------

We measured Streptomyces abundance using a Streptomyces specific rpoB gene primer set. To remove non-Streptomyces OTUs and because we do not know which Strepotmyces species is the host, we also sequenced the rpoB amplicons to get abundances of each Streptomyces OTU.

### Getting count values based on standards

``` r
rpoB_plate1.df = read_xlsx("/Users/sambarnett/Documents/Buckley Lab/Phage_project/PCR_results/rpoB_plate1_Run2_txtFiles/rpoB_plate1_run2_results.xlsx")

## Standard curve (remove standards 7 and 8 due to primer dimers and very little amplification if any)
rpoB_plate1_std.df = rpoB_plate1.df %>%
  filter(grepl("Std", Content)) %>%
  filter(!(Content %in% c("Std-01", "Std-07", "Std-08"))) %>%
  mutate(control_conc = as.numeric(Sample),
         Cq = as.numeric(Cq))

rpoB_plate1.lm = lm(Cq~log10(control_conc), data=rpoB_plate1_std.df)
rpoB_plate1.lm
```

    ## 
    ## Call:
    ## lm(formula = Cq ~ log10(control_conc), data = rpoB_plate1_std.df)
    ## 
    ## Coefficients:
    ##         (Intercept)  log10(control_conc)  
    ##              10.290               -3.763

``` r
rpoB_plate1_slope = rpoB_plate1.lm$coefficients[2]
rpoB_plate1_intercept = rpoB_plate1.lm$coefficients[1]

ggplot(data=rpoB_plate1_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point() +
  geom_abline(slope = rpoB_plate1_slope, intercept = rpoB_plate1_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
## Get concentrations
rpoB_plate1_unknowns.df = rpoB_plate1.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-rpoB_plate1_intercept)/rpoB_plate1_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc)

ggplot(data=rpoB_plate1_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point(size=3, shape=1) +
  geom_point(data=rpoB_plate1_unknowns.df, aes(x=log10(mean_conc), y=mean_Cq), size=3, color="red") +
  geom_abline(slope = rpoB_plate1_slope, intercept = rpoB_plate1_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-8-2.png)

``` r
rpoB_plate2.df = read_xlsx("/Users/sambarnett/Documents/Buckley Lab/Phage_project/PCR_results/rpoB_Plate2_Run2_txtFiles/rpoB_Plate2_run2_results.xlsx")

## Standard curve
rpoB_plate2_std.df = rpoB_plate2.df %>%
  filter(grepl("Std", Content)) %>%
  filter(!(Content %in% c("Std-01", "Std-07", "Std-08"))) %>%
  mutate(control_conc = as.numeric(Sample),
         Cq = as.numeric(Cq))

rpoB_plate2.lm = lm(Cq~log10(control_conc), data=rpoB_plate2_std.df)
rpoB_plate2.lm
```

    ## 
    ## Call:
    ## lm(formula = Cq ~ log10(control_conc), data = rpoB_plate2_std.df)
    ## 
    ## Coefficients:
    ##         (Intercept)  log10(control_conc)  
    ##              10.035               -3.791

``` r
rpoB_plate2_slope = rpoB_plate2.lm$coefficients[2]
rpoB_plate2_intercept = rpoB_plate2.lm$coefficients[1]


ggplot(data=rpoB_plate2_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point() +
  geom_abline(slope = rpoB_plate2_slope, intercept = rpoB_plate2_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
## Get concentrations
rpoB_plate2_unknowns.df = rpoB_plate2.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-rpoB_plate2_intercept)/rpoB_plate2_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc)

ggplot(data=rpoB_plate2_std.df, aes(x=log10(control_conc), y=Cq)) +
  geom_point(size=3, shape=1) +
  geom_point(data=rpoB_plate2_unknowns.df, aes(x=log10(mean_conc), y=mean_Cq), size=3, color="red") +
  geom_abline(slope = rpoB_plate2_slope, intercept = rpoB_plate2_intercept, color="red")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-9-2.png)

### Examining raw rpoB abundance over time

``` r
rpoB_unknowns.df = rbind(mutate(rpoB_plate1_unknowns.df, plate = "plate1"),
                         mutate(rpoB_plate2_unknowns.df, plate = "plate2")) %>%
  filter(Sample != "Kit") %>%
  mutate(Treatment = ifelse(Sample %in% c("WCR1", "WCR2", "WCR3"), "Water only", "Substrate added"),
         Day = ifelse(Sample %in% c("WCR1", "WCR2", "WCR3"), 48, Day),
         Replicate = as.numeric(Replicate))

rpoB_unknowns.sum = rpoB_unknowns.df %>%
  group_by(Day, Treatment) %>%
  summarize(sum_mean_log_conc = mean(mean_log_conc),
            sd_log_conc = sd(mean_log_conc),
            n_samples = n()) %>%
  as.data.frame() %>%
  mutate(SE_log_conc = sd_log_conc/sqrt(n_samples)) %>%
  rename(mean_log_conc = sum_mean_log_conc)

ggplot(data=rpoB_unknowns.df, aes(x=Day, y=mean_log_conc, color=Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_log_conc-sd_log_conc, ymax = mean_log_conc+sd_log_conc)) +
  geom_line(data=rpoB_unknowns.sum) +
  labs(x="Time since substrate addition (Day)", y="Log10 rpoB concentration")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
ggplot(data=rpoB_unknowns.sum, aes(x=Day, y=mean_log_conc, color=Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_log_conc-SE_log_conc, ymax = mean_log_conc+SE_log_conc)) +
  geom_line() +
  labs(x="Time since substrate addition (Day)", y="Log10 rpoB concentration")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-10-2.png)

``` r
ggplot(data=rpoB_unknowns.sum, aes(x=Day, y=10^mean_log_conc, color=Treatment)) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean_log_conc-SE_log_conc, ymax = mean_log_conc+SE_log_conc)) +
  geom_line() +
  labs(x="Time since substrate addition (Day)", y="rpoB concentration")
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-10-3.png)

Next by estimated rpoB counts

``` r
# 25.666 ng/µl initial concentration (Janani's measurement)
# 35.761974 ng/µl initial concentration (Sam's measurement)
# 9854443 bp length
# 650 ng/nmol bp

genome_equivalents = ((35.76/(650*9854443))*1E-9)*6.02214076E23

rpoB_plate1_unknowns.df = rpoB_plate1.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-rpoB_plate1_intercept)/rpoB_plate1_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc) %>%
  mutate(count = mean_conc*genome_equivalents*100*4,
         plate = "Plate1")

rpoB_plate2_unknowns.df = rpoB_plate2.df %>%
  filter(grepl("Unkn", Content)) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  mutate(log_conc = (Cq-rpoB_plate2_intercept)/rpoB_plate2_slope) %>%
  group_by(Sample, Content) %>%
  summarize(mean_log_conc = mean(log_conc),
            sd_log_conc = sd(log_conc),
            mean_Cq = mean(Cq)) %>%
  as.data.frame %>%
  tidyr::separate(Sample, into=c("Day", "Replicate"), sep="R", remove=FALSE) %>%
  mutate(Day = as.numeric(gsub("D", "", Day)),
         mean_conc = 10^mean_log_conc) %>%
  mutate(count = mean_conc*genome_equivalents*100*4,
         plate = "Plate2")

rpoB_plate_unknowns.df = rbind(rpoB_plate1_unknowns.df, rpoB_plate2_unknowns.df) %>%
  filter(Sample != "Kit") %>%
  mutate(Treatment = ifelse(Sample %in% c("WCR1", "WCR2", "WCR3"), "Water only", "Substrate added"),
         Day = ifelse(Sample %in% c("WCR1", "WCR2", "WCR3"), 48, Day),
         Replicate = as.numeric(Replicate))

rpoB_plate_unknowns.sum = rpoB_plate_unknowns.df %>%
  group_by(Day, Treatment) %>%
  summarize(mean_count = mean(count),
            sd_count = sd(count),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_count = sd_count/sqrt(n_samples))

rpoB_count.plot = ggplot(data=rpoB_plate_unknowns.sum, aes(x=Day, y=mean_count/1E6, color=Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_count/1E6-SE_count/1E6, ymax = mean_count/1E6+SE_count/1E6)) +
  geom_line() +
  scale_color_manual(values= c("Substrate added" = "black", "Water only" = "grey")) +
  labs(x="Time since substrate addition (Day)", y="rpoB count per g soil\n(x1,000,000)") +
  lims(y=c(0, 5.5)) +
  theme_bw()
rpoB_count.plot
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-11-1.png)

### Streptomyces OTU abundances

Reading in the sequncing data

``` r
# The OTU table and metadata
rpob_OTU.biom = import_biom("/Users/sambarnett/Documents/Buckley Lab/Phage_project/final.otutab.biom")
rpob_metadata.df = read.table("/Users/sambarnett/Documents/Buckley Lab/Phage_project/rpoB_metadata.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  left_join(rpoB_plate_unknowns.df, by = c("Sample", "Day", "Replicate"))
rownames(rpob_metadata.df) = rpob_metadata.df$Sample

rpob_OTU_adj.df = data.frame(rpob_OTU.biom) %>%
  tibble::rownames_to_column(var="OTU") %>%
  tidyr::gather(key="Sample", value="read_count", -OTU) %>%
  group_by(Sample) %>%
  mutate(total_reads = sum(read_count)) %>%
  ungroup %>%
  mutate(rel_abund = read_count/total_reads) %>%
  left_join(select(rpob_metadata.df, Sample, count), by = "Sample") %>%
  mutate(abs_abund = rel_abund*count) %>%
  select(OTU, Sample, rel_abund, total_reads, abs_abund)

# The OTU taxonomy table
rpob_tax.df = read.table("/Users/sambarnett/Documents/Buckley Lab/Phage_project/final.otu.fasta.blca.out", header=FALSE, sep="\t") %>%
  rename(OTU = V1, Taxonomy = V2) %>%
  tidyr::separate(Taxonomy, into=c("Domain", "DQ", "Phylum", "PQ", "Class", "CQ", "Order", "OQ", "Family", "FQ", "Genus", "GQ", "Species", "SQ"), sep=";", convert=TRUE) %>%
  mutate(Domain = ifelse(DQ < 80 | is.na(DQ), "Unclassified", gsub("superkingdom:", "", Domain)),
         Phylum = ifelse(PQ < 80 | is.na(PQ), "Unclassified", gsub("phylum:", "", Phylum)),
         Class = ifelse(CQ < 80 | is.na(CQ), "Unclassified", gsub("class:", "", Class)),
         Order = ifelse(OQ < 80 | is.na(OQ), "Unclassified", gsub("order:", "", Order)),
         Family = ifelse(FQ < 80 | is.na(FQ), "Unclassified", gsub("family:", "", Family)),
         Genus = ifelse(GQ < 80 | is.na(GQ), "Unclassified", gsub("genus:", "", Genus)),
         Species = ifelse(SQ < 80 | is.na(SQ), "Unclassified", gsub("species:", "", Species))) %>%
  group_by(OTU) %>%
  ungroup %>%
  select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species)

print(paste("There are", nrow(filter(rpob_tax.df, Genus == "Streptomyces")), "Streptomyces OTUs"))
```

    ## [1] "There are 500 Streptomyces OTUs"

``` r
kable(rpob_OTU_adj.df %>%
        left_join(rpob_tax.df, by="OTU") %>%
        filter(rel_abund > 0) %>%
        mutate(Streptomyces = ifelse(Genus == "Streptomyces", 1, 0)) %>%
        group_by(Sample, total_reads) %>%
        summarize(n_OTU = n(), n_streptomyces_OTUs = sum(Streptomyces)) %>%
        ungroup())
```

| Sample |  total\_reads|  n\_OTU|  n\_streptomyces\_OTUs|
|:-------|-------------:|-------:|----------------------:|
| D14R1  |          2469|     332|                    170|
| D14R2  |          4309|     362|                    174|
| D14R3  |          3132|     282|                    128|
| D14R4  |          1867|     218|                     94|
| D14R5  |          4492|     481|                    197|
| D1R1   |          2486|     383|                    192|
| D1R2   |          2689|     383|                    177|
| D1R3   |          2155|     344|                    206|
| D1R4   |          1436|     263|                    152|
| D1R5   |          2954|     413|                    222|
| D30R1  |          3252|     446|                    183|
| D30R2  |          3179|     457|                    169|
| D30R3  |          3017|     416|                    238|
| D30R4  |          2128|     295|                    148|
| D30R5  |          1732|     272|                    122|
| D3R1   |          1769|     223|                    134|
| D3R2   |          5342|     469|                    264|
| D3R3   |          3440|     328|                    200|
| D3R4   |           796|     175|                    125|
| D3R5   |          2373|     245|                    132|
| D48R1  |             1|       1|                      1|
| D48R2  |          1472|     272|                    118|
| D48R3  |           683|     237|                     99|
| D48R4  |          1345|     310|                    174|
| D48R5  |          1662|     310|                    132|
| D6R1   |          2300|     279|                    141|
| D6R2   |          1914|     250|                    128|
| D6R3   |          8847|     465|                    222|
| D6R4   |          3658|     339|                    170|
| D6R5   |          3737|     303|                    144|
| WCR1   |           947|     240|                    107|
| WCR2   |           832|     262|                     89|
| WCR3   |           598|     229|                     87|

Tracking abundance over time

``` r
# First, plot all streptomyces OTUs combined
streptomyces_count.sum = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample, Day, Treatment) %>%
  summarize(sum_abs_abund = sum(abs_abund)) %>%
  ungroup %>%
  group_by(Day, Treatment) %>%
  summarize(mean_sum_abs_abund = mean(sum_abs_abund),
            sd_sum_abs_abund = sd(sum_abs_abund),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_sum_abs_abund = sd_sum_abs_abund/sqrt(n_samples))

streptomyces_count.plot = ggplot(data=streptomyces_count.sum, aes(x=Day, y=mean_sum_abs_abund/1E6)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = mean_sum_abs_abund/1E6-SE_sum_abs_abund/1E6, ymax = mean_sum_abs_abund/1E6+SE_sum_abs_abund/1E6), width=0, size=0.5) +
  geom_point(size=2, fill="white", shape=21) +
  #scale_shape_manual(values= c("Substrate added" = 21, "Water only" = 22)) +
  labs(x="Time since substrate addition (Day)", y="rpoB count per g soil\n(x1,000,000)") +
  theme_bw()

# Next, plot the top 10 most abundant stretpomyces OTUs separately
streptomyces_rank.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(OTU) %>%
  summarize(max_abs_abund = max(abs_abund)) %>%
  ungroup %>%
  arrange(-max_abs_abund) %>%
  mutate(Rank = row_number())


top_streptomyces_count.sum = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU %in% filter(streptomyces_rank.df, Rank <11)$OTU) %>%
  group_by(OTU, Day, Treatment) %>%
  summarize(mean_abs_abund = mean(abs_abund),
            sd_abs_abund = sd(abs_abund),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_abs_abund = sd_abs_abund/sqrt(n_samples),
         OTU = factor(OTU, levels=filter(streptomyces_rank.df, Rank <11)$OTU))

topOTU.col = paultol_colors(10)
names(topOTU.col) = filter(streptomyces_rank.df, Rank <11)$OTU

top_streptomyces_count.plot = ggplot(data=top_streptomyces_count.sum, 
                                     aes(x=Day, y=mean_abs_abund/1E6, color=OTU, fill=OTU)) +
  geom_line(color="black", size=1) +
  geom_line(size=0.75) +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.50, color="black") +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.45) +
  geom_point(size=2, color="black", shape=21) +
  scale_color_manual(values= topOTU.col) +
  scale_fill_manual(values= topOTU.col) +
  #scale_shape_manual(values= c("Substrate added" = 21, "Water only" = 22)) +
  labs(x="Time since substrate addition (Day)", y="rpoB count per g soil\n(x1,000,000)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title = element_text(size=8)) +
  guides(fill=guide_legend(ncol=5, title.position="top", title.hjust = 0.5),
         color=guide_legend(ncol=5, title.position="top", title.hjust = 0.5))

no1_top_streptomyces_count.plot = ggplot(data=filter(top_streptomyces_count.sum, OTU != "OTU_1"), 
                                         aes(x=Day, y=mean_abs_abund/1E6, color=OTU, fill=OTU)) +
  geom_line(color="black", size=1) +
  geom_line(size=0.75) +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.50, color="black") +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.45) +
  geom_point(size=2, color="black", shape=21) +
  scale_color_manual(values= topOTU.col) +
  scale_fill_manual(values= topOTU.col) +
  #scale_shape_manual(values= c("Substrate added" = 21, "Water only" = 22)) +
  labs(x="Time since substrate addition (Day)", y="rpoB count per g soil\n(x1,000,000)") +
  theme_bw() +
  guides(fill=guide_legend(ncol=1, title.position="top", title.hjust = 0.5),
         color=guide_legend(ncol=1, title.position="top", title.hjust = 0.5))

facet_top_streptomyces_count.plot = ggplot(data=top_streptomyces_count.sum, 
                                           aes(x=Day, y=mean_abs_abund/1E6, color=OTU, fill=OTU)) +
  geom_line(color="black", size=1) +
  geom_line(size=0.75) +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.50, color="black") +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.45) +
  geom_point(size=2, color="black", shape=21) +
  scale_color_manual(values= topOTU.col) +
  scale_fill_manual(values= topOTU.col) +
  #scale_shape_manual(values= c("Substrate added" = 21, "Water only" = 22)) +
  labs(x="Time since substrate addition (Day)", y="rpoB count per g soil\n(x1,000,000)") +
  facet_wrap(~OTU, ncol=1, strip.position="right", scales = "free_y") +
  theme_bw() +
  theme(strip.text = element_text(size=7),
        legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title = element_text(size=8)) +
  guides(fill=guide_legend(ncol=5, title.position="top", title.hjust = 0.5),
         color=guide_legend(ncol=5, title.position="top", title.hjust = 0.5))
```

Plotting the Streptomyces OTU counts

``` r
rpoB.leg = g_legend(top_streptomyces_count.plot)

comb.theme = theme(axis.text = element_text(size=7),
                   axis.title = element_blank(),
                   legend.position = "none")

comb.plot1 = cowplot::plot_grid(phage_count.plot + lims(y=c(0,NA)) + comb.theme + annotate("richtext", x=35, y=12, label="vOTU", size=9*5/14, hjust=0.5), 
                                streptomyces_count.plot + lims(y=c(0,NA)) + comb.theme + annotate("richtext", x=35, y=3.8, label="Combined<br>*Streptomyces* OTU<sub>*rpoB*</sub>", size=9*5/14, hjust=0.5), 
                                top_streptomyces_count.plot + lims(y=c(0,NA)) + comb.theme + annotate("richtext", x=35, y=2.05, label="Top 10<br>*Streptomyces* OTU<sub>*rpoB*</sub>", size=9*5/14, hjust=0.5),
                                ncol=1, align="v", labels=c("A", "B", "C"), label_size = 9)
comb.plot2 = cowplot::plot_grid(comb.plot1,
                                facet_top_streptomyces_count.plot + lims(y=c(0,NA)) + comb.theme,
                                ncol=2, labels=c("", "D"), label_size = 9)

x.grob <- textGrob("Time since C addition (days)", gp=gpar(fontsize=8))
y.grob <- textGrob("Counts per g soil (x1,000,000)", gp=gpar(fontsize=8), rot=90)

FigS1.plot = grid.arrange(arrangeGrob(comb.plot2, bottom=x.grob, left=y.grob))
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
#Figs2.plot = cowplot::plot_grid(comb.plot3, rpoB.leg, nrow=2, rel_heights = c(1,0.15))
FigS1.plot
```

    ## TableGrob (1 x 1) "arrange": 1 grobs
    ##   z     cells    name            grob
    ## 1 1 (1-1,1-1) arrange gtable[arrange]

``` r
#ggsave(FigS1.plot, filename = "/Users/sambarnett/Documents/Buckley Lab/Phage_project/Manuscript/Figures/FigS1.tiff", 
#       device = "tiff", width = 7, height = 7, units = "in")
```

What are the top 10 most abundant Streptomyces?

``` r
kable(streptomyces_rank.df %>%
        left_join(rpob_tax.df, by="OTU") %>%
        filter(Rank <= 10))
```

| OTU      |  max\_abs\_abund|  Rank| Domain   | Phylum         | Class          | Order            | Family            | Genus        | Species      |
|:---------|----------------:|-----:|:---------|:---------------|:---------------|:-----------------|:------------------|:-------------|:-------------|
| OTU\_1   |        3081196.1|     1| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_6   |         407323.4|     2| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_32  |         250183.4|     3| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | NRRL WC-3626 |
| OTU\_25  |         231352.4|     4| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_58  |         184265.3|     5| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_13  |         161011.4|     6| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_37  |         161011.4|     7| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_2   |         156523.8|     8| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |
| OTU\_113 |         144581.6|     9| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | collinus     |
| OTU\_139 |         143783.2|    10| Bacteria | Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |

For publication I want just the count figures of the phage and OTU\_1. I'll add on phage:bacterium ratio later.

``` r
OTU1_count.plot = ggplot(data=filter(top_streptomyces_count.sum, OTU == "OTU_1"), 
                         aes(x=Day, y=mean_abs_abund/1E6)) +
  geom_line(color="black", size=1) +
  geom_errorbar(aes(ymin = mean_abs_abund/1E6-SE_abs_abund/1E6, ymax = mean_abs_abund/1E6+SE_abs_abund/1E6), width = 0, size=0.50, color="black") +
  geom_point(size=2, color="black", shape=21, fill="white") +
  lims(y=c(0, NA)) +
  labs(x="Time since substrate addition (Day)", y=expression("OTU"["rpoB"]~1)) +
  theme_bw() +
  guides(fill=guide_legend(ncol=1, title.position="top", title.hjust = 0.5),
         color=guide_legend(ncol=1, title.position="top", title.hjust = 0.5))

comb.theme = theme(axis.text = element_text(size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size=8),
                   legend.position = "none")

phage_OTU1.plot = cowplot::plot_grid(phage_count.plot + comb.theme + labs(y="vOTU"), 
                                     OTU1_count.plot + comb.theme, 
                                     ncol=1, align="v", labels=c("A", "B"), label_size = 9)

x.grob <- textGrob("Time since C addition (Days)", gp=gpar(fontsize=8))
y.grob <- textGrob("qPCR based count per g soil (x1,000,000)", gp=gpar(fontsize=8), rot=90)

phage_OTU1.plot = grid.arrange(arrangeGrob(phage_OTU1.plot, bottom=x.grob, left=y.grob))
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-16-1.png)

### Measure linear growth and decay

First combined Streptomyces OTUs

``` r
# Growth between day 1 and 6
streptomyces_days_1_3.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample, Day, Treatment) %>%
  summarize(sum_abs_abund = sum(abs_abund)) %>%
  ungroup %>%
  mutate(ln_sum_abs_abund = log(sum_abs_abund)) %>%
  filter(Day <= 3)
streptomyces_days_1_3.model = lm(ln_sum_abs_abund ~ Day, data=streptomyces_days_1_3.df)
summary(streptomyces_days_1_3.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_sum_abs_abund ~ Day, data = streptomyces_days_1_3.df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5437 -0.3346 -0.2209  0.2664  0.7826 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  14.3121     0.3600  39.756 1.66e-09 ***
    ## Day           0.2167     0.1687   1.285     0.24    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5029 on 7 degrees of freedom
    ## Multiple R-squared:  0.1908, Adjusted R-squared:  0.07525 
    ## F-statistic: 1.651 on 1 and 7 DF,  p-value: 0.2397

``` r
# Decay between day 6 and 48
streptomyces_days_3_48.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample, Day, Treatment) %>%
  summarize(sum_abs_abund = sum(abs_abund)) %>%
  ungroup %>%
  mutate(ln_sum_abs_abund = log(sum_abs_abund)) %>%
  filter(Day >= 3)
streptomyces_days_3_48.model = lm(ln_sum_abs_abund ~ Day, data=streptomyces_days_3_48.df)
summary(streptomyces_days_3_48.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_sum_abs_abund ~ Day, data = streptomyces_days_3_48.df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.53421 -0.20364 -0.03996  0.24134  0.60479 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 14.972856   0.106537 140.542   <2e-16 ***
    ## Day         -0.006714   0.004450  -1.509    0.147    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3183 on 20 degrees of freedom
    ## Multiple R-squared:  0.1022, Adjusted R-squared:  0.05728 
    ## F-statistic: 2.276 on 1 and 20 DF,  p-value: 0.147

Next OTU\_1

``` r
# Growth between day 1 and 6
rpoB_OTU_1_days_1_6.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU == "OTU_1") %>%
  mutate(ln_abs_abund = log(abs_abund)) %>%
  filter(Day <= 6)
rpoB_OTU_1_days_1_6.model = lm(ln_abs_abund ~ Day, data=rpoB_OTU_1_days_1_6.df)
summary(rpoB_OTU_1_days_1_6.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_abs_abund ~ Day, data = rpoB_OTU_1_days_1_6.df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7914 -0.4025 -0.1247  0.4034  1.2252 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.68058    0.29875  42.446  1.9e-14 ***
    ## Day          0.33055    0.07519   4.396 0.000871 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5978 on 12 degrees of freedom
    ## Multiple R-squared:  0.6169, Adjusted R-squared:  0.585 
    ## F-statistic: 19.33 on 1 and 12 DF,  p-value: 0.0008713

``` r
# Growth between day 6 and 48
rpoB_OTU_1_days_6_48.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU == "OTU_1") %>%
  mutate(ln_abs_abund = log(abs_abund)) %>%
  filter(Day >= 6)
rpoB_OTU_1_days_6_48.model = lm(ln_abs_abund ~ Day, data=rpoB_OTU_1_days_6_48.df)
summary(rpoB_OTU_1_days_6_48.model)
```

    ## 
    ## Call:
    ## lm(formula = ln_abs_abund ~ Day, data = rpoB_OTU_1_days_6_48.df)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.76749 -0.20819  0.07891  0.27695  0.53266 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 14.685676   0.155059  94.710  < 2e-16 ***
    ## Day         -0.023916   0.005867  -4.076 0.000879 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3686 on 16 degrees of freedom
    ## Multiple R-squared:  0.5094, Adjusted R-squared:  0.4788 
    ## F-statistic: 16.61 on 1 and 16 DF,  p-value: 0.0008795

``` r
rpob_OTU_adj.sum = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU == "OTU_1") %>%
  group_by(Day) %>%
  summarize(mean_abs_abund = mean(abs_abund)) %>%
  ungroup

# How much of a decay is there between days 6 and 48 (%)?
(filter(rpob_OTU_adj.sum, Day==6)$mean_abs_abund - filter(rpob_OTU_adj.sum, Day==48)$mean_abs_abund)/filter(rpob_OTU_adj.sum, Day==6)$mean_abs_abund
```

    ## [1] 0.5656289

phage:bacterial ratio
---------------------

Often a phage-host system has a phage:bacterium ratio of ~10. Lets see if that is the case here and if it changes over time.

``` r
# First, all streptomyces OTUs combined
streptomyces_phage_ratio.sum = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample, Day, Treatment, Replicate) %>%
  summarize(Streptomyces_count = sum(abs_abund)) %>%
  ungroup %>%
  left_join(phage_plate_unknowns.df %>%
              select(Sample, Day, Replicate, count) %>%
              rename(vOTU_count = count),
            by = c("Sample", "Day", "Replicate")) %>%
  mutate(PBratio = vOTU_count/Streptomyces_count) %>%
  group_by(Day) %>%
  summarize(mean_PBratio = mean(PBratio),
            sd_PBratio = sd(PBratio),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_PBratio = sd_PBratio/sqrt(n_samples))

streptomyces_phage_ratio.plot = ggplot(data=streptomyces_phage_ratio.sum, aes(x=Day, y=mean_PBratio)) +
  geom_line(color="black", size=1) +
  geom_errorbar(aes(ymin = mean_PBratio-SE_PBratio, ymax = mean_PBratio+SE_PBratio), width = 0, size=0.50, color="black") +
  geom_point(size=2, color="black", shape=21, fill="white") +
  labs(x="Time since substrate addition (Day)", y="vOTU:Streptomyces rpoB count ratio") +
  theme_bw()


# Next, just the top 10 most abundant stretpomyces OTUs separately
top_streptomyces_phage_ratio.sum = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU %in% filter(streptomyces_rank.df, Rank <11)$OTU) %>%
  select(OTU, Sample, Day, Treatment, Replicate, abs_abund) %>%
  rename(OTU_count = abs_abund) %>%
  left_join(phage_plate_unknowns.df %>%
              select(Sample, Day, Replicate, count) %>%
              rename(vOTU_count = count),
            by = c("Sample", "Day", "Replicate")) %>%
  mutate(PBratio = vOTU_count/OTU_count) %>%
  group_by(OTU, Day) %>%
  summarize(mean_PBratio = mean(PBratio),
            sd_PBratio = sd(PBratio),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_PBratio = sd_PBratio/sqrt(n_samples))

OTU_1_phage_ratio.plot = ggplot(data=filter(top_streptomyces_phage_ratio.sum, OTU == "OTU_1"), aes(x=Day, y=mean_PBratio)) +
  geom_line(color="black", size=1) +
  geom_errorbar(aes(ymin = mean_PBratio-SE_PBratio, ymax = mean_PBratio+SE_PBratio), width = 0, size=0.50, color="black") +
  geom_point(size=2, color="black", shape=21, fill="white") +
  labs(x="Time since substrate addition (Day)", y="vOTU:OTU_1 rpoB count ratio") +
  theme_bw()

# Plot together
shared_theme = theme(axis.text = element_text(size=7),
                     axis.title = element_blank(),
                     plot.title = element_text(size=8, hjust=0.5))

comb_ratio.plot = cowplot::plot_grid(streptomyces_phage_ratio.plot + ggtitle("Combined Streptomyces OTUs") + shared_theme, 
                                     OTU_1_phage_ratio.plot + ggtitle("Streptomyces OTU_1") + shared_theme,
                                     ncol=1, labels=c("A", "B"), label_size = 9)

x.grob <- textGrob("Time since substrate addition (Day)", gp=gpar(fontsize=8))
y.grob <- textGrob("vOTU : bacterial OTU ratio", gp=gpar(fontsize=8), rot=90)

comb_ratio.plot = grid.arrange(arrangeGrob(comb_ratio.plot, bottom=x.grob, left=y.grob))
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-19-1.png)

Now combine with abundance figure

``` r
comb.theme = theme(axis.text = element_text(size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_markdown(size=8),
                   legend.position = "none")

comb_phage_OTU1_ratio.plot = cowplot::plot_grid(phage_count.plot + lims(y=c(0,NA)) + comb.theme + labs(y="vOTU count<br>per g soil (x 1,000,000)"), 
                                                OTU1_count.plot + lims(y=c(0,NA)) + comb.theme + labs(y="OTU<sub>*rpoB*</sub> 1 count<br>per gram soil (x 1,000,000)"), 
                                                OTU_1_phage_ratio.plot + lims(y=c(0,NA)) + comb.theme + labs(y="vOTU : OTU<sub>*rpoB*</sub> 1"),
                                                ncol=1, align="v", labels=c("A", "B", "C"), 
                                                label_size = 9)

x.grob <- textGrob("Time since C addition (days)", gp=gpar(fontsize=8))

Fig3.plot = grid.arrange(arrangeGrob(comb_phage_OTU1_ratio.plot, bottom=x.grob))
```

![](Phage_abundance_dynamics_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
#ggsave(Fig3.plot, filename = "/Users/sambarnett/Documents/Buckley Lab/Phage_project/Manuscript/Figures/Fig3.tiff", 
#       device = "tiff", width = 3.14961, height = 5, units = "in")
```

Does this change over time

``` r
OTU1_phage_ratio.df = rpob_OTU_adj.df %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by = "OTU") %>%
  filter(OTU == "OTU_1") %>%
  select(OTU, Sample, Day, Treatment, Replicate, abs_abund) %>%
  rename(OTU_count = abs_abund) %>%
  left_join(phage_plate_unknowns.df %>%
              select(Sample, Day, Replicate, count) %>%
              rename(vOTU_count = count),
            by = c("Sample", "Day", "Replicate")) %>%
  mutate(PBratio = vOTU_count/OTU_count)

# Change between day 1 and 6
OTU1_phage_ratio_days_1_6.model = lm(PBratio ~ Day, data=filter(OTU1_phage_ratio.df, Day <= 6))
summary(OTU1_phage_ratio_days_1_6.model)
```

    ## 
    ## Call:
    ## lm(formula = PBratio ~ Day, data = filter(OTU1_phage_ratio.df, 
    ##     Day <= 6))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8428 -0.7191 -0.3809  0.5570  2.7292 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -0.3588     0.8018  -0.447    0.663    
    ## Day           1.1623     0.2018   5.759 9.03e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.605 on 12 degrees of freedom
    ## Multiple R-squared:  0.7343, Adjusted R-squared:  0.7122 
    ## F-statistic: 33.17 on 1 and 12 DF,  p-value: 9.033e-05

More info about the Streptomyces OTUs
-------------------------------------

Relative abundance of OTU\_1 in the Streptomyces population

``` r
OTU1_abd.df = data.frame(rpob_OTU.biom) %>%
  tibble::rownames_to_column(var="OTU") %>%
  tidyr::gather(key="Sample", value="read_count", -OTU) %>%
  group_by(Sample) %>%
  mutate(total_reads = sum(read_count)) %>%
  ungroup %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by="OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample) %>%
  mutate(total_streptomyces_reads = sum(read_count)) %>%
  ungroup %>%
  mutate(rel_streptomyces_abund = read_count/total_streptomyces_reads) %>%
  filter(OTU == "OTU_1") %>%
  group_by(OTU, Day, Phylum, Class, Order, Family, Genus, Species) %>%
  summarize(mean_rel_streptomyces_abund = mean(rel_streptomyces_abund),
            sd_rel_streptomyces_abund = sd(rel_streptomyces_abund),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_rel_streptomyces_abund = sd_rel_streptomyces_abund/sqrt(n_samples))
kable(OTU1_abd.df)
```

| OTU    |  Day| Phylum         | Class          | Order            | Family            | Genus        | Species      |  mean\_rel\_streptomyces\_abund|  sd\_rel\_streptomyces\_abund|  n\_samples|  SE\_rel\_streptomyces\_abund|
|:-------|----:|:---------------|:---------------|:-----------------|:------------------|:-------------|:-------------|-------------------------------:|-----------------------------:|-----------:|-----------------------------:|
| OTU\_1 |    1| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.1649412|                     0.0157184|           5|                     0.0070295|
| OTU\_1 |    3| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.5048249|                     0.0640213|           4|                     0.0320106|
| OTU\_1 |    6| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.6644323|                     0.0595561|           5|                     0.0266343|
| OTU\_1 |   14| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.6240577|                     0.1563392|           5|                     0.0699170|
| OTU\_1 |   30| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.4896220|                     0.1725022|           5|                     0.0771453|
| OTU\_1 |   48| Actinobacteria | Actinobacteria | Streptomycetales | Streptomycetaceae | Streptomyces | Unclassified |                       0.3981013|                     0.2487297|           3|                     0.1436042|

Relative abundance of Streptomyces in all rpoB

``` r
streptomyces_abd.df = data.frame(rpob_OTU.biom) %>%
  tibble::rownames_to_column(var="OTU") %>%
  tidyr::gather(key="Sample", value="read_count", -OTU) %>%
  group_by(Sample) %>%
  mutate(total_reads = sum(read_count)) %>%
  ungroup %>%
  filter(total_reads > 1000) %>%
  left_join(rpob_metadata.df, by = "Sample") %>%
  left_join(rpob_tax.df, by="OTU") %>%
  filter(Genus == "Streptomyces") %>%
  group_by(Sample, Day, total_reads) %>%
  summarize(total_streptomyces_reads = sum(read_count)) %>%
  ungroup %>%
  mutate(rel_abund = total_streptomyces_reads/total_reads*100) %>%
  group_by(Day) %>%
  summarize(mean_rel_abund = mean(rel_abund),
            sd_rel_abund = sd(rel_abund),
            n_samples = n()) %>%
  ungroup %>%
  mutate(SE_rel_abund = sd_rel_abund/sqrt(n_samples))
kable(streptomyces_abd.df)
```

|  Day|  mean\_rel\_abund|  sd\_rel\_abund|  n\_samples|  SE\_rel\_abund|
|----:|-----------------:|---------------:|-----------:|---------------:|
|    1|          75.75958|        3.692654|           5|        1.651405|
|    3|          87.08687|        2.420148|           4|        1.210074|
|    6|          87.42275|        2.886824|           5|        1.291027|
|   14|          84.29753|        5.042983|           5|        2.255291|
|   30|          76.24935|        5.147348|           5|        2.301964|
|   48|          73.05891|        2.738966|           3|        1.581343|
