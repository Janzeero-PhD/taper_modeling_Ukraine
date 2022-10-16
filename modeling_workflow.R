rm(list = ls(all.names = T))
path <- "insert your path to downloaded and unpacked folder DB\\DB\\" # download and unpack DB.zip with all data
setwd(path)
library(nlme)
library(dplyr)
library(ggplot2)
# -------------------------------------------------------------------------------------------

trsp <- "tico"
trsp.df <- read.csv(paste0(path, toupper(trsp), "/", trsp, "_taper.csv"))

# Bootstrap to evaluate parameters CV of the Kozak's (2004) taper function ---------------
mod <- nls_bootstrap(trsp.df, iter = 1000, random_seed = 10, save_table = F)
mod

# Mixed-effect modeling assuming that two parameters with higher CV have a random effect
mixmod <- nlme(dob ~ Kozak_2004(D = d, H = h, hi, a0, a1, a2, b1, b2, b3, b4, b5, b6), 
               data = trsp.df,  
               fixed = a0 + a1 + a2 + b1 + b2 + b3 + b4 + b5 + b6 ~ 1,
               random = b4+ b6 ~ 1 | Plot/MD, # select two most variable parameters from the mod object
               correlation = corCAR1(form = ~ 1 | Plot/MD), # dealing with correlation structure (continuous lags)
               start = c(a0 = mod[1,1], a1 = mod[1,2], a2 = mod[1,3], 
                         b1 = mod[1,4], b2 = mod[1,5], b3 = mod[1,6], b4 = mod[1,7], b5 = mod[1,8], b6 = mod[1,9]),  # Start values
               weights = varPower(form = ~ d),  # Weights to model heteroscedasticity
               method = "REML",  # Method used for fitting: ML or REML
               control = nlmeControl(tolerance = 1e-01, msMaxIter = 150))  # Tolerance reduced to get the convergence
summary(mixmod)
acf(residuals(mixmod, type = "normalized"))

# Plot selected model profile along with fitted fixed- and random-effect model
plot_single.fn(sample_index = 10, dataset = trsp.df, mixmod, save = F)

# Plot IQ range of random effects
plot_ranef.fn(D = 30, H = 24, mixmod, ran_coef1 = "b4", ran_coef2 = "b6", save = T)

# Overlay several species profiles
# pisy       quro       frex       cabe       potr       bepe       algl       tico
# "#e54c00"  "#7f3300"  "#2b7161"  "#63736F"  "#0d7b07"  "#00ccff"  "#6f248a"  "#ffe800"
#                       "#4169E1"                                              "#FFA500"
plot_several_sp.fn(D = 48, H = 30, 
                   sp_list = c("pisy", "quro", "frex", "cabe",  "potr", "bepe", "algl", "tico"), 
                   color_list = c("#e54c00", "#7f3300", "#4169E1", "#63736F", "#0d7b07", "#00ccff", "#6f248a", "#FFA500"), 
                   mixmod_list = list(mixmod_pisy, mixmod_quro, mixmod_frex, mixmod_cabe, mixmod_potr, mixmod_bepe, mixmod_algl, mixmod_tico),
                   save = T)

# Plot residuals
plot_residuals(trsp, trsp.df, mixmod, save_plot = T)

# Plot observed and predicted volumes
pred_vs_obs_volumes(trsp, mixmod, lab_step = 0.5, save_plot = T, save_table = F)

# Plot form factors
formfact_vs_Vtabl(trsp, xlims = c(8,52), ylims = c(0.3, 0.65), save_plot = T)

# Compare MAB and RMSE of volume estimation using taper function vs published volume equations
predV_vs_volumeEQ(trsp, mixmod, save_table = F)

# Data set statistics
dataset_stat(trsp.df)

# ------------------------------------- Functions -------------------------------------------
Kozak_2004 <- function(D, H, hi, a0, a1, a2, b1, b2, b3, b4, b5, b6) {
  # Notation from: Kozak, A. (2004). My last words on taper equations. 
  # The Forestry Chronicle, 80(4), 507–515. https://doi.org/10.5558/tfc80507-4
  p <- 1.3/H
  zi = hi/H
  Qi = 1.0 - zi^(1/3)
  Xi = Qi/(1-p^(1/3))
  a0*D^a1*H^a2*Xi^(b1*zi^4 + b2*(1/exp(D/H)) + b3*Xi^(0.1) + b4*(1/D) + b5*H^Qi + b6*Xi)
}

nls_bootstrap <- function(dataset, iter, random_seed, save_table){
  filename <- paste0(path, toupper(trsp), "/", trsp, "_niter_", iter, "_nls_stat.csv")
  loop <- function(df, df2, n_iter){
    for(i in 1:n_iter){
      sampled <- trsp.df[sample(1:nrow(dataset), replace = T),]
      mod <- try(nls(dob ~ Kozak_2004(D = d, H = h, hi, a0, a1, a2, b1, b2, b3, b4, b5, b6), 
                     data = sampled,
                     start = list(a0 = 0.9, a1 = 0.95, a2 = 0.046, 
                                  b1 = 0.36, b2 = -0.050, b3 = 0.52, b4 = -0.59, b5 = 0.02, b6 = 0.03)))
      if(class(mod) == "try-error") next
      df[i,] <- unname(coef(mod))
      df2[i,] <- unname(coef(summary(mod))[,4])
    }
    list(na.omit(df), na.omit(df2))
  }
  
  if(!file.exists(filename)) {
    print(paste("NLS stat", toupper(trsp), "will be caclulated"))
    coeff <- coeff.p <- array(dim = c(iter,9))
    set.seed(random_seed)
    coeff <- loop(coeff, coeff.p, iter)
    
    if(nrow(coeff[[1]]) < iter){
      
      while(nrow(coeff[[1]]) < iter){
        n_fail <- iter - nrow(coeff[[1]])
        coeff_ <- coeff_.p <- array(dim = c(n_fail,9))
        coeff_ <- loop(coeff_, coeff_.p, n_fail)
        coeff[[1]] <- rbind(coeff[[1]], coeff_[[1]])
        coeff[[2]] <- rbind(coeff[[2]], coeff_[[2]])
        n_fail <- iter - nrow(coeff[[1]])
      }
    }

    mn.coef <- apply(coeff[[1]], 2, mean, na.rm = T)
    sd.coef <- apply(coeff[[1]], 2, sd, na.rm = T)
    cv.coef <- abs(sd.coef/mn.coef*100)
    p.value <- apply(coeff[[2]], 2, function(x) sum(x < 0.05)/iter*100)
    stat <- data.frame(rbind(mn.coef, sd.coef, cv.coef, p.value))
    stat[5,] <- apply(coeff[[1]], 2, length)
    row.names(stat) <- c("Mean", "SD", "CV", "% p value < 0.05 ", "No of resamples")
    names(stat) <- c("a0", "a1", "a2", "b1", "b2", "b3", "b4", "b5", "b6")
    if(save_table == T){
      write.csv(stat, file = filename, row.names = T, quote = F)
      }
    } else {
    print(paste("The NLS stat in file", paste0(trsp, "_niter_", iter, "_nls_stat.csv"), "has been found"))
    stat <- read.csv(filename, row.names = "X") # reads the input data
    }
  stat
}

plot_single.fn <- function(sample_index, dataset, mixmod, save){
  sample_list <- sort(unique(dataset$MD))
  sampleID <- sample_list[sample_index]
  dataset <- dataset[dataset$MD == sampleID,]

  md = substitute(sampleID); V = dataset$Vob[1]; d = dataset$d[1]; h = dataset$h[1]
  f = dataset$fob[1]; q1 = dataset$q1[1]; q2 = dataset$q2[1]; q3 = dataset$q3[1]
  
  coef_fixed <- unname(mixmod$coefficients$fixed)
  a0 = coef_fixed[1]; a1 = coef_fixed[2]; a2 = coef_fixed[3]
  b1 = coef_fixed[4]; b2 = coef_fixed[5]; b3 = coef_fixed[6]
  b4 = coef_fixed[7]; b5 = coef_fixed[8]; b6 = coef_fixed[9]
  coef_random <- as.numeric(coef(mixmod)[sample_index,])
  a0_ = coef_random[1]; a1_ = coef_random[2]; a2_ = coef_random[3]
  b1_ = coef_random[4]; b2_ = coef_random[5]; b3_ = coef_random[6]
  b4_ = coef_random[7]; b5_ = coef_random[8]; b6_ = coef_random[9]
  
  p <- ggplot(data = dataset, aes(x = hi, y = dob)) + 
    geom_line(aes(colour = "black")) + geom_point(aes(colour = "black"), pch = 21, fill = 'white') +
    geom_line(aes(x = hi, y = Kozak_2004(d, h, hi, a0,a1,a2,b1,b2,b3,b4,b5,b6), size = 'fixed'), colour = "magenta") +
    geom_line(aes(x = hi, y = Kozak_2004(d, h, hi, a0_,a1_,a2_,b1_,b2_,b3_,b4_,b5_,b6_), size = "random"), colour = "blue")+
    scale_colour_manual(name = NULL, values = c("black" = "black"), breaks = c("black"),
                        labels = c("black" = "Stem outside bark"))  +
    scale_size_manual(name = NULL, values = c("fixed" = 0.5, "random" = 0.5),
                      breaks = c("fixed", "random"),
                      labels = c("fixed" = "Fixed-effect model",
                                 "random" = "Random-effect model")) +
    guides(size = guide_legend(override.aes = list(colour = list("magenta", "blue")))) +
    scale_x_continuous(breaks = c(0, 1, seq(3, dataset$h[1]+2, 2)),
                       labels = c(0, 1, seq(3, dataset$h[1]+2, 2)),
                       limits = c(0, dataset$h[1] + 1),  name = "Height, m") +
    scale_y_continuous(breaks = seq(0, max(Kozak_2004(D = d, H = h, hi = 0, a0,a1,a2,b1,b2,b3,b4,b5,b6), dataset[dataset$hi == 0, "dob"]) + 2, 5),
                       limits = c(0, max(Kozak_2004(D = d, H = h, hi = 0, a0,a1,a2,b1,b2,b3,b4,b5,b6), dataset[dataset$hi == 0, "dob"]) + 2),  name = "Diameter, cm") +
    labs(title = paste('Sample tree ', md), subtitle = bquote(italic('d = ') ~ .(d) ~ 'cm' ~
                                                           italic('    h = ') ~ .(h) ~ 'm' ~
                                                           italic('    V = ') ~ .(V) ~ 'm'**3 ~
                                                           italic('    f = ') ~ .(f) ~
                                                           italic('    q'[2]) ~ ' =' ~ .(q2))) +
    theme_bw() +
    theme(legend.position = c(1, 1), legend.justification=c(1, 1),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = 'white', colour = 'grey60'),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, "cm"),
          panel.grid.minor = element_blank())
  p
  if (save == T){
    ggsave(filename = paste0( "SampleTree_", sampleID, ".png"), 
           device = 'png', width = 14, height = 8, units = 'cm', dpi = 600) 
  } else p
}

plot_several_sp.fn <- function(D, H, sp_list, color_list, mixmod_list, save){
  coef_fixed <- unname(mixmod_list[[1]]$coefficients$fixed)
  p <- ggplot() +
    geom_line(aes_(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1),
                                                     a0 = coef_fixed[1], a1 = coef_fixed[2], a2 = coef_fixed[3],
                                                     b1 = coef_fixed[4], b2 = coef_fixed[5], b3 = coef_fixed[6],
                                                     b4 = coef_fixed[7], b5 = coef_fixed[8], b6 = coef_fixed[9]), colour = color_list[1]), size = 0.25) +
    scale_x_continuous(breaks = c(0, 1, seq(3, H+1, 2)),
                       labels = c(0, 1, seq(3, H+1, 2)),
                       limits = c(0, H+1),  name = "Height, m") +
    scale_y_continuous(breaks = seq(0, max(Kozak_2004(D, H, seq(0, H, 0.1),
                                                      a0 = coef_fixed[1], a1 = coef_fixed[2], a2 = coef_fixed[3],
                                                      b1 = coef_fixed[4], b2 = coef_fixed[5], b3 = coef_fixed[6],
                                                      b4 = coef_fixed[7], b5 = coef_fixed[8], b6 = coef_fixed[9])) + 5, 5),
                       limits = c(0, max(Kozak_2004(D, H, seq(0, H, 0.1),
                                                    a0 = coef_fixed[1], a1 = coef_fixed[2], a2 = coef_fixed[3],
                                                    b1 = coef_fixed[4], b2 = coef_fixed[5], b3 = coef_fixed[6],
                                                    b4 = coef_fixed[7], b5 = coef_fixed[8], b6 = coef_fixed[9])) + 5),  name = "Diameter, cm") +
    labs(subtitle = bquote(italic('D = ') ~ .(D) ~ 'cm' ~ italic('    H = ') ~ .(H) ~ 'm'))
  if(length(sp_list) > 1){
    for(i in 2:length(sp_list)){
      coef_fixed <- unname(mixmod_list[[i]]$coefficients$fixed)
      p <- p + geom_line(aes_(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1),
                                                                a0 = coef_fixed[1], a1 = coef_fixed[2], a2 = coef_fixed[3],
                                                                b1 = coef_fixed[4], b2 = coef_fixed[5], b3 = coef_fixed[6],
                                                                b4 = coef_fixed[7], b5 = coef_fixed[8], b6 = coef_fixed[9]), colour = color_list[i]), size = 0.25)
    }
  }
  p <- p + scale_colour_manual(name = NULL, values = color_list, breaks = color_list, labels = toupper(sp_list)) +
    theme_bw() +
    theme(legend.position = c(1, 1), legend.justification=c(1, 1),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = 'white', colour = 'grey60'),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, "cm"),
          panel.grid.minor = element_blank())
  p
  
  if(save == T){
    ggsave(filename = paste0(path, toupper(trsp), "/", trsp, "_profile_compare.png"),
           device = 'png', width = 12, height = 7, units = 'cm', dpi = 600)
    p
  } else p
}


plot_ranef.fn <- function(D, H, mixmod, ran_coef1 = "a2", ran_coef2 = "b6", save){
  coef_fixed <- unname(mixmod$coefficients$fixed)
  # quantiles of random-effect parameters 
  coef1_quantiles <- quantile(ranef(mixmod)[[2]][,1], probs = c(0.25, 0.75))
  coef2_quantiles <- quantile(ranef(mixmod)[[2]][,2], probs = c(0.25, 0.75))
  
  coef_list <- c("a0","a1","a2","b1","b2","b3","b4","b5","b6")
  ranef1_index <- which(coef_list == ran_coef1)
  ranef2_index <- which(coef_list == ran_coef2)
  
  min1_min2 <- max1_max2 <- min1_max2 <- max1_min2 <- rep(0,9)
  min1_min2[c(ranef1_index, ranef2_index)] <- c(coef1_quantiles[1], coef2_quantiles[1])
  max1_max2[c(ranef1_index, ranef2_index)] <- c(coef1_quantiles[2], coef2_quantiles[2])
  min1_max2[c(ranef1_index, ranef2_index)] <- c(coef1_quantiles[1], coef2_quantiles[2])
  max1_min2[c(ranef1_index, ranef2_index)] <- c(coef1_quantiles[2], coef2_quantiles[1])
  
  coef_rndm1 <- coef_fixed + min1_min2
  coef_rndm2 <- coef_fixed + max1_max2
  coef_rndm3 <- coef_fixed + min1_max2
  coef_rndm4 <- coef_fixed + max1_min2
  
  # how coefficients will be shown in Fig.
  ran_coef1 <- strsplit(ran_coef1, "")
  ran_coef1_term <- ran_coef1[[1]][1];  ran_coef1_index <- as.numeric(ran_coef1[[1]][2])
  ran_coef2 <- strsplit(ran_coef2, "")
  ran_coef2_term <- ran_coef2[[1]][1];  ran_coef2_index <- as.numeric(ran_coef2[[1]][2])

    p <- ggplot() + 
    geom_line(aes(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1), 
                                                     a0 = coef_rndm1[1], a1 = coef_rndm1[2], a2 = coef_rndm1[3],
                                                     b1 = coef_rndm1[4], b2 = coef_rndm1[5], b3 = coef_rndm1[6],
                                                     b4 = coef_rndm1[7], b5 = coef_rndm1[8], b6 = coef_rndm1[9]), colour = "black"), size = 0.25) +
    geom_line(aes(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1), 
                                                     a0 = coef_rndm2[1], a1 = coef_rndm2[2], a2 = coef_rndm2[3],
                                                     b1 = coef_rndm2[4], b2 = coef_rndm2[5], b3 = coef_rndm2[6],
                                                     b4 = coef_rndm2[7], b5 = coef_rndm2[8], b6 = coef_rndm2[9]), colour = "black"), size = 0.25) +
    geom_line(aes(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1), 
                                                     a0 = coef_rndm3[1], a1 = coef_rndm3[2], a2 = coef_rndm3[3],
                                                     b1 = coef_rndm3[4], b2 = coef_rndm3[5], b3 = coef_rndm3[6],
                                                     b4 = coef_rndm3[7], b5 = coef_rndm3[8], b6 = coef_rndm3[9]), colour = "black"), size = 0.25) +
    geom_line(aes(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1), 
                                                     a0 = coef_rndm4[1], a1 = coef_rndm4[2], a2 = coef_rndm4[3],
                                                     b1 = coef_rndm4[4], b2 = coef_rndm4[5], b3 = coef_rndm4[6],
                                                     b4 = coef_rndm4[7], b5 = coef_rndm4[8], b6 = coef_rndm4[9]), colour = "black"), size = 0.25) +
    geom_line(aes(x = seq(0, H, 0.1), y = Kozak_2004(D, H, seq(0, H, 0.1), 
                                                     a0 = coef_fixed[1], a1 = coef_fixed[2], a2 = coef_fixed[3],
                                                     b1 = coef_fixed[4], b2 = coef_fixed[5], b3 = coef_fixed[6],
                                                     b4 = coef_fixed[7], b5 = coef_fixed[8], b6 = coef_fixed[9]), colour = "magenta"), size = 0.25) +
    
    scale_colour_manual(name = NULL, values = c("magenta", "black"), breaks = c("magenta", "black"),
                        labels = c("magenta" = "Fixed-efect model", "black" = bquote("Model with random effects estimated for"~ 
                                                                                 .(ran_coef1_term)[.(ran_coef1_index)] ~ "and" ~
                                                                                 .(ran_coef2_term)[.(ran_coef2_index)])))  +
    scale_x_continuous(breaks = c(0, 1, seq(3, H+1, 2)),
                       labels = c(0, 1, seq(3, H+1, 2)),
                       limits = c(0, H+1),  name = "Height, m") +
    scale_y_continuous(breaks = seq(0, max(Kozak_2004(D, H, seq(0, H, 0.1), 
                                                      a0 = coef_rndm1[1], a1 = coef_rndm1[2], a2 = coef_rndm1[3],
                                                      b1 = coef_rndm1[4], b2 = coef_rndm1[5], b3 = coef_rndm1[6],
                                                      b4 = coef_rndm1[7], b5 = coef_rndm1[8], b6 = coef_rndm1[9])) + 5, 5),
                       limits = c(0, max(Kozak_2004(D, H, seq(0, H, 0.1), 
                                                    a0 = coef_rndm1[1], a1 = coef_rndm1[2], a2 = coef_rndm1[3],
                                                    b1 = coef_rndm1[4], b2 = coef_rndm1[5], b3 = coef_rndm1[6],
                                                    b4 = coef_rndm1[7], b5 = coef_rndm1[8], b6 = coef_rndm1[9])) + 7),  name = "Diameter, cm") +
    labs(title = toupper(trsp)) +
    theme_bw() +
    theme(legend.position = c(1, 1), legend.justification=c(1, 1),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = 'white', colour = 'grey60'),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.8, 'lines'),
          legend.spacing.y = unit(0, "cm"),
          panel.grid.minor = element_blank())
  p
  if (save == T){
    ggsave(filename = paste0(path, toupper(trsp), "/", trsp, "_ranef.png"), 
           device = 'png', width = 12, height = 7, units = 'cm', dpi = 600) 
    p
  } else p
}

plot_residuals <- function(trsp, dataset, mixmod, save_plot){
  coef_fixed <- unname(mixmod$coefficients$fixed)
  a0 = coef_fixed[1]; a1 = coef_fixed[2]; a2 = coef_fixed[3]
  b1 = coef_fixed[4]; b2 = coef_fixed[5]; b3 = coef_fixed[6]
  b4 = coef_fixed[7]; b5 = coef_fixed[8]; b6 = coef_fixed[9]
  
  dataset <- dataset %>% mutate(cuts_d = as.numeric(as.character(cut(d, seq(2, 102, 4), seq(4, 100, 4), right = F)))) %>%
    mutate(cuts_h = as.numeric(as.character(cut(h, seq(4.5, 43.5, 1), seq(5, 43, 1), right = F))))
  dataset.fit <- dataset %>% rowwise() %>% 
    mutate(dob_fit = round(Kozak_2004(D = d, H = h, hi = hi, a0,a1,a2,b1,b2,b3,b4,b5,b6), 1)) %>% 
    mutate(delta_dob = dob_fit - dob) %>%
    mutate(he = hi/h) %>%
    mutate(delta_de = delta_dob/dob*100)
  
  dataset.fit_mean <- dataset.fit %>% ungroup() %>% 
    mutate(cuts_he = as.numeric(as.character(cut(he, seq(-0.05, 1.05, 0.1), seq(0, 1, 0.1), right = F)))) %>%
    group_by(cuts_he)  %>% summarise(delta_dob = mean(delta_dob)) 
  
  p <- ggplot() + 
    geom_abline(intercept = 0, slope = 0) +
    geom_hline(yintercept = 0) +
    geom_point(data = dataset.fit, aes(x = he, y = delta_dob, fill = "#619CFF"), pch = 21, alpha = 0.5) +
    geom_point(data = dataset.fit_mean, aes(x = cuts_he, y = delta_dob, fill = "red"), size = 3, pch = 21) +
    scale_x_continuous(name = "Relative height", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "Rediduals, ñm", breaks = seq(-10,10,5), labels = seq(-10,10,5), limits = c(-12.5,12.5)) +
    labs(title = toupper(trsp)) +
    scale_fill_manual(name = NULL, values = c("#619CFF", "red"), breaks = c("#619CFF", "red"),
                      labels = c("#619CFF" = "Diameter residuals",
                                 "red" = "Mean residual values calculated for 0.1h bins")) +
    theme(panel.background = element_rect(colour = 'grey60', fill = NA),
          panel.grid.major = element_line(colour = 'grey80'),
          panel.grid.minor = element_line(colour = NA),
          axis.line = element_line(colour = 'grey20'),
          legend.position = "bottom",
          legend.key = element_rect(fill = 'white')) 
  p
  
  if (save_plot == T){
    ggsave(filename = paste0(path, toupper(trsp), "/", trsp, "_residuals.png"), device = 'png', 
           width = 12, height = 9, units = 'cm', dpi = 600) 
    p
  } else p
}

LogVol.Simps <- function(d, h, L, U, stump = 0.2, mixmod, k = 100){
  coef_fixed <- unname(mixmod$coefficients$fixed)
  a0 = coef_fixed[1]; a1 = coef_fixed[2]; a2 = coef_fixed[3]
  b1 = coef_fixed[4]; b2 = coef_fixed[5]; b3 = coef_fixed[6]
  b4 = coef_fixed[7]; b5 = coef_fixed[8]; b6 = coef_fixed[9]
  
  df <- data.frame()
  section = (U-L)/k
  # browser()
  if(L < stump){ # shifts the Lover end to the stump point if it is bellow
    h0 = stump
    ht = L + section
    hm = h0 + (ht-h0)/2
    g0 <- 0.7854*(Kozak_2004(D = d, H = h, hi = h0, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gm <- 0.7854*(Kozak_2004(D = d, H = h, hi = hm, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gt <- 0.7854*(Kozak_2004(D = d, H = h, hi = ht, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    df[1,1] <- round((ht-h0)/6*(g0+4*gm+gt),10)
   } else { 
    h0 = L
    ht = L + section
    hm = h0 + (ht-h0)/2
    g0 <- 0.7854*(Kozak_2004(D = d, H = h, hi = h0, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gm <- 0.7854*(Kozak_2004(D = d, H = h, hi = hm, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gt <- 0.7854*(Kozak_2004(D = d, H = h, hi = ht, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    df[1,1] <- round(section/6*(g0+4*gm+gt),10)
  }
  
  for (n in 2:k){
    h0 = L+(n-1)*section
    ht = L+n*section
    hm = h0 + (ht-h0)/2

    g0 <- 0.7854*(Kozak_2004(D = d, H = h, hi = h0, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gm <- 0.7854*(Kozak_2004(D = d, H = h, hi = hm, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    gt <- 0.7854*(Kozak_2004(D = d, H = h, hi = ht, a0,a1,a2,b1,b2,b3,b4,b5,b6))^2/10000
    df[n,1] <- round(section/6*(g0+4*gm+gt),10)
  }
  df
}

pred_vs_obs_volumes <- function(trsp, mixmod, lab_step, save_plot, save_table = T){
  filename <- paste0(path, toupper(trsp), "/", trsp, "_Vprd_vs_Vobs.csv")
  
  if (!file.exists(filename)) {
    print(paste("Stem volume for", toupper(trsp), "will be caclulated"))
    trsp.df <- trsp.df %>% filter(hi == 1.3) %>% rowwise() %>% 
      mutate(Vpred = sum(LogVol.Simps(d = d, h = h, L = 0.0, U = h, stump = 0.0, k = 100, 
                                      mixmod = mixmod))) %>%
      mutate(fpred = round(Vpred/d^2/h/0.7854*10000,3))
    if(save_table == T){
      write.csv(trsp.df, file = filename, row.names = F, quote = F)
    }
  } else {
    print(paste("The file ", paste0(trsp, "_Vprd_vs_Vobs.csv"), "has been found"))
    trsp.df <- read.csv(filename) # reads the input data
  }
  
  p <- ggplot() + 
    geom_point(data = trsp.df, aes(x = Vob, y = Vpred), fill = "#619CFF", pch = 21, alpha = 0.5) + 
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(breaks = seq(0, max(trsp.df$Vob, na.rm = T), lab_step), 
                       limits = (c(0, max(trsp.df$Vob, na.rm = T))),
                       name = bquote("Observed volume,"~"m"**3)) +
    scale_y_continuous(breaks = seq(0, max(trsp.df$Vob, na.rm = T), lab_step),
                       limits = (c(0, max(trsp.df$Vob, na.rm = T))),
                       name = bquote("Predicted volume,"~"m"**3)) +
    labs(title = toupper(trsp)) +
    
    theme(panel.background = element_rect(colour = 'grey60', fill = NA),
          panel.grid.major = element_line(colour = 'grey80'),
          panel.grid.minor = element_line(colour = NA),
          axis.line = element_line(colour = 'grey20'),
          legend.key = element_rect(fill = 'white')) +
    theme(legend.position = c(0, 1), legend.justification=c(0, 1),
          legend.background = element_rect(fill = 'white', colour = 'grey60'),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.8, 'lines'))
  p
  if (save_plot == T){
    ggsave(filename = paste0(path, toupper(trsp), "/", trsp, "_Vprd_vs_Vobs", '.png'), 
           device = 'png', width = 8, height = 8, units = 'cm', dpi = 300) 
  } else p
  p
}


formfact_vs_Vtabl <- function(trsp, xlims = c(8, 60), ylims = c(0.35, 0.7), save_plot = F){
  subpath <- paste0(toupper(trsp), "\\", trsp)
  mature_vt <-  read.csv(paste0(path, subpath, "_volume_tables.csv"))
  trsp.df <- trsp.df %>% mutate(cuts_d = as.numeric(as.character(cut(d, seq(2, 102, 4), seq(4, 100, 4), right = F)))) %>%
    mutate(cuts_h = as.numeric(as.character(cut(h, seq(3.5, 41.5, 1), seq(4, 41, 1), right = F))))
  
  mature_vt <- mature_vt %>% mutate(fob_vt2 = round(Vob/d^2/h/0.7854*10000,3)) 
  merged <- merge(x = trsp.df, y = mature_vt[,-3], by.x = c("cuts_d", "cuts_h"), by.y = c("d", "h"), all.x = T)
  merged <- data.frame(merged)
  # ------------
  # actual values of form factors
  trsp.group_d <- merged %>% group_by(cuts_d) %>% summarise(fob_avr = mean(fob)) 
  # predicted form factors using taper equations
  trsp.df_md <- merged %>% filter(hi == 1.3) %>% rowwise() %>% 
    mutate(Vpred = sum(LogVol.Simps(d = d, h = h, L = 0.0, U = h, stump = 0.0, k = 100, mixmod = mixmod))) %>%
    mutate(fpred = round(Vpred/d^2/h/0.7854*10000,3)) %>% ungroup()
  trsp.group.m_d <- trsp.df_md %>% group_by(cuts_d) %>% 
    summarise_at(c("fpred", "fob_vt2"), mean, na.rm = TRUE) 
  
  p <- ggplot() +
      geom_point(data = trsp.df_md, aes(x = d, y = fob, size = "1"), pch = 16, colour = "black", alpha = 0.5) +
      geom_point(data = trsp.group_d, aes(x = as.numeric(cuts_d), y = fob_avr, size = "2"), fill = alpha("white", 0.5), pch = 21) +
      geom_line(data = trsp.group.m_d, aes(x = as.numeric(cuts_d), y = fpred, colour = "3")) +
      geom_line(data = trsp.group.m_d, aes(x = as.numeric(cuts_d), y = fob_vt2, colour = "4")) +
      
      scale_size_manual(name = "Observed form factors", values = c("1" = 1, "2" = 3), 
                        breaks = c("1", "2"), 
                        labels = c("Individual trees",
                                   "Mean values for\n4-cm diameter bins")) +
      scale_colour_manual(name = "Predicted mean values\nfor 4-cm diameter bins", 
                          values = c("3" = "magenta", "4" = "blue"), 
                          breaks = c("3", "4"),
                          labels = c("Taper equation",
                                     "Volume equation")) +
      scale_x_continuous(name = "Diameter, cm", breaks = seq(4, 80, 4), labels = seq(4, 80, 4), limits = xlims) +
      scale_y_continuous(name = "Cylindrical form factor", breaks = seq(0.3, 0.7, 0.05), labels = seq(0.3, 0.7, 0.05), limits = ylims) +
      guides(size = guide_legend(override.aes = list(fill = list( "red", "white")), nrow = 2, title.position="top"),
             colour = guide_legend(nrow = 2, title.position="top")) +
      labs(title = toupper(trsp)) +
      theme(panel.background = element_rect(colour = 'grey60', fill = NA),
            panel.grid.major = element_line(colour = 'grey80'),
            panel.grid.minor = element_line(colour = NA),
            axis.line = element_line(colour = 'grey20'),
            legend.key = element_rect(fill = 'white'),
            legend.position = 'bottom') 
  p
  if (save_plot == T){
    ggsave(filename = paste0(path, toupper(trsp), "/", trsp, "_form_factors.png"), device = 'png', 
           width = 12, height = 11, units = 'cm', dpi = 600) 
  } else p
  p
}

predV_vs_volumeEQ <- function(trsp, mixmod, save_table = F){
  filename <- paste0(path, toupper(trsp), "/", trsp, "_Vprd_vs_Vobs.csv")
  
  if (!file.exists(filename)) {
    print(paste("Stem volume for", toupper(trsp), "will be caclulated"))
    trsp.df <- trsp.df %>% filter(hi == 1.3) %>% rowwise() %>% 
      mutate(Vpred = sum(LogVol.Simps(d = d, h = h, L = 0.0, U = h, stump = 0.0, k = 100, 
                                      mixmod = mixmod))) %>%
      mutate(fpred = round(Vpred/d^2/h/0.7854*10000,3))
    if(save_table == T){
      write.csv(trsp.df, file = filename, row.names = F, quote = F)
    }
  } else {
    print(paste("The file ", paste0(trsp, "_Vprd_vs_Vobs.csv"), "has been found"))
    trsp.df <- read.csv(filename) # reads the input data
  }
  
  if(trsp == "pisy") V = function(d,h){exp(7.767-0.04235*log(d+8)-0.6374*log(h+2)+0.02158*(h+2))*d^2*h*7.854*10^(-8)}
  if(trsp == "quro") V = function(d,h){(0.3812+0.4955*d^(-0.4811))*d^2*h*7.854*10^(-5)}
  if(trsp == "abal") V = function(d,h){(0.4118+1.176/(h-0.289))*d^2*h*7.854*10^(-5)}
  if(trsp == "cabe") V = function(d,h){(0.103+0.5889*d^(-0.1702))*d^2*h*7.854*10^(-5)}
  if(trsp == "piab") V = function(d,h){(438-2.3*h+0.0934*h^2-(d-40))/(0.163+0.00874*d)*d^2*h*7.854*10^(-8)}
  if(trsp == "potr") V = function(d,h){0.3081*d^1.8708*h^1.1932*10^(-4)}
  if(trsp == "bepe") V = function(d,h){0.5631*d^1.755*h^1.073*10^(-4)}
  if(trsp == "frex") V = function(d,h){(-0.1634+0.8110*d^(-0.07708))*d^2*h*7.854*10^(-5)}
  if(trsp == "fasy") V = function(d,h){(0.5007-0.001507*d+2.106/d^2)*d^2*h*7.854*10^(-5)}
  if(trsp == "algl") V = function(d,h){(0.4105+1.285/(d+1.084))*d^2*h*7.854*10^(-5)}
  if(trsp == "tico") V = function(d,h){(-4.166+4.849*d^(-0.01390))*d^2*h*7.854*10^(-5)}
  
  MAB <- function(obs, prd){
    resid <- abs(trsp.df[prd] - trsp.df[obs])
    n <- nrow(trsp.df)
    mab <- sum(resid)/n
    mab
  }
  
  RMSE <- function(obs, prd){
    resid2 <- (trsp.df[prd] - trsp.df[obs])^2
    n <- nrow(trsp.df)
    rmse <- (sum(resid2)/n)^0.5
    rmse
  }
  
  MPB <- function(obs, prd){
    resid <- abs(trsp.df[prd] - trsp.df[obs])
    n <- nrow(trsp.df)
    mpb <- sum(resid)/sum(trsp.df[obs])*100
    mpb
  }
    
  trsp.df$Veq <- V(trsp.df$d, trsp.df$h)
  out <- array(dim = c(5,3))
  out[1,1] <- mean(trsp.df$Vpred);  out[1,2] <- mean(trsp.df$Veq); out[1,3] <- mean(trsp.df$Vob)
  out[2,1] <- MAB("Vob", "Vpred");  out[2,2] <- MAB("Vob", "Veq")
  out[3,1] <- RMSE("Vob", "Vpred"); out[3,2] <- RMSE("Vob", "Veq")
  out[4,1] <- MPB("Vob", "Vpred");  out[4,2] <- MPB("Vob", "Veq")
  out[5,] <- rep(nrow(trsp.df),3)
  # write.csv(trsp.df, file = filename, row.names = F, quote = F)
  out <- data.frame(out)
  names(out) <- c("Taper FN", "Volume EQ", "Observed")
  row.names(out) <- c("Mean", "MAB", "RMSE", "MPB", "N trees")
  out
}

dataset_stat <- function(dataset){
  min <- apply(dataset[,c("Age", "d", "h", "Vob", "fob")], 2, min, na.rm = T)
  max <- apply(dataset[,c("Age", "d", "h", "Vob", "fob")], 2, max, na.rm = T)
  mean <- apply(dataset[,c("Age", "d", "h", "Vob", "fob")], 2, mean, na.rm = T)
  out <- rbind(mean, min, max)
  names(out) <- c("Age", "d", "h", "Vob", "fob")
  list(out, "Number of sample plots:", length(unique(dataset$Plot)),
       "Number of sample trees:", length(unique(dataset$MD)))
} 
  
