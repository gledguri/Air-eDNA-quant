# Libraries -----------------------------------------------------------------------------------

library(ggridges)
library(rstan)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(here)
library(MASS)
select <- dplyr::select

# Functions -----------------------------------------------------------------------------------
scientific <- function(x) {
  ifelse(
    x == 0, 
    expression(0), 
    str2expression(paste0(ifelse(x < 0, "-", ""), "10^", round(log10(abs(x)), 5)))
  )
}

scientific_10 <- function(x) {
  c <- log10(x)
  expo <- floor(c)
  base <- round(x / 10^expo,3)
  formatted <- paste0(base,'%*%',"10^", expo)
  str2expression(formatted)
}

round_by10 <- function(x) {
  c <- log10(x)
  expo <- floor(c)
  base <- round(x / 10^expo)
  return(base*10^expo)
}

extract_param <- function (model, par) 
{
  fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, 
                                                                   par = par)
  fit <- fit$summary
  return(fit %>% unlist() %>% as.data.frame %>% round(., 6))
}

join_ext_param <- function(stanmod,par){
  extract_param(stanmod,par) %>% 
    rownames_to_column('param') %>% 
    mutate(g_idx = str_extract(param, "\\d+")) %>% 
    mutate(g_idx = as.numeric(g_idx))
}

join_ext_param_num <- function(stanmod,par){
  extract_param(stanmod,par) %>%  
    rownames_to_column('x') %>% 
    mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>% 
    mutate(g_idx=str_split_fixed(x,'\\[', 3)[,2]) %>% 
    mutate(g_idx = str_extract(g_idx, "\\d+")) %>%
    mutate(g_idx = as.numeric(g_idx)) %>%
    dplyr::dplyr::select(-x)
}

# Load data -----------------------------------------------------------------------------------
stan_post <- readRDS(here('Output','stanMod_output_3.rds'))
stan_data <- readRDS(here('Output','stan_data_input_3.rds'))


# ESS -----------------------------------------------------------------------------------------
stan_post_df <- extract_param(stan_post)

# Determine the range of ESS
ess_range <- stan_post_df$n_eff %>% 
  range(na.rm = TRUE) %>% log() %>% 
  floor() %>% as.data.frame() %>% 
  setNames('ESS_range') %>% 
  mutate(ESS_range = if_else(row_number(.) == 2, ESS_range + 1, ESS_range)) %>% 
  mutate(ESS_range=exp(ESS_range)) %>% 
  pull(ESS_range)


ess_breaks <-
seq(from = sqrt(ess_range[1]), to = sqrt(ess_range[2]), length.out = 4) %>% 
  as.data.frame() %>% 
  setNames('ESS_breaks') %>% 
  mutate(ESS_breaks=ESS_breaks^2) %>% 
  mutate(ESS_breaks=round_by10(ESS_breaks)) %>% 
  pull(ESS_breaks)



# Diagnostics & Goodness-of-Fit
p1 <-
stan_post_df %>% 
  ggplot()+
  geom_histogram(aes(x=n_eff),binwidth = 2,fill = "skyblue", color = "black")+
  scale_x_sqrt(labels=scientific_10,limits=ess_range,breaks=ess_breaks) +
  theme_bw()+
  scale_y_sqrt()+
  xlab('ESS (Effective Sampling Size)')+
  ylab('Number of\nparameters')+
  theme(
    axis.title.x = element_text(size=13),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=14)
  )


# Rhat ----------------------------------------------------------------------------------------

# Determine the range of Rhat
rhat_range <-
stan_post_df$Rhat %>% 
  as.data.frame() %>% 
  setNames('Rhat') %>% 
  mutate(nparam = nrow(.),
         sd=sd(Rhat, na.rm = TRUE),
         max = max(Rhat, na.rm = TRUE),
         min = min(Rhat, na.rm = TRUE)) %>% 
  summarise(nparam=first(nparam),
            sd = first(sd),
            max = first(max),
            min = first(min)) %>% 
  mutate(binwidth=(sd*((nparam+100)^(1/3)))/100)
  

rhat_breaks <-
seq(from = rhat_range$min, to = rhat_range$max, length.out = 4) %>%
  tibble(values = .) %>%
  mutate(values = ifelse(abs(values - 1) == min(abs(values - 1)), 1, values)) %>% 
  mutate(values = round(values,round(log10(rhat_range$binwidth)*-1))) %>% 
  pull(values)

p2 <-
stan_post_df %>% 
  ggplot()+
  geom_histogram(aes(x=Rhat),binwidth = rhat_range$binwidth,fill = "skyblue", color = "black")+
  theme_bw()+
  scale_x_continuous(limits=c(rhat_range$min,rhat_range$max),labels = rhat_breaks,breaks = rhat_breaks)+
  scale_y_sqrt()+
  xlab('Rhat')+
  theme(
    axis.title.x = element_text(size=13),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    axis.title.y = element_blank()
  )

pp1 <- cowplot::plot_grid(p1,p2,nrow = 1)

pp1
# ggsave(here('plots','Diagnostic_fig_1_convergance.jpg'),sfxxx1,height = 6,width = 12)

# Tree depth ----------------------------------------------------------------------------------
sampler_params <- get_sampler_params(stan_post, inc_warmup = T)

warmup_it <- stan_post@stan_args[[1]]$warmup
samp_it <- stan_post@stan_args[[1]]$iter

# Extract tree depth for each chain
tre <- lapply(sampler_params, function(x) x[, "treedepth__"])
acc <- lapply(sampler_params, function(x) x[, "accept_stat__"])
sts <- lapply(sampler_params, function(x) x[, "stepsize__"])
div <- lapply(sampler_params, function(x) x[, "divergent__"])
log_lik <- get_logposterior(stan_post, inc_warmup = TRUE)

tree_depths_df <- 
  tibble(Iteration = rep(1:length(tre[[1]]),length(tre)),
    Chain = rep(1:length(tre),each = length(tre[[1]])),
         Tree_Depth = unlist(tre))
acc_df <- 
  tibble(Iteration = rep(1:length(acc[[1]]),length(acc)),
         Chain = rep(1:length(acc),each = length(acc[[1]])),
         Tree_Depth = unlist(acc))
sts_df <- 
  tibble(Iteration = rep(1:length(sts[[1]]),length(sts)),
         Chain = rep(1:length(sts),each = length(sts[[1]])),
         Tree_Depth = unlist(sts))
div_df <- 
  tibble(Iteration = rep(1:length(div[[1]]),length(div)),
         Chain = rep(1:length(div),each = length(div[[1]])),
         Tree_Depth = unlist(div))

lik_df <- 
  tibble(Iteration = rep(1:length(log_lik[[1]]),length(log_lik)),
         Chain = rep(1:length(log_lik),each = length(log_lik[[1]])),
         Tree_Depth = unlist(log_lik))

chain_col_map <- setNames(
  c('#F0944C',
    '#A3A55F',
    '#b4e5dd',
    '#C47688',
    '#B2B9BC',
    '#2E92A2'),c(1:6))

  
p1 <-
ggplot() +
  geom_line(data=tree_depths_df %>% filter(Iteration>=warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5) +
  geom_line(data=tree_depths_df %>% filter(Iteration<warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5,alpha=0.6) +
  labs(x = "Iteration",y = "Tree depth",color = "Chain") +
  theme_bw()+
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = warmup_it, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.2)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"),
    legend.position = 'none')

p2 <-
ggplot() +
  geom_point(data=div_df %>% filter(Iteration>=warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain))) +
  geom_point(data=div_df %>% filter(Iteration<warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),alpha=0.6) +
  labs(x = "Iteration",y = "Divergence\ntransition",color = "Chain") +
  scale_y_continuous(limits = c(-0.25,1.25), breaks=c(0,1),labels = c('No','Yes'))+
  theme_bw()+
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = warmup_it, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.2)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"),
    legend.position = 'none')


p3 <-
ggplot() +
  geom_line(data=sts_df %>% filter(Iteration>=warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5) +
  geom_line(data=sts_df %>% filter(Iteration<warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5,alpha=0.6) +
  labs(x = "Iteration",y = "Step size",color = "Chain") +
  theme_bw()+
  ylim(quantile(sts_df$Tree_Depth,0.0005),quantile(sts_df$Tree_Depth,0.9995))+
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = warmup_it, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.2)+
  theme(
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"),
    legend.position = 'none')

p4 <-
ggplot() +
  geom_line(data=acc_df %>% filter(Iteration>=warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5) +
  geom_line(data=acc_df %>% filter(Iteration<warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5,alpha=0.6) +
  labs(x = "Iteration",y = "Accepctance rate",color = "Chain") +
  theme_bw()+
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = warmup_it, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.2)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=13),
    axis.title.y = element_text(size=14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"),
    legend.position = 'none')


lik_range <-
lik_df %>% filter(Iteration>=warmup_it) %>% 
  mutate(lik_mean=mean(Tree_Depth),
         lik_min=min(Tree_Depth),
         lik_max=max(Tree_Depth)) %>% 
  summarise(lik_mean=first(lik_mean),
            lik_min=first(lik_min),
            lik_max=first(lik_max)) %>% 
  mutate(diff=lik_max-lik_min)

precision <- (log10(abs(lik_range$lik_mean)) - log10(lik_range$diff)) %>% round()
if(precision<0){precision = 0}

scientific_10 <- function(x) {
  c <- log10(abs(x))
  expo <- floor(c)
  base <- round(x / 10^expo,precision)
  formatted <- paste0(base,'%*%',"10^", expo)
  str2expression(formatted)
}

p5 <-
ggplot() +
  geom_line(data=lik_df %>% filter(Iteration>=warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5) +
  geom_line(data=lik_df %>% filter(Iteration<warmup_it),aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 0.5,alpha=0.6) +
  labs(x = "Iteration",y = "Posterior likelihood",color = "Chain") +
  theme_bw()+
  scale_y_continuous(labels=scientific_10,limits=c(lik_range$lik_min,lik_range$lik_max))+
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = warmup_it, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.2)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=11),
    axis.title.y = element_text(size=14),
    legend.title = element_text(size = 16),
    legend.position = 'none')

legend <-
ggplot() +
  geom_line(data=lik_df,aes(x = Iteration, y = Tree_Depth, color = as.factor(Chain)),linewidth = 1) +
  theme_bw()+
  labs(x = "Iteration",y = "Accepctance rate",color = "Chain") +
  scale_color_manual(values = chain_col_map)+
  geom_rect(aes(xmin = 1, xmax = 200, ymin = -Inf, ymax = Inf, fill = "Warmup"),color = "black", alpha = 0.3) +
  geom_rect(aes(xmin = 200, xmax = 1000, ymin = -Inf, ymax = Inf, fill = "Sampling"), color = "black",alpha = 0.3) +
  scale_fill_manual(name = "Phase", values = c("Warmup" = "grey70", "Sampling" = "white")) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"))

leg <- cowplot::get_legend(legend)

pp2 <- cowplot::plot_grid(p5,p4,p1,p2,p3,ncol = 1,align = 'v',rel_heights = c(5,5,5,2.2,6))

pp3 <- cowplot::plot_grid(pp1,pp2,ncol = 1,rel_heights = c(1,3),labels = c('A','B'),label_x = 1)
df1 <- cowplot::plot_grid(pp3,leg,ncol = 2,rel_widths = c(7,1.5))

df1

ggsave(here('Plots','Diagnostic_Fig_1.jpg'),df1,height = 13,width = 9)


# # Prior Sensitivity ---------------------------------------------------------------------------


X_STATE <-
  extract_param(stan_post,'X_STATE') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::dplyr::select(-x) %>% mutate(param='X_STATE')

eta <-
  extract_param(stan_post,'eta') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::dplyr::select(-x) %>% mutate(param='eta')

omega <-
  extract_param(stan_post,'omega') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::dplyr::select(-x) %>% mutate(param='omega')

epsilon <-
  extract_param(stan_post,'epsilon') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::dplyr::select(-x) %>% mutate(param='epsilon')

tau <-
  extract_param(stan_post,'tau') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::dplyr::select(-x) %>% mutate(param='tau')

delta <-
  extract_param(stan_post,'delta') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='delta') %>% 
  slice(c(1:12,25:36))

rho <-
  extract_param(stan_post,'tau_raw') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='rho') %>% 
  slice(c(1,4))

phi <-
  extract_param(stan_post,'logit_phi') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='phi')

gamma_0 <-
  extract_param(stan_post,'gamma_0') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='gamma_0')

gamma_1 <-
  extract_param(stan_post,'gamma_1') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='gamma_1')

beta_0 <-
  extract_param(stan_post,'beta_0') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='beta_0')

beta_1 <-
  extract_param(stan_post,'beta_1') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  dplyr::select(-x) %>% mutate(param='beta_1')



X_STATE_prior <- data.frame(mean=rnorm(2e4,100,100)) %>%
	filter(mean>0) %>% 
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='X_STATE')

eta_prior <- data.frame(mean=rnorm(1e4,-2,4)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='eta')

omega_prior <- data.frame(mean=rnorm(1e4,6,1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='omega')

phi_prior <- data.frame(mean=rnorm(1e4,-1,1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='phi')

beta_0_prior <- data.frame(mean=rnorm(1e4,40,1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='beta_0')

beta_1_prior <- data.frame(mean=rnorm(1e4,-3,1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='beta_1')

gamma_0_prior <- data.frame(mean=rnorm(1e4,1,0.1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='gamma_0')

gamma_1_prior <- data.frame(mean=rnorm(1e4,0,0.1)) %>%
  filter(mean<0) %>% 
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='gamma_1')

tau_prior <- data.frame(mean=rnorm(2e4,0,1)) %>%
  filter(mean>0) %>% 
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='tau')

epsilon_prior <- data.frame(mean=rnorm(1e4,0,tau_prior$mean)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='epsilon')

rho_prior <- data.frame(mean=rnorm(1e4,0,1)) %>%
  filter(mean>0) %>% 
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='rho')

delta_prior <- data.frame(mean=rnorm(1e4,0,rho_prior$mean)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='delta')


post_jitter <-
bind_rows(X_STATE,eta,omega,epsilon,beta_0,beta_1,gamma_0,gamma_1,phi,tau,rho,delta) %>% 
  group_by(param) %>%
  mutate(nudge_y = (row_number() - 1) * 0.05,
         count = n()) %>% 
  ungroup() %>% 
  group_by(param) %>%
  mutate(nudge_yp = nudge_y/max(nudge_y)) %>%
  ungroup() %>% 
  mutate(nudge_y = if_else(count>10,nudge_yp,nudge_y)) %>%
  # mutate(nudge_y = if_else(is.na(nudge_y),0,nudge_y)) %>%
  # mutate(nudge_y = nudge_y-0.4) %>%
  as.data.frame()
  
  # mutate(nudge_y=nudge_y-0.4)
priors <- bind_rows(X_STATE_prior,eta_prior,omega_prior,epsilon_prior,beta_0_prior,beta_1_prior,gamma_0_prior,gamma_1_prior,phi_prior,tau_prior,rho_prior,delta_prior)
# jitter <- position_jitter(width = 0, height = 0.4, seed = 123)


pp1 <-
post_jitter %>% filter(param%in%c('beta_0','beta_1')) %>%
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('beta_0','beta_1')) %>% pull(nudge_y))) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('beta_0','beta_1')) %>% pull(nudge_y)), width = 0.02, alpha = 1)+
    geom_point(data=priors %>% filter(param%in%c('beta_0','beta_1')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.2))+
    geom_errorbar(data=priors%>% filter(param%in%c('beta_0','beta_1')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.2)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "beta_0" = expression(beta[0]),
    "beta_1" = expression(beta[1])))

pp2 <-
post_jitter %>% filter(param%in%c('gamma_0','gamma_1')) %>% 
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('gamma_0','gamma_1')) %>% pull(nudge_y)),alpha=1) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('gamma_0','gamma_1')) %>% pull(nudge_y)),width = 0.02, alpha = 1)+
    geom_point(data=priors %>% filter(param%in%c('gamma_0','gamma_1')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.1))+
    geom_errorbar(data=priors%>% filter(param%in%c('gamma_0','gamma_1')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.1)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "gamma_0" = expression(gamma[0]),
    "gamma_1" = expression(gamma[1])))

pp3 <-
post_jitter %>% filter(param%in%c('epsilon','tau')) %>% 
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('epsilon','tau')) %>% pull(nudge_y)),
              alpha=1) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('epsilon','tau')) %>% pull(nudge_y)),
                width = 0.02, alpha = 0.6)+
    geom_point(data=priors %>% filter(param%in%c('epsilon','tau')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.2))+
    geom_errorbar(data=priors%>% filter(param%in%c('epsilon','tau')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.2)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "tau" = expression(tau),
    "epsilon" = expression(epsilon)))

pp4 <-
post_jitter %>% filter(param%in%c('delta','rho')) %>% 
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('delta','rho')) %>% pull(nudge_y)),
              alpha=1) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('delta','rho')) %>% pull(nudge_y)),
                width = 0.02, alpha = 0.6)+
    geom_point(data=priors %>% filter(param%in%c('delta','rho')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.2))+
    geom_errorbar(data=priors%>% filter(param%in%c('delta','rho')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.2)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "delta" = expression(delta),
    "rho" = expression(rho),
    "tau" = expression(tau)
    ))

pp5 <-
post_jitter %>% filter(param%in%c('phi','omega')) %>% 
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('phi','omega')) %>% pull(nudge_y)),
              alpha=1) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('phi','omega')) %>% pull(nudge_y)),
                width = 0.02, alpha = 0.6)+
    geom_point(data=priors %>% filter(param%in%c('phi','omega')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.1))+
    geom_errorbar(data=priors%>% filter(param%in%c('phi','omega')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.1)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "omega" = expression(omega),
    "phi" = expression(phi)))

pp6 <-
post_jitter %>% filter(param%in%c('eta','X_STATE')) %>% 
  ggplot() +
  geom_jitter(aes(x = mean, y = param),
              position = position_nudge(y = post_jitter %>% filter(param%in%c('eta','X_STATE')) %>% pull(nudge_y)),
              alpha=1) +
  geom_errorbar(aes(x = mean, y = param, xmin = `2.5%`, xmax = `97.5%`), 
                position = position_nudge(y = post_jitter %>% filter(param%in%c('eta','X_STATE')) %>% pull(nudge_y)),
                width = 0.02, alpha = 0.6)+
    geom_point(data=priors %>% filter(param%in%c('eta','X_STATE')),aes(x=mean,y=param),color='tomato2',size=3,
               position = position_nudge(y = -0.1))+
    geom_errorbar(data=priors%>% filter(param%in%c('eta','X_STATE')),aes(y = param,xmin = lower, xmax = upper),color='red', width = 0.1,
                  position = position_nudge(y = -0.1)) +
    facet_wrap(~param,scales = "free",ncol = 2)+
    theme_bw()+
    theme(strip.text = element_blank(),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=17),
          axis.title.x =element_blank(),
          axis.title.y = element_blank())+
  scale_y_discrete(labels = c(
    "X_STATE" = expression('X'),
    "eta" = expression(eta)))

legend <-
  data_frame(mean=c(1,1),param=c(1,1),lik=c('Prior','Posterior')) %>%
  ggplot() +
  geom_point(aes(x = mean, y = param,colour = lik),size=3) +
  geom_line(aes(x = mean, y = param,colour = lik),size=1.5)+
  theme_bw()+
  labs(colour = "Likelihood")+
  scale_color_manual(values=c('Prior'='tomato2','Posterior'='black'))+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.7, "lines"))
leg <- cowplot::get_legend(legend)

p1 <- cowplot::plot_grid(pp1,pp2,pp3,pp4,pp5,pp6,ncol = 1,align='v')

df2 <- cowplot::plot_grid(p1,leg,ncol = 2,rel_widths = c(6,1))

df2
ggsave(here('Plots','Diagnostic_Fig_2.jpg'),df2,height = 12,width = 9)


stan_data$N




# Posterior Predictive Checks -----------------------------------------------------------------

st_pred_pd <- data.frame(S_q=seq(-2,13,by=0.1)) %>% 
  mutate(phi=join_ext_param(stan_post,'logit_phi') %>% pull(mean) %>% inverselogit()) %>%
  mutate(phi_lo=join_ext_param(stan_post,'logit_phi') %>% pull(`2.5%`) %>% inverselogit()) %>%
  mutate(phi_up=join_ext_param(stan_post,'logit_phi') %>% pull(`97.5%`) %>% inverselogit()) %>%
  mutate(psi_pred=1-exp(exp(S_q)*phi*-1)) %>% 
  mutate(psi_pred_lo=1-exp(exp(S_q)*phi_lo*-1)) %>% 
  mutate(psi_pred_up=1-exp(exp(S_q)*phi_up*-1))
  
st_data_pd <- cbind(stan_data$S_q,stan_data$Z_qst) %>% 
  as.data.frame() %>% 
  setNames(c('S_q','Z')) %>% 
  mutate(phi=join_ext_param(stan_post,'logit_phi') %>% pull(mean) %>% inverselogit()) %>%
  mutate(phi_lo=join_ext_param(stan_post,'logit_phi') %>% pull(`2.5%`) %>% inverselogit()) %>%
  mutate(phi_up=join_ext_param(stan_post,'logit_phi') %>% pull(`97.5%`) %>% inverselogit()) %>%
  group_by(S_q) %>% 
  mutate(Z_mean=mean(Z)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  mutate(psi_pred=1-exp(exp(S_q)*phi*-1)) %>% 
  mutate(psi_pred_lo=1-exp(exp(S_q)*phi_lo*-1)) %>% 
  mutate(psi_pred_up=1-exp(exp(S_q)*phi_up*-1))

W_data_pd <-
cbind(stan_data$Z_qen_wat,stan_data$j_qen_wat_idx) %>% 
  as.data.frame() %>% 
  setNames(c('Z','g_idx')) %>% 
  left_join(.,
            join_ext_param(stan_post,'log_W') %>% 
              rename(log_W='mean') %>% 
              dplyr::select(log_W,g_idx),by='g_idx') %>% 
  mutate(phi=join_ext_param(stan_post,'logit_phi') %>% pull(mean) %>% inverselogit()) %>%
  mutate(phi_lo=join_ext_param(stan_post,'logit_phi') %>% pull(`2.5%`) %>% inverselogit()) %>%
  mutate(phi_up=join_ext_param(stan_post,'logit_phi') %>% pull(`97.5%`) %>% inverselogit()) %>%
  group_by(g_idx) %>% 
  mutate(Z_mean=mean(Z)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  mutate(psi_pred=1-exp(exp(log_W)*phi*-1)) %>% 
  mutate(psi_pred_lo=1-exp(exp(log_W)*phi_lo*-1)) %>% 
  mutate(psi_pred_up=1-exp(exp(log_W)*phi_up*-1))


A_data_pd <-
cbind(stan_data$Z_qen_air,stan_data$j_qen_air_idx) %>% 
  as.data.frame() %>% 
  setNames(c('Z','g_idx')) %>% 
  left_join(.,
            join_ext_param(stan_post,'log_A') %>% 
              rename(log_A='mean') %>% 
              dplyr::select(log_A,g_idx),by='g_idx') %>% 
  mutate(phi=join_ext_param(stan_post,'logit_phi') %>% pull(mean) %>% inverselogit()) %>%
  mutate(phi_lo=join_ext_param(stan_post,'logit_phi') %>% pull(`2.5%`) %>% inverselogit()) %>%
  mutate(phi_up=join_ext_param(stan_post,'logit_phi') %>% pull(`97.5%`) %>% inverselogit()) %>%
  group_by(g_idx) %>% 
  mutate(Z_mean=mean(Z)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  mutate(psi_pred=1-exp(exp(log_A)*phi*-1)) %>% 
  mutate(psi_pred_lo=1-exp(exp(log_A)*phi_lo*-1)) %>% 
  mutate(psi_pred_up=1-exp(exp(log_A)*phi_up*-1))


pred_cm <-
data.frame(S_q=rep(seq(-2,14,by=0.1),5),
           p=rep(c(1:5),each=length(seq(-2,14,by=0.1)))) %>% 
  mutate(gamma_0=join_ext_param(stan_post,'gamma_0') %>% pull(mean)) %>% 
  mutate(gamma_1=join_ext_param(stan_post,'gamma_1') %>% pull(mean)) %>% 
  mutate(beta_1=join_ext_param(stan_post,'beta_1') %>% pull(mean)) %>% 
  left_join(.,
            join_ext_param_num(stan_post,'beta_0') %>% 
              dplyr::select(mean,g_idx) %>% 
              rename(beta_0='mean'),
            by=c('p'='g_idx')) %>% 
  mutate(pred_mu=beta_0+beta_1*S_q) %>% 
  mutate(pred_sigma=exp(gamma_0+(gamma_1 * S_q))) %>% 
  rowwise() %>%
  mutate(pred_ct = rnorm(1, pred_mu, pred_sigma)) %>%
  ungroup() %>% as.data.frame()


st_data_cm <-
cbind(stan_data$S_q_p,stan_data$R_qst,stan_data$plate_st_idx) %>% 
  as.data.frame() %>% 
  setNames(c('S_q','Y','p')) %>% 
  mutate(gamma_0=join_ext_param(stan_post,'gamma_0') %>% pull(mean)) %>% 
  mutate(gamma_1=join_ext_param(stan_post,'gamma_1') %>% pull(mean)) %>% 
  mutate(beta_1=join_ext_param(stan_post,'beta_1') %>% pull(mean)) %>% 
  left_join(.,
            join_ext_param_num(stan_post,'beta_0') %>% 
              dplyr::select(mean,g_idx) %>% 
              rename(beta_0='mean'),
            by=c('p'='g_idx')) %>% 
  mutate(mu=beta_0+(beta_1*S_q)) %>% 
  mutate(sigma=exp(gamma_0+(gamma_1*S_q))) %>% 
  mutate(pred_ct=mu) %>% 
  mutate(pred_ct_lo=mu+(sigma^2)) %>% 
  mutate(pred_ct_up=mu-(sigma^2))

W_data_cm <-
cbind(stan_data$R_qen_wat,stan_data$j_qen_wat_p_idx,stan_data$plate_en_wat_idx)  %>% 
  as.data.frame() %>% 
  setNames(c('Y','idx','p')) %>% 
  left_join(.,
            join_ext_param(stan_post,'log_W') %>% 
              rename(log_W='mean') %>% 
              dplyr::select(log_W,g_idx),by=c('idx'='g_idx')) %>% 
  mutate(gamma_0=join_ext_param(stan_post,'gamma_0') %>% pull(mean)) %>% 
  mutate(gamma_1=join_ext_param(stan_post,'gamma_1') %>% pull(mean)) %>% 
  mutate(beta_1=join_ext_param(stan_post,'beta_1') %>% pull(mean)) %>% 
  left_join(.,
            join_ext_param_num(stan_post,'beta_0') %>% 
              dplyr::select(mean,g_idx) %>% 
              rename(beta_0='mean'),
            by=c('p'='g_idx')) %>% 
  mutate(mu=beta_0+(beta_1*log_W)) %>% 
  mutate(sigma=exp(gamma_0+(gamma_1*log_W))) %>% 
  mutate(pred_ct=mu) %>% 
  mutate(pred_ct_lo=mu+(sigma^2)) %>% 
  mutate(pred_ct_up=mu-(sigma^2))


A_data_cm <-
cbind(stan_data$R_qen_air,stan_data$j_qen_air_p_idx,stan_data$plate_en_air_idx)  %>% 
  as.data.frame() %>% 
  setNames(c('Y','idx','p')) %>% 
  left_join(.,
            join_ext_param(stan_post,'log_A') %>% 
              rename(log_A='mean') %>% 
              dplyr::select(log_A,g_idx),by=c('idx'='g_idx')) %>% 
  mutate(gamma_0=join_ext_param(stan_post,'gamma_0') %>% pull(mean)) %>% 
  mutate(gamma_1=join_ext_param(stan_post,'gamma_1') %>% pull(mean)) %>% 
  mutate(beta_1=join_ext_param(stan_post,'beta_1') %>% pull(mean)) %>% 
  left_join(.,
            join_ext_param_num(stan_post,'beta_0') %>% 
              dplyr::select(mean,g_idx) %>% 
              rename(beta_0='mean'),
            by=c('p'='g_idx')) %>% 
  mutate(mu=beta_0+(beta_1*log_A)) %>% 
  mutate(sigma=exp(gamma_0+(gamma_1*log_A))) %>% 
  mutate(pred_ct=mu) %>% 
  mutate(pred_ct_lo=mu+(sigma^2)) %>% 
  mutate(pred_ct_up=mu-(sigma^2))


N_data <- bind_cols(stan_data$N %>% as.data.frame() %>% setNames('N'),
  extract_param(stan_post,'lambda') %>% dplyr::select(mean,`2.5%`,`97.5%`) %>% setNames(c('lambda','lambda_lo','lambda_up'))) %>% 
  group_by(row_number(.)) %>% 
  mutate(N_pred=median(rnegbin(1e4,lambda,20)),
         N_pred_lo=median(rnegbin(1e4,lambda_lo,20)),
         N_pred_up=median(rnegbin(1e4,lambda_up,20))) 
  # mutate(N_pred=median(rpois(1e4,lambda)),
  #        N_pred_lo=median(rpois(1e4,lambda_lo)),
  #        N_pred_up=median(rpois(1e4,lambda_up))) 

legend <- ggplot(data.frame(x = rep(1, 4),y = rep(1, 4),group = factor(rep(1:4))), 
       aes(x = x, y = y, group = group, color = group,linetype = group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("black", "deepskyblue2", "tomato", "darkorchid"), 
                     labels = c("Standard samples", "Water samples", "Air samples", "Visual samples")) +
  scale_linetype_manual(values = c(1,1,1,1), 
                        labels = c("Standard samples", "Water samples", "Air samples", "Visual samples")) +
  labs(title = '', 
       color = "Observation method", 
       linetype = "Observation method") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))

leg <- cowplot::get_legend(legend)

p1 <-
  ggplot()+
  geom_point(data=st_data_pd %>% distinct(S_q,Z_mean,psi_pred),
             aes(x=Z_mean,y=psi_pred), alpha=0.6)+
  geom_errorbar(data=st_data_pd %>% distinct(S_q,Z_mean,psi_pred,psi_pred_lo,psi_pred_up),
                aes(y=psi_pred,ymin = psi_pred_lo,ymax = psi_pred_up,x=Z_mean),color='black',alpha=0.4,width=0.02)+
  geom_point(data=A_data_pd %>% distinct(g_idx,Z_mean,psi_pred),aes(x=Z_mean,y=psi_pred), alpha=0.6,color='tomato')+
  geom_errorbar(data=A_data_pd %>% distinct(g_idx,Z_mean,psi_pred,psi_pred_lo,psi_pred_up),
                aes(y=psi_pred,ymin = psi_pred_lo,ymax = psi_pred_up,x=Z_mean),color='tomato',alpha=0.24,width=0.02)+
  geom_point(data=W_data_pd %>% distinct(g_idx,Z_mean,psi_pred),aes(x=Z_mean,y=psi_pred),alpha=0.6,color='deepskyblue2')+
  geom_errorbar(data=W_data_pd %>% distinct(g_idx,Z_mean,psi_pred,psi_pred_lo,psi_pred_up),
                aes(y=psi_pred,ymin = psi_pred_lo,ymax = psi_pred_up,x=Z_mean),color='deepskyblue2',alpha=0.4,width=0.02)+
  geom_abline(intercept = 0,slope=1,lty=2)+
  labs(x='Observed probability of detection',y='Predicted probability\nof detection')+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13))

p2 <- A_data_cm %>% 
  ggplot()+
  geom_point(data=st_data_cm,aes(y=pred_ct,x=Y),color='black',alpha=0.4)+
  geom_errorbar(data=st_data_cm,aes(y=pred_ct,ymin = pred_ct_lo,ymax = pred_ct_up,x=Y),color='black',alpha=0.4)+
  geom_point(data=W_data_cm,aes(y=pred_ct,x=Y),color='deepskyblue2',alpha=0.4)+
  geom_errorbar(data=W_data_cm,aes(y=pred_ct,ymin = pred_ct_lo,ymax = pred_ct_up,x=Y),color='deepskyblue2',alpha=0.4)+
  geom_point(aes(y=pred_ct,x=Y),color='tomato',alpha=0.4)+
  geom_errorbar(aes(y=pred_ct,ymin = pred_ct_lo,ymax = pred_ct_up,x=Y),color='tomato',alpha=0.4)+
  geom_abline(intercept = 0,slope=1,lty=2)+
  labs(x='Observed Ct values',y='Predicted Ct values')+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13))

p3 <- N_data %>% 
  ggplot()+
  geom_point(aes(x=N,y=N_pred),color='darkorchid')+
  geom_errorbar(aes(y=N_pred,ymin = N_pred_lo,ymax = N_pred_up,x=N),color='darkorchid')+
  geom_abline(intercept = 0,slope=1,lty=2)+
  labs(x='Observed N values',y='Predicted N values')+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13))


pp1 <- cowplot::plot_grid(p1,p2,p3,ncol = 1,align = 'v',labels=c('A','B','C'))
df3 <- cowplot::plot_grid(pp1,leg,ncol = 2,rel_widths = c(5,2))

df3
ggsave(here('Plots','Diagnostic_Fig_3.jpg'),df3,width = 7,height = 12)


# qPCR assay efficiency -----------------------------------------------------------------------

legend <- ggplot(data.frame(x = rep(1, 3),y = rep(1, 3),group = factor(rep(1:3))), 
                 aes(x = x, y = y, group = group, color = group,linetype = group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("black", "deepskyblue2", "tomato"), 
                     labels = c("Standard", "Water", "Air")) +
  scale_linetype_manual(values = c(1,1,1,1), 
                        labels = c("Standard", "Water", "Air")) +
  labs(title = '', 
       color = "Samples", 
       linetype = "Samples") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))

leg <- cowplot::get_legend(legend)

p4 <-
ggplot()+
  geom_line(data=st_pred_pd,aes(x=exp(S_q),y=psi_pred),size=1,alpha=0.5)+
  geom_line(data=st_pred_pd,aes(x=exp(S_q),y=psi_pred_lo),size=1,lty=2,alpha=0.3)+
  geom_line(data=st_pred_pd,aes(x=exp(S_q),y=psi_pred_up),size=1,lty=2,alpha=0.3)+
  geom_jitter(data=st_data_pd,aes(x=exp(S_q),y=Z),width = 0.1,height = 0.05,alpha=0.5)+
  geom_point(data=st_data_pd,aes(x=exp(S_q),y=Z_mean),size=3,color='black',alpha=0.1)+
  geom_jitter(data=W_data_pd,aes(x=exp(log_W),y=Z),width = 0.1,height = 0.05,color='deepskyblue2',alpha=0.7)+
  geom_point(data=W_data_pd,aes(x=exp(log_W),y=Z_mean),size=3,color='deepskyblue2')+
  geom_jitter(data=A_data_pd,aes(x=exp(log_A),y=Z),width = 0.1,height = 0.05,color='tomato',alpha=0.7)+
  geom_point(data=A_data_pd,aes(x=exp(log_A),y=Z_mean),size=3,color='tomato')+
  scale_x_log10(labels=scientific)+
  labs(x='DNA concentration',y='Probability of detection')+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13))

p5 <-
pred_cm %>% mutate(p=as.factor(p)) %>%
  ggplot()+
  geom_line(aes(x=exp(S_q),y=pred_mu,colour = p),alpha=0.4,lty=2)+
  scale_color_manual(values=rep('black',5))+
  geom_point(data=st_data_cm,aes(y=pred_ct,x=exp(S_q),alpha=0.4))+
  geom_point(data=W_data_cm,aes(y=pred_ct,x=exp(log_W)),color='deepskyblue2',alpha=0.9)+
  geom_errorbar(data=W_data_cm,aes(y=pred_ct,ymin = pred_ct_lo,ymax = pred_ct_up,x=exp(log_W)),color='deepskyblue2',alpha=0.4)+
  geom_point(data=A_data_cm,aes(y=pred_ct,x=exp(log_A)),color='tomato',alpha=0.9)+
  geom_errorbar(data=A_data_cm,aes(y=pred_ct,ymin = pred_ct_lo,ymax = pred_ct_up,x=exp(log_A)),color='tomato',alpha=0.4)+
  theme_bw()+
  labs(x='DNA concentration',y='Ct')+
  scale_x_log10(labels=scientific)+
  theme(legend.position = 'none',
    axis.title = element_text(size=14),
        axis.text = element_text(size=13))

pp1 <- cowplot::plot_grid(p4,p5,ncol = 1,align = 'v',labels=c('A','B'))
df3 <- cowplot::plot_grid(pp1,leg,ncol = 2,rel_widths = c(5,1))

ggsave(here('Plots','Supplementary_Fig_1.jpg'),df3,width = 10,height = 10)
