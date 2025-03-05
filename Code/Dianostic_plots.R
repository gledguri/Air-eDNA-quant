library(ggridges)
library(rstan)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

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
  select(-x) %>% mutate(param='X_STATE')

eta <-
  extract_param(stan_post,'eta') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='eta')

omega <-
  extract_param(stan_post,'omega') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='omega')

epsilon <-
  extract_param(stan_post,'epsilon') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='epsilon')

tau <-
  extract_param(stan_post,'tau') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='tau')

delta <-
  extract_param(stan_post,'bio_rep_RE') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='delta') %>% 
  slice(c(1:12,25:36))

rho <-
  extract_param(stan_post,'tau_bio_rep') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='rho') %>% 
  slice(c(1,4))

phi <-
  extract_param(stan_post,'logit_phi') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='phi')

gamma_0 <-
  extract_param(stan_post,'gamma_0') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='gamma_0')

gamma_1 <-
  extract_param(stan_post,'gamma_1') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='gamma_1')

beta_0 <-
  extract_param(stan_post,'beta_0') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='beta_0')

beta_1 <-
  extract_param(stan_post,'beta_1') %>%
  rownames_to_column('x') %>%
  mutate(param=str_split_fixed(x,'\\[', 3)[,1]) %>%
  select(-x) %>% mutate(param='beta_1')



X_STATE_prior <- data.frame(mean=rgamma(1e4,10,1)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='X_STATE')

eta_prior <- data.frame(mean=rnorm(1e4,-2,4)) %>%
  mutate(lower=quantile(mean,probs = 0.025)) %>%
  mutate(upper=quantile(mean,probs = 0.975)) %>%
  summarise(mean=mean(mean),lower=mean(lower),upper=mean(upper)) %>%
  mutate(param='eta')

omega_prior <- data.frame(mean=rnorm(1e4,0,1)) %>%
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

ggsave(here('Plots','Diagnostic_Fig_2.jpg'),df2,height = 12,width = 9)


