devtools::install_github("gledguri/QM")
library(QM)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(rstan);options(mc.cores = parallel::detectCores())
library(here)

M2 <- load_model('M2')

qpcr_co <- readRDS(here('qPCR_code_coho_chinook','qpcr_coho.rds')) %>% 
	filter(FilterType=='Water filter MCE'|Task=='Standard') %>% 
	mutate(Task=if_else(Task=='Standard','STANDARD','UNKNOWN')) %>% 
	filter(!grepl('field negative',Location))

qpcr_ch <- readRDS(here('qPCR_code_coho_chinook','qpcr_chin.rds')) %>% 
	filter(FilterType=='Water filter MCE'|Task=='Standard') %>% 
	mutate(Task=if_else(Task=='Standard','STANDARD','UNKNOWN')) %>% 
	filter(!grepl('field negative',Location)) %>% 
	mutate(plate_idx='1')

stan_data_M2_co <- prep_stan_M2(
	qpcr_data = qpcr_co,
	sample_type = "Task",
	Ct = "Cq",
	sample_name_column = "Date",
	standard_concentration = "Quantity",
	plate_index = 'plate_idx')

stan_data_M2_ch <- prep_stan_M2(
	qpcr_data = qpcr_ch,
	sample_type = "Task",
	Ct = "Cq",
	sample_name_column = "Date",
	standard_concentration = "Quantity",
	plate_index = 'plate_idx')

M2_output_co <- Run_Model(stan_object = M2, stan_data = stan_data_M2_co)
M2_output_ch <- Run_Model(stan_object = M2, stan_data = stan_data_M2_ch)

co_plot_dat <- extract_est_conc(M2_output_co) %>% as_tibble() %>% 
	mutate(Year=substr(Sample_name,0,4)) %>% 
	mutate(Month=substr(Sample_name,5,6)) %>% 
	mutate(Day=substr(Sample_name,7,8)) %>% 
	mutate(Date_raw=paste0(Year,'-',Month,'-',Day)) %>% 
	mutate(Date=as.Date(Date_raw)) 

ch_plot_dat <- extract_est_conc(M2_output_ch) %>% as_tibble() %>% 
	mutate(Year=substr(Sample_name,0,4)) %>% 
	mutate(Month=substr(Sample_name,5,6)) %>% 
	mutate(Day=substr(Sample_name,7,8)) %>% 
	mutate(Date_raw=paste0(Year,'-',Month,'-',Day)) %>% 
	mutate(Date=as.Date(Date_raw)) 

co_plot_dat %>% 
	ggplot()+
	geom_point(aes(x=Date,y=C_est_log),color='deepskyblue3')+
	geom_smooth(aes(x=Date,y=C_est_log),se=F,span = 0.2,color='deepskyblue3',lty=2)+
	geom_errorbar(aes(x=Date,ymin=`C_est_log_2.5%CI`,ymax=`C_est_log_97.5%CI`),color='deepskyblue3')+
	geom_point(data=ch_plot_dat, aes(x=Date,y=C_est_log),color='tomato2')+
	geom_smooth(data=ch_plot_dat,aes(x=Date,y=C_est_log),se=F,span = 0.2,color='tomato2',lty=2)+
	geom_errorbar(data=ch_plot_dat,aes(x=Date,ymin=`C_est_log_2.5%CI`,ymax=`C_est_log_97.5%CI`),color='tomato2')+
	theme_bw()

