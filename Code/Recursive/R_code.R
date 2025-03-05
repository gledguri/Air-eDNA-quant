library(here)
library(dplyr)
library(ggplot2)
library(QM);#devtools::install_github("gledguri/QM")
library(rstan);options(mc.cores = parallel::detectCores())
library(tibble)
library(stringr)
library(readr)
library(purrr)
library(ggu.base.fun)

make_plate_outlook <- function(df,samp,Well){
	temp <- 
		df %>% 
		rename(Label=all_of(samp)) %>% 
		rename(Well=all_of(Well)) %>% 
		mutate(Letter=substr(Well,0,1)) %>% 
		mutate(No=substr(Well,2,3)) %>% 
		select(Label,Letter,No) %>% 
		pivot_wider(
			names_from = No,
			values_from = Label) %>% 
		arrange(Letter) %>%
		select('Letter',as.character(as.numeric(colnames(.)) %>% sort())) 
	return(temp)
}

make_index <- function(data, index_variable, index_name) {
	data %>%
		group_by(across(all_of(index_variable))) %>%
		mutate(!!index_name := cur_group_id()) %>%
		ungroup()
}
combine_index <- function(data, index_variables, combined_index_name) {
	data %>% 
		mutate(across(all_of(index_variables), 
									~ paste0(
										tolower(LETTERS[((. - 1) %/% (26^2)) %% 26 + 1]), # First letter
										tolower(LETTERS[((. - 1) %/% 26) %% 26 + 1]),     # Second letter
										tolower(LETTERS[(. - 1) %% 26 + 1])              # Third letter
									), 
									.names = "letter_{col}")) %>% # Convert to letters
		rowwise() %>% # Row-wise operation for combining letters
		mutate(combined = paste(across(starts_with("letter_")), collapse = "_")) %>%
		ungroup() %>%
		group_by(combined) %>%
		mutate(!!combined_index_name := cur_group_id()) %>% # Create numeric index
		ungroup() %>%
		select(-starts_with("letter_"), -combined) # Remove intermediary letter-based columns
}
scientific <- function(x) {
	c <- log10(x)
	formatted <- paste0("10^", round(c, 2))
	str2expression(formatted)
}

inverselogit <- function (x) 
{
	return(1/(1 + exp(-x)))
}

file_paths <- list.files('Data', pattern = '\\.csv$', full.names = TRUE)
file_paths_meta <- file_paths[grepl('metadata',file_paths)]
file_paths_data <- file_paths[!(file_paths%in%file_paths_meta)]

meta <- lapply(seq_along(file_paths_meta), function(idx) {
	read_csv(file_paths_meta[idx]) %>% mutate(plate = idx)}) %>%
	bind_rows()

qpcr <- lapply(seq_along(file_paths_data), function(idx) {
	read_csv(file_paths_data[idx]) %>% mutate(plate = idx)}) %>%
	bind_rows()

qpcr <- qpcr %>% 
	mutate(Sample=ifelse(Sample=='cohoSHEL','coho_shel',Sample)) %>% 
	mutate(Sample=ifelse(Sample=='chinookSHEL','chinook_shel',Sample))

qpcr %>% filter(Task=='Standard') %>% 
	distinct(plate,Sample,Target) %>% 
	count(Sample,Target)


st_qpcr <- qpcr %>% filter(Task=='Standard')
en_qpcr <- qpcr %>% filter(Task=='Unknown')

# Correcting for annotations

qpcr_coho_1 <- qpcr %>% 
	filter(Sample=='coho_shel') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	mutate(Quantity=if_else(plate==7&Quantity==1000,Quantity*1e6,Quantity)) %>%
	mutate(Quantity=if_else(plate==7&Quantity==10000,Quantity/10,Quantity)) %>%
	mutate(Quantity=if_else(plate==7&Quantity==1000000000,Quantity/1e5,Quantity)) %>%
	mutate(Quantity=if_else(plate==7&Well=='E4',Quantity*10,Quantity)) %>%
	mutate(Quantity=if_else(plate==7&Well=='F5',Quantity*10,Quantity)) %>% 
	as.data.frame()


# Coho ----------------------------------------------------------------------------------------


qpcr_coho_2 <- qpcr %>% 
	filter(Sample=='cohoDUDA') %>% 
	rename(Target_2='Sample') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame() %>% 
	left_join(meta %>% filter(plate==1) %>% 
							mutate(plate=9),by=c('Well','plate'))
#Indexing
qpcr_coho_2 <- qpcr_coho_2 %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>% 
	combine_index(c('j_idx','t_idx','f_idx'),'idx')

qpcr_coho_2_st <- qpcr_coho_2 %>% filter(Task=='Standard') %>% 
	make_index("plate", "plate_idx") %>% 
	mutate(Ct = replace(Cq, Cq == "Undetermined", NA)) %>% 
	mutate(Ct = as.numeric(Ct)) %>% mutate(pres = 1, pres = replace(pres, is.na(Ct), 0)) %>% 
	mutate(st_conc = as.numeric(Quantity))


## Assay efficiency ----------------------------------------------------------------------------
stan_data_M1 <- list(N_st_q = nrow(qpcr_coho_2_st), 
										 N_st_qp = nrow(qpcr_coho_2_st %>% filter(pres==1)), 
										 N_plate = length(unique(qpcr_coho_2_st$plate_idx)),
										 plate_idx = as.numeric(qpcr_coho_2_st %>% filter(pres==1) %>% pull(plate_idx)),
										 Z_qst = as.integer(qpcr_coho_2_st$pres), 
										 S_q = log10(qpcr_coho_2_st$st_conc), 
										 R_qst = as.numeric(qpcr_coho_2_st %>% filter(pres==1) %>% pull(Ct)), 
										 S_q_p = qpcr_coho_2_st %>% filter(pres==1) %>%
										 	mutate(st_conc=log10(st_conc)) %>% pull(st_conc))

stanMod_1 <- stan(
	file = here('Code','assay_opt_plate.stan'),
	data = stan_data_M1,
	iter = 4000,
	chains = 4)

param <- ss_param_extract(stanMod_1) %>% rownames_to_column('rownames')

samp_point <- cbind(stan_data_M1$S_q, stan_data_M1$Z_qst) %>% as.data.frame() %>% 
	setNames(c("C", "pres")) %>% 
	mutate(!!!setNames(as.list(param$mean), param$rownames))

df <- param %>% 
	filter(grepl("eta_0\\[", rownames)) %>% 
	pull(rownames) %>% 
	imap_dfr(~{
		samp_point %>%
			mutate(Ct = .data[[.x]] + eta_1 * C, plate = gsub("eta_0\\[|\\]", "", .x)) %>%
			mutate(plate = as.numeric(plate)) %>%  # Convert plate to numeric
			select(!!sym(.x), C, pres, eta_1, Ct, plate) %>%
			rename(eta_0 = !!sym(.x)) 
	})

coho_st_plot <- df %>% 
	filter(plate%in%c(4,5)) %>% 
	mutate(plate=as.factor(plate)) %>%
	ggplot() + 
	geom_line(aes(x = 10^C, y = Ct,colour = plate))+
	geom_point(data=qpcr_coho_2_st %>% filter(pres==1) %>% 
						 	filter(plate_idx%in%c(4,5)) %>% 
						 	mutate(plate_idx=as.factor(plate_idx)),
						 aes(x = Quantity, y = Ct,color=plate_idx))+
	scale_x_log10()+
	theme_bw()
coho_st_plot

## qPCR model ----------------------------------------------------------------------------------
# Prepare the data for going into the model
stan_data_M2 <- prep_stan_M2(data = qpcr_coho_2 %>% 
														 	mutate(Task = if_else(Task=='Standard','STANDARD','UNKNOWN')),
														 sample_type = "Task",
														 Ct = "Cq",
														 sample_name_column = "idx",
														 standard_concentration = "Quantity")
# # Run the model in analoge mode
stanMod_2 <- stan(
	file = here('Code','qPCR_mod.stan'),
	data = stan_data_M2,
	iter = 4000,
	chains = 4)

qpcr_coho_2_post <- qpcr_coho_2 %>% left_join(
	extract_param(stanMod_2,'C_q') %>% 
		rownames_to_column('param') %>% 
		mutate(idx=as.numeric(str_extract(param, "\\d+"))) %>% 
		select(-param) %>% 
		rename(Conc='mean'),
	by='idx') %>% 
	filter(idx!=123) %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	select(-Date) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	mutate(time=as.Date(time))

qpcr_coho_2_post %>% 
	# filter(grepl('Air',Sample)) %>% 
	filter(grepl('SAL',Sample)) %>% 
	# filter((Location=='Wallace 156A')) %>%
	# filter(!(Location=='Fisheries pond')) %>%
	# filter(!(grepl('negative',Location))) %>% 
	ggplot()+
	geom_jitter(aes(y=10^(Conc),x=time))+
	facet_grid(FilterType~Location)+
	geom_abline(intercept = -1.5,slope=0)+
	geom_abline(intercept = 1,slope=0,color='red')+
	# facet_grid(~FilterType)+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()


# Chinook -------------------------------------------------------------------------------------

qpcr_chin_1 <- qpcr %>% 
	filter(Sample=='chinook_shel') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame()

qpcr_chin_2 <- qpcr %>% 
	filter(Sample=='chinookDUDA') %>% as.data.frame() %>% 
	rename(Target_2='Sample') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame() %>% 
	left_join(meta %>% filter(plate==1) %>% 
							mutate(plate=8),by=c('Well','plate'))

qpcr_chin_2 %>% filter(Task=='Negative Control')

# qpcr_chin_2 <- qpcr_chin_2 %>% 
# 	filter(!Task=='Negative Control')

#Indexing
qpcr_chin_2 <- qpcr_chin_2 %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>% 
	combine_index(c('j_idx','t_idx','f_idx'),'idx')

qpcr_chin_2_st <- qpcr_chin_2 %>% filter(Task=='Standard') %>% 
	make_index("plate", "plate_idx") %>% 
	mutate(Ct = replace(Cq, Cq == "Undetermined", NA)) %>% 
	mutate(Ct = as.numeric(Ct)) %>% mutate(pres = 1, pres = replace(pres, is.na(Ct), 0)) %>% 
	mutate(st_conc = as.numeric(Quantity))



## Assay efficiency ----------------------------------------------------------------------------
stan_data_M1 <- list(N_st_q = nrow(qpcr_chin_2_st), 
										 N_st_qp = nrow(qpcr_chin_2_st %>% filter(pres==1)), 
										 N_plate = length(unique(qpcr_chin_2_st$plate_idx)),
										 plate_idx = as.numeric(qpcr_chin_2_st %>% filter(pres==1) %>% pull(plate_idx)),
										 Z_qst = as.integer(qpcr_chin_2_st$pres), 
										 S_q = log10(qpcr_chin_2_st$st_conc), 
										 R_qst = as.numeric(qpcr_chin_2_st %>% filter(pres==1) %>% pull(Ct)), 
										 S_q_p = qpcr_chin_2_st %>% filter(pres==1) %>%
										 	mutate(st_conc=log10(st_conc)) %>% pull(st_conc))

stanMod_1 <- stan(
	file = here('Code','assay_opt_plate.stan'),
	data = stan_data_M1,
	iter = 4000,
	chains = 4)

param <- ss_param_extract(stanMod_1) %>% rownames_to_column('rownames')

samp_point <- cbind(stan_data_M1$S_q, stan_data_M1$Z_qst) %>% as.data.frame() %>% 
	setNames(c("C", "pres")) %>% 
	mutate(!!!setNames(as.list(param$mean), param$rownames))

df <- param %>% 
	filter(grepl("eta_0\\[", rownames)) %>% 
	pull(rownames) %>% 
	imap_dfr(~{
		samp_point %>%
			mutate(Ct = .data[[.x]] + eta_1 * C, plate = gsub("eta_0\\[|\\]", "", .x)) %>%
			mutate(plate = as.numeric(plate)) %>%  # Convert plate to numeric
			select(!!sym(.x), C, pres, eta_1, Ct, plate) %>%
			rename(eta_0 = !!sym(.x)) 
	})

chin_st_plot <- df %>% 
	mutate(plate=as.factor(plate)) %>%
	ggplot() + 
	geom_line(aes(x = 10^C, y = Ct,colour = plate))+
	geom_point(data=qpcr_chin_2_st %>% filter(pres==1) %>% 
						 	mutate(plate_idx=as.factor(plate_idx)),
						 aes(x = Quantity, y = Ct,color=plate_idx))+
	scale_x_log10()+
	theme_bw()

chin_st_plot

cowplot::plot_grid(coho_st_plot,chin_st_plot,align = 'h')


## qPCR model ----------------------------------------------------------------------------------
# Prepare the data for going into the model
stan_data_M2 <- prep_stan_M2(data = qpcr_chin_2 %>% 
														 	mutate(Task = if_else(Task=='Standard','STANDARD','UNKNOWN')),
														 sample_type = "Task",
														 Ct = "Cq",
														 sample_name_column = "idx",
														 standard_concentration = "Quantity")

stanMod_2 <- stan(
	file = here('Code','qPCR_mod.stan'),
	data = stan_data_M2,
	iter = 4000,
	chains = 4)


qpcr_chin_2_post <- qpcr_chin_2 %>% left_join(
	extract_param(stanMod_2,'C_q') %>% 
		rownames_to_column('param') %>% 
		mutate(idx=as.numeric(str_extract(param, "\\d+"))) %>% 
		select(-param) %>% 
		rename(Conc='mean'),
	by='idx') %>% 
	filter(idx!=123) %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	select(-Date) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	mutate(time=as.Date(time))

qpcr_chin_2_post %>% 
	# filter(grepl('Air',Sample)) %>% 
	filter(grepl('SAL',Sample)) %>% 
	# filter((Location=='Wallace 156A')) %>%
	# filter(!(Location=='Fisheries pond')) %>%
	# filter(!(grepl('negative',Location))) %>% 
	ggplot()+
	geom_jitter(aes(y=10^(Conc),x=time))+
	facet_grid(FilterType~Location)+
	# facet_grid(~FilterType)+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()


# Corrected annotations

qpcr_st_sockeye <- qpcr %>% 
	filter(grepl('sockeye',Sample)) %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5,5),Quantity*100,Quantity)) %>%
	as.data.frame()
	

# # Prepare the data for going into the model
# stan_data_M1 <- prep_stan_M1(data = qpcr_coho_2_st, 
# 														 Ct = "Cq",
# 														 standard_concentration = "Quantity")

# Prepare the data for going into the model
stan_data_M2 <- prep_stan_M2(data = qpcr_coho_2 %>% 
														 	mutate(Task = if_else(Task=='Standard','STANDARD','UNKNOWN')),
														 sample_type = "Task",
														 Ct = "Cq",
														 sample_name_column = "idx",
														 standard_concentration = "Quantity")

# # Run the model in analoge mode
stanMod_1 <- stan(
	file = here('Code','assay_opt_plate.stan'),
	data = stan_data_M1,
	iter = 4000,
	chains = 4)

ss_param_extract(stanMod_1)

# # Run the model in analoge mode
stanMod_2 <- stan(
	file = here('Code','qPCR_mod.stan'),
	data = stan_data_M2,
	iter = 4000,
	chains = 4)

qpcr_coho_2_post <- qpcr_coho_2 %>% left_join(
	extract_param(stanMod_2,'C_q') %>% 
		rownames_to_column('param') %>% 
		mutate(idx=as.numeric(str_extract(param, "\\d+"))) %>% 
		select(-param) %>% 
		rename(Conc='mean'),
	by='idx') %>% 
	filter(idx!=123) %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	select(-Date) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	mutate(time=as.Date(time))

qpcr_coho_2_post %>% 
	filter(grepl('Air',Sample)) %>% 
	# filter((Location=='Wallace 156A')) %>%
	# filter(!(Location=='Fisheries pond')) %>%
	# filter(!(grepl('negative',Location))) %>% 
	ggplot()+
	geom_jitter(aes(y=10^(Conc),x=time))+
	facet_grid(FilterType~Location)+
	# facet_grid(~FilterType)+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()


# Ready made functions for ploting outputs
plot_ss_param(stan_data=stan_data_M1,ss_param=ss_param)+ggtitle('sample= coho_shel & target= chi1269')
plot_ss_param_pd(stan_data=stan_data_M1,ss_param=ss_param_extract(stanMod_1))+ggtitle('sample= coho_shel & target= chi1269')
plot_ss_param_cm(stan_data=stan_data_M1,ss_param=ss_param_extract(stanMod_1))


# 
# Checks
extract_param(stanMod_1,'sigma_st')
extract_param(stanMod_1,'mu_st')

extract_param(stanMod_1,'alpha_0')
extract_param(stanMod_1,'alpha_1')
extract_param(stanMod_1,'eta_0')
extract_param(stanMod_1,'eta_1')
extract_param(stanMod_1,'gamma_0')
extract_param(stanMod_1,'gamma_1')
extract_param(stanMod_1,'theta_st') %>% mutate(mm=inverselogit(mean))


data_frame <- data.frame(
	S_q_p = stan_data_M1$S_q_p,
	R_qst = stan_data_M1$R_qst
)

matrix(data = rep(-1:6, each = 100), ncol = 1, nrow = 800) %>% 
	as.data.frame() %>% 
	setNames('C') %>% 
	mutate(
		gamma_0 = extract_param(stanMod_1, 'gamma_0')$mean,
		gamma_1 = extract_param(stanMod_1, 'gamma_1')$mean,
		eta_0 = extract_param(stanMod_1, 'eta_0')$mean,
		eta_1 = extract_param(stanMod_1, 'eta_1')$mean,
		sigma = exp(gamma_0 + (gamma_1 * C)),
		mu_st = eta_0 + (eta_1 * C)) %>% 
	rowwise() %>% 
	mutate(Ct = rnorm(1, mu_st, sigma)) %>% 
	ungroup() %>% as.data.frame() %>% 
	ggplot()+
	geom_point(aes(y=Ct,x=C))+
	geom_point(data=data_frame,aes(y=R_qst,x=S_q_p+0.1),color='red')+
	ylim(18,47)

matrix(data = rep(-1:6, each = 100), ncol = 1, nrow = 800) %>% 
	as.data.frame() %>% 
	setNames('C') %>% 
	mutate(
		gamma_0 = extract_param(stanMod_1, 'gamma_0')$mean,
		gamma_1 = extract_param(stanMod_1, 'gamma_1')$mean,
		eta_0 = extract_param(stanMod_1, 'eta_0')$mean,
		eta_1 = extract_param(stanMod_1, 'eta_1')$mean,
		sigma = exp(gamma_0 + (gamma_1 * C)),
		mu_st = eta_0 + (eta_1 * C)) %>% 
	rowwise() %>% 
	mutate(Ct = rnorm(1, mu_st, sigma)) %>% 
	ungroup() %>% as.data.frame() %>% 
	ggplot()+
	geom_point(aes(y=log(sigma),x=C))


pd <- data.frame(
	Z_qst = stan_data_M1$Z_qst,
	S_q = stan_data_M1$S_q
)
pd %>% group_by(S_q) %>% 
	summarise(pd=mean(Z_qst))
