# Libraries -----------------------------------------------------------------------------------

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


# Functions -----------------------------------------------------------------------------------
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


# Load data -----------------------------------------------------------------------------------

file_paths <- list.files('Data', pattern = '\\.csv$', full.names = TRUE)
file_paths_meta <- file_paths[grepl('metadata',file_paths)]
file_paths_count <- file_paths[grepl('count',file_paths)]
file_paths_data <- file_paths[!(file_paths%in%file_paths_meta|file_paths%in%file_paths_count)]

coun <- lapply(seq_along(file_paths_count), function(idx) {
	read_csv(file_paths_count[idx]) %>% mutate(target = idx)}) %>%
	bind_rows()

meta <- lapply(seq_along(file_paths_meta), function(idx) {
	read_csv(file_paths_meta[idx]) %>% mutate(csv = idx)}) %>%
	bind_rows()

qpcr <- lapply(seq_along(file_paths_data), function(idx) {
	read_csv(file_paths_data[idx]) %>% mutate(plate = idx)}) %>%
	bind_rows()

# Rename the targets to the same target
qpcr <- qpcr %>% 
	mutate(Sample=ifelse(Sample=='cohoSHEL','coho_shel',Sample)) %>% 
	mutate(Sample=ifelse(Sample=='chinookSHEL','chinook_shel',Sample))

# Checking the different target names
qpcr %>% filter(Task=='Standard') %>% 
	distinct(plate,Sample,Target) %>% 
	count(Sample,Target)

qpcr_coho <- qpcr %>% 
	filter(Sample=='cohoDUDA') %>% 
	rename(Target_2='Sample') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame() %>% 
	rename(Plate_idx=plate) %>% 
	left_join(meta %>% filter(csv==1),by=c('Well','Plate_idx'))

# Indexing
qpcr_coho <- qpcr_coho %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>% 
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>% 
	make_index("Plate_idx", "plate_idx")

qpcr_coho_st <- qpcr_coho %>% filter(Task=='Standard') %>% 
	mutate(Ct = replace(Cq, Cq == "Undetermined", NA)) %>% 
	mutate(Ct = as.numeric(Ct)) %>% mutate(pres = 1, pres = replace(pres, is.na(Ct), 0)) %>% 
	mutate(st_conc = as.numeric(Quantity))

qpcr_coho_en <- qpcr_coho %>% filter(Task=='Unknown') %>% 
	mutate(Ct = replace(Cq, Cq == "Undetermined", NA)) %>% 
	mutate(Ct = as.numeric(Ct)) %>% mutate(pres = 1, pres = replace(pres, is.na(Ct), 0)) %>% 
	mutate(st_conc = as.numeric(Quantity))

stan_data_M1 <- list(N_st_q = nrow(qpcr_coho_st), 
										 N_st_qp = nrow(qpcr_coho_st %>% filter(pres==1)), 
										 N_plate = length(unique(qpcr_coho_st$plate_idx)),
										 plate_idx = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(plate_idx)),
										 Z_qst = as.integer(qpcr_coho_st$pres), 
										 S_q = log10(qpcr_coho_st$st_conc), 
										 R_qst = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(Ct)), 
										 S_q_p = qpcr_coho_st %>% filter(pres==1) %>%
										 	mutate(st_conc=log10(st_conc)) %>% pull(st_conc))

stan_data_M2 <- list(N_st_q = nrow(qpcr_coho_st), 
										 N_en_q = nrow(qpcr_coho_en),
										 N_st_qp = nrow(qpcr_coho_st %>% filter(pres==1)), 
										 N_en_qp = nrow(qpcr_coho_en %>% filter(pres==1)),
										 N_plate = length(unique(qpcr_coho_st$plate_idx)),
										 N_j = length(unique(qpcr_coho_en$idx)),
										 plate_st_idx = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(plate_idx)),
										 plate_en_idx = as.numeric(qpcr_coho_en %>% filter(pres==1) %>% pull(plate_idx)),
										 j_qen_idx = qpcr_coho_en$idx, 
										 j_qen_p_idx = qpcr_coho_en %>% filter(pres==1) %>% pull(idx), 
										 Z_qst = as.integer(qpcr_coho_st$pres), 
										 Z_qen = as.integer(qpcr_coho_en$pres), 
										 S_q = log10(qpcr_coho_st$st_conc), 
										 R_qst = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(Ct)), 
										 R_qen = as.numeric(qpcr_coho_en %>% filter(pres==1) %>% pull(Ct)),
										 S_q_p = qpcr_coho_st %>% filter(pres==1) %>%
										
										 	mutate(st_conc=log10(st_conc)) %>% pull(st_conc),
										label_M2 = qpcr_coho_en %>% distinct(idx, Sample) %>% 
										
											arrange(idx) %>% as.data.frame()
										)

stanMod_2 <- stan(
	file = here('Code','qPCR_mod_plate.stan'),
	data = stan_data_M2,
	iter = 4000,
	chains = 4)

param <- ss_param_extract(stanMod_2) %>% rownames_to_column('rownames')

samp_point <- cbind(stan_data_M2$S_q, stan_data_M2$Z_qst) %>% as.data.frame() %>% 
	setNames(c("C", "pres")) %>% 
	mutate(!!!setNames(as.list(param$mean), param$rownames))

# Standard curves
param %>% 
	filter(grepl("eta_0\\[", rownames)) %>% 
	pull(rownames) %>% 
	imap_dfr(~{
		samp_point %>%
			mutate(Ct = .data[[.x]] + eta_1 * C, plate = gsub("eta_0\\[|\\]", "", .x)) %>%
			mutate(plate = as.numeric(plate)) %>%  # Convert plate to numeric
			select(!!sym(.x), C, pres, eta_1, Ct, plate) %>%
			rename(eta_0 = !!sym(.x)) }) %>% 
	filter(plate%in%c(1:5)) %>% 
	mutate(plate=as.factor(plate)) %>%
	ggplot() + 
	geom_line(aes(x = 10^C, y = Ct,colour = plate))+
	geom_point(data=qpcr_coho_st %>% filter(pres==1) %>% 
						 	filter(plate_idx%in%c(4,5)) %>% 
						 	mutate(plate_idx=as.factor(plate_idx)),
						 aes(x = Quantity, y = Ct,color=plate_idx))+
	scale_x_log10()+
	theme_bw()

# Create posterior dataframe
qpcr_coho_post <- qpcr_coho_en %>% left_join(
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

count_fish <- coun %>% filter(Target=='Coho') %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	select(-Date) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	mutate(time=as.Date(time)) %>% 
	mutate(Location='Issaquah waterfall') %>% 
	mutate(FilterType='Observation')
	
qpcr_coho_post %>% 
	filter(time> as.Date('2024-10-1')) %>%
	filter(grepl('Air',Sample)) %>%
	# filter(grepl('SAL',Sample)) %>% 
	filter(!(Location=='Issaquah park bridge')) %>%
	filter(!(Location=='Wallace 156A')) %>%
	filter(!(Location=='Fisheries pond')) %>%
	# filter(!(grepl('negative',Location))) %>% 
	ggplot()+
	geom_point(aes(y=10^(Conc),x=time))+
	# geom_errorbar(aes(ymin=10^(`2.5%`),ymax=10^(`97.5%`),x=time), width = 0.4)+
	facet_wrap(FilterType~Location,ncol=1)+
	# geom_abline(intercept = -1.5,slope=0)+
	# geom_abline(intercept = 1,slope=0,color='red')+
	# facet_grid(~FilterType)+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4),
								limits=c(10^-2, 10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()



qpcr_coho_post %>% arrange(time) %>% 
	filter(time> as.Date('2024-10-1')) %>% 
	filter(grepl('Air',Sample)) %>%
	filter(!(Location=='Wallace 156A')) %>%
	filter(!(Location=='Fisheries pond')) %>%
	group_by(Location,FilterType) %>% 
	mutate(CS=cumsum(10^(Conc)),
				 TotalSum = sum(10^(Conc)), 
				 ProportionalCS = CS / TotalSum) %>% 
	ungroup() %>% 
	as.data.frame() %>% 
	ggplot()+
	geom_line(aes(y=ProportionalCS,x=time))+
	geom_line(data=x,aes(y=ProportionalCS,x=time),color='grey50',lty=2)+
	facet_grid(FilterType~Location)+
	ggtitle('Air samples')+
	theme_bw()


x <- expand_grid(count_fish,
	new_loca = unique(qpcr_coho_post$Location)[c(1,6)]) %>%
	mutate(Location = new_loca) %>%
	select(-new_loca) %>% 
expand_grid(.,
	new_filt = unique(qpcr_coho_post$FilterType)[c(5,4,3,1)]) %>%
	mutate(Obs = new_filt) %>%
	select(-new_filt,-FilterType) %>% 
	arrange(time) %>% 
	filter(!(Location=='Wallace 156A')) %>%
	filter(!(Location=='Fisheries pond')) %>%
	group_by(Location,Obs) %>% 
	mutate(CS=cumsum((Total)),
				 TotalSum = sum((Total)), 
				 ProportionalCS = CS / TotalSum) %>% 
	ungroup()
	
ggplot()+
	geom_line(data=x,aes(y=ProportionalCS,x=time))+
	facet_wrap(Location~Obs)+
	ggtitle('Air samples')+
	theme_bw()	


	ggplot()+
	geom_point(aes(y=(Total),x=time))+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4),
								limits=c(10^-2, 10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()

count_fish <- count_fish %>% 
		mutate(cum_count=cumsum(Total))

standata_count <- list(
	time=length(count_fish$cum_count),
	y_cum_sum=count_fish$cum_count,
	y=count_fish$Total
)

stanMod_count <- stan(
	file = here('Code','Count_model.stan'),
	data = standata_count,
	iter = 4000,
	chains = 4)
# 
# count_fish %>% pull(Total) %>% barplot()
# extract_param(stanMod_count,'lambda_s') %>% pull(mean) %>% barplot()
# extract_param(stanMod_count,'y_step') %>% pull(mean) %>% barplot()

ggplot(count_fish) + 
	geom_bar(stat = "identity",aes(y=Total, x=time))

n=nrow(count_fish)
custom_matrix <-
	diag(1, n, n) + 
	# rbind(rep(0,n),diag(1, n - 1, n))+
	rbind(diag(1, n, n)[-1,],rep(0,n))+
	rbind(diag(1, n, n)[-(1:2),],rep(0,n),rep(0,n))
	
custom_matrix <- custom_matrix %>% as.data.frame() %>% 
	mutate_all(~ replace(., . == 0, NA))

count_fish_matrix <- count_fish$Total*custom_matrix 

count_fish_matrix_lf <- count_fish_matrix %>% 
	pivot_longer(names_to = 'idx',
							 cols = everything(),
							 values_to = 'count') %>% 
	filter(!is.na(count)) %>% 
	mutate(idx=substr(idx,2,nchar(.))) %>% 
	mutate(idx=as.numeric(idx)) %>% 
	arrange(idx) %>% as.data.frame()


standata_count_2 <- list(
	time=length(count_fish$cum_count),
	N=length(count_fish_matrix_lf$idx),
	idx=count_fish_matrix_lf$idx,
	y=count_fish_matrix_lf$count)

stanMod_count <- stan(
	file = here('Code','Count_model_2.stan'),
	data = standata_count_2,
	iter = 4000,
	chains = 4)
# 
count_fish %>% pull(Total) %>% barplot()

count_fish %>% 
	cbind(.,extract_param(stanMod_count_2,'lambda') %>% select(mean) %>% rename(count_smooth_2='mean')) %>% 
	cbind(.,extract_param(stanMod_count_3,'lambda') %>% select(mean) %>% rename(count_smooth_3='mean')) %>% 
	ggplot() + 
	# geom_bar(stat = "identity",aes(y=Total, x=time))+
	geom_smooth(aes(y=Total, x=time),se=F,color='red')+
	geom_smooth(aes(y=count_smooth_2, x=time),se=F,color='orange',lty=2)+
	geom_smooth(aes(y=count_smooth_3, x=time),se=F,color='blue',lty=2)+
	geom_point(data=qpcr_coho_post %>% filter(time> as.Date('2024-10-1')) %>% 
							filter(grepl('Air',Sample)) %>%
							# filter(grepl('SAL',Sample)) %>% 
							filter(!(Location=='Issaquah park bridge')) %>%
							filter(!(Location=='Wallace 156A')) %>%
							filter(!(Location=='Fisheries pond')) %>%
						 	filter(FilterType=='MCE DI water'),
						 	aes(x=time,y=(10^(Conc))*0.3))+
						 	# filter(FilterType=='Gelatin'),
						 	# aes(x=time,y=(10^(Conc))*15))+
						 	# filter(FilterType=='PTFE'),
						 	# aes(x=time,y=(10^(Conc))*5))+
						 	# filter(FilterType=='MCE Air'),
	theme_bw()
# 

qpcr_coho_post %>% filter(time> as.Date('2024-10-1')) %>% 
	filter(grepl('Air',Sample)) %>%
	# filter(grepl('SAL',Sample)) %>% 
	filter(!(Location=='Issaquah park bridge')) %>%
	filter(!(Location=='Wallace 156A')) %>%
	filter(!(Location=='Fisheries pond')) %>%
	# filter(FilterType=='MCE DI water') %>% 
	filter(FilterType=='Gelatin') %>% 
	# filter(FilterType=='PTFE') %>% 
	# filter(FilterType=='MCE Air') %>% 
	
	ggplot()+
	geom_point(aes(y=10^(Conc),x=time))+
	geom_errorbar(aes(ymin=10^(`2.5%`),ymax=10^(`97.5%`),x=time), width = 0.4)+
	facet_wrap(FilterType~Location,ncol=2)+
	# geom_abline(intercept = -1.5,slope=0)+
	# geom_abline(intercept = 1,slope=0,color='red')+
	# facet_grid(~FilterType)+
	scale_y_log10(labels=scientific,breaks=c(10^-2,10^0,10^2,10^4),
								limits=c(10^-2, 10^4))+
	ylab(expression('Log'[10]*' DNA concentration'))+
	ggtitle('Air samples')+
	theme_bw()
	
extract_param(stanMod_count_2,'lambda') %>% pull(mean) %>% barplot()
extract_param(stanMod_count_3,'lambda') %>% pull(mean) %>% barplot()

count_fish_effort <- count_fish %>% 
	mutate(Total=if_else(time=='2024-12-12',1,Total)) %>% 
	filter(Total>0) %>% 
	mutate(time_days = as.numeric(time - lag(time))) %>% 
	mutate(time_days=if_else(is.na(time_days),1,time_days)) %>% 
	mutate(Total=if_else(time=='2024-12-12',0,Total))

standata_count <- list(
	time=nrow(count_fish_effort),
	days=count_fish_effort$time_days,
	lambda_p = c(1,0.1),
	lambda_rp = c(1,0.1),
	beta_p = c(1,1),
	epsilon_p = c(0,1),
	y=count_fish_effort$Total
)

stanMod_count <- stan(
	file = here('Code','Count_model_5.stan'),
	data = standata_count,
	iter = 4000,
	control=list(max_treedepth=12),
	chains = 4)
stanMod_count
# 


n=nrow(count_fish_effort)
custom_matrix <-
	diag(1, n, n) + 
	# rbind(rep(0,n),diag(1, n - 1, n))+
	rbind(diag(1, n, n)[-1,],rep(0,n))

custom_matrix <- custom_matrix %>% as.data.frame() %>% 
	mutate_all(~ replace(., . == 0, NA))

count_fish_matrix <- count_fish_effort$Total*custom_matrix 
count_fish_matrix_days <- count_fish_effort$time_days*custom_matrix

count_fish_matrix_lf <- count_fish_matrix %>% 
	pivot_longer(names_to = 'idx',
							 cols = everything(),
							 values_to = 'count') %>% 
	filter(!is.na(count)) %>% 
	mutate(idx=substr(idx,2,nchar(.))) %>% 
	mutate(idx=as.numeric(idx)) %>% 
	arrange(idx) %>% as.data.frame()

count_fish_matrix_days_lf <- count_fish_matrix_days %>% 
	pivot_longer(names_to = 'idx',
							 cols = everything(),
							 values_to = 'days') %>% 
	filter(!is.na(days)) %>% 
	mutate(idx=substr(idx,2,nchar(.))) %>% 
	mutate(idx=as.numeric(idx)) %>% 
	arrange(idx) %>% as.data.frame()

count_fish_effort_lf_2 <- 
	count_fish_matrix_lf %>% 
	cbind(.,count_fish_matrix_days_lf %>% select(-idx))

standata_count <- list(
	time=max(count_fish_effort_lf_2$idx),
	N=nrow(count_fish_effort_lf_2),
	days=count_fish_effort_lf_2$days,
	idx=count_fish_effort_lf_2$idx,
	y=count_fish_effort_lf_2$count)

stanMod_count <- stan(
	file = here('Code','Count_model_4.stan'),
	data = standata_count,
	iter = 4000,
	# control=list(max_treedepth=12),
	chains = 4)
stanMod_count
# 



