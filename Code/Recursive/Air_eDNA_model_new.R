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

join_ext_param <- function(stanmod,par){
	extract_param(stanmod,par) %>% 
	rownames_to_column('param') %>% 
	mutate(g_idx = str_extract(param, "\\d+")) %>% 
	mutate(g_idx = as.numeric(g_idx))
}

# Load data -----------------------------------------------------------------------------------

# file_paths <- list.files('Data', pattern = '\\.csv$', full.names = TRUE)
# file_paths_meta <- file_paths[grepl('metadata',file_paths)]
# file_paths_count <- file_paths[grepl('count',file_paths)]
# file_paths_data <- file_paths[!(file_paths%in%file_paths_meta|file_paths%in%file_paths_count)]
meta_data_files <- here('Data','Metadata') %>% list.files()
file_paths_meta <- here('Data','Metadata',meta_data_files)

count_data_files <- here('Data','Count_data') %>% list.files()
file_paths_count <- here('Data','Count_data',count_data_files)

qPCR_data_files <- here('Data','qPCR_data') %>% list.files()
file_paths_data <- here('Data','qPCR_data',qPCR_data_files)

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
	select(-Target) %>% 
	rename(Target='Sample') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame() %>% 
	rename(Plate_idx=plate) %>% 
	left_join(meta,by=c('Well','Plate_idx','Target'))

exclude_chin_duda <- c('A13','A14','B13','B14','C13','C14','D13','D14','E13','E14','F13','F14','G13','G14','H13','H14','I13','J13','K13','L13','M13','N13','O13','P13')
qpcr_chin <- qpcr %>% 
	filter(Sample=='chinookDUDA') %>% 
	select(-Target) %>% 
	rename(Target='Sample') %>% 
	mutate(Quantity=if_else(plate%in%c(3,4,5),Quantity*100,Quantity)) %>%
	as.data.frame() %>% 
	rename(Plate_idx=plate) %>% 
	left_join(meta,by=c('Well','Plate_idx','Target')) %>% 
	filter(!(Well%in%exclude_chin_duda&
					 	Plate_idx==3&
					 	Target=='chinookDUDA'&
					 	Task=='Unknown'))

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

count_fish <- coun %>% filter(Target=='Coho') %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	select(-Date) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	mutate(time=as.Date(time)) %>% 
	mutate(Location='Issaquah waterfall') %>% 
	mutate(FilterType='Observation')

count_fish_effort <- count_fish %>% 
	mutate(Total=if_else(time=='2024-12-12',1,Total)) %>% 
	filter(Total>0) %>% 
	mutate(time_days = as.numeric(time - lag(time))) %>% 
	mutate(time_days=if_else(is.na(time_days),1,time_days)) %>% 
	mutate(Total=if_else(time=='2024-12-12',0,Total))


qpcr_coho_wat_en <- qpcr_coho_en %>% 
	filter(Location!='Issaquah park bridge') %>% 
	filter(Location!='Wallace 156A') %>% 
	filter(Location!='Fisheries pond') %>% 
	filter(!grepl('field negative',Location)) %>% 
	filter(Location!='Issaquah Bridge, Confluence Park, by hand') %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	filter(time>as.Date('2024-10-1')) %>% as.data.frame() %>% 
	filter(FilterType=='Water filter MCE') %>% 
	arrange(time) %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>%
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>% arrange(j_idx) %>% 
	as.data.frame()

qpcr_coho_en_prep <-
	qpcr_coho_en %>% 
	filter(Location!='Issaquah park bridge') %>% 
	filter(Location!='Wallace 156A') %>% 
	filter(Location!='Fisheries pond') %>% 
	filter(!grepl('field negative',Location)) %>% 
	mutate(Year=substr(Date,0,4)) %>% 
	mutate(Month=substr(Date,5,6)) %>% 
	mutate(Day=substr(Date,7,8)) %>% 
	mutate(time=paste(Year,Month,Day,sep = '-')) %>% 
	filter(time>as.Date('2024-10-1')) %>% 
	filter(FilterType!='Water filter MCE')
  
	### Filtering only the good Gelatin filters ###
	#filter(!(Well%in%c('E10','E11','E12')&FilterType=='Gelatin'&time=='2024-10-24')) %>%
	#filter(!(Well%in%c('M4','H16','H17','H18')&FilterType=='Gelatin'&time=='2024-11-07')) %>%
	#filter(!(Well%in%c('H22','H23','H24')&FilterType=='Gelatin'&time=='2024-11-21'))

qpcr_coho_air_en_a <- qpcr_coho_en_prep %>% filter(FilterType=='Gelatin') %>% 
	arrange(time) %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>%
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>%
	combine_index(c('t_idx','f_idx'),'ft_idx') %>% 
	arrange(j_idx) %>% 
	as.data.frame()
qpcr_coho_air_en_b <- qpcr_coho_en_prep %>% filter(FilterType=='PTFE') %>% 
	arrange(time) %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>%
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>%
	combine_index(c('t_idx','f_idx'),'ft_idx') %>% 
	arrange(j_idx) %>% 
	as.data.frame()
qpcr_coho_air_en_c <- qpcr_coho_en_prep %>% filter(FilterType=='MCE Air') %>% 
	arrange(time) %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>%
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>%
	combine_index(c('t_idx','f_idx'),'ft_idx') %>% 
	arrange(j_idx) %>% 
	as.data.frame()
qpcr_coho_air_en_d <- qpcr_coho_en_prep %>% filter(FilterType=='MCE DI water') %>% 
	arrange(time) %>% 
	make_index("Sample", "j_idx") %>% 
	make_index("Date", "t_idx") %>% 
	make_index("FilterType", "f_idx") %>%
	combine_index(c('j_idx','t_idx','f_idx'),'idx') %>%
	combine_index(c('t_idx','f_idx'),'ft_idx') %>% 
	arrange(j_idx) %>% 
	as.data.frame()

count_fish_effort_selected <- 
	count_fish_effort %>% filter(time%in%unique(qpcr_coho_en_prep$time))


for (i in 1:4){
if(i==1){qpcr_coho_air_en <- qpcr_coho_air_en_a}
else if(i == 2){qpcr_coho_air_en <- qpcr_coho_air_en_b}
else if(i == 3){qpcr_coho_air_en <- qpcr_coho_air_en_c}
else if(i == 4){qpcr_coho_air_en <- qpcr_coho_air_en_d}

stan_data <- list(N_st_q = nrow(qpcr_coho_st),
									N_en_wat_q = nrow(qpcr_coho_wat_en),
									N_en_air_q = nrow(qpcr_coho_air_en),
									# 
									N_st_qp = nrow(qpcr_coho_st %>% filter(pres==1)),
									N_en_wat_qp = nrow(qpcr_coho_wat_en %>% filter(pres==1)),
									N_en_air_qp = nrow(qpcr_coho_air_en %>% filter(pres==1)),
									# 
									N_filt = length(unique(qpcr_coho_air_en$FilterType)),
									N_plate = length(unique(qpcr_coho_st$plate_idx)),
									N_j_wat = length(unique(qpcr_coho_wat_en$idx)),
									N_j_air = length(unique(qpcr_coho_air_en$idx)),
									N_ft_air = qpcr_coho_air_en %>% distinct(FilterType,time,f_idx) %>% nrow(),
									# 
									plate_st_idx = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(plate_idx)),
									plate_en_wat_idx = as.numeric(qpcr_coho_wat_en %>% filter(pres==1) %>% pull(plate_idx)),
									plate_en_air_idx = as.numeric(qpcr_coho_air_en %>% filter(pres==1) %>% pull(plate_idx)),
									# 
									j_qen_wat_idx = qpcr_coho_wat_en$idx,
									j_qen_air_idx = qpcr_coho_air_en$idx,
									j_qen_wat_p_idx = qpcr_coho_wat_en %>% filter(pres==1) %>% pull(idx),
									j_qen_air_p_idx = qpcr_coho_air_en %>% filter(pres==1) %>% pull(idx),
									# 
									Z_qst = as.integer(qpcr_coho_st$pres),
									Z_qen_wat = as.integer(qpcr_coho_wat_en$pres),
									Z_qen_air = as.integer(qpcr_coho_air_en$pres),
									# 
									R_qst = as.numeric(qpcr_coho_st %>% filter(pres==1) %>% pull(Ct)),
									R_qen_wat = as.numeric(qpcr_coho_wat_en %>% filter(pres==1) %>% pull(Ct)),
									R_qen_air = as.numeric(qpcr_coho_air_en %>% filter(pres==1) %>% pull(Ct)),
									# 
									S_q = log(qpcr_coho_st$st_conc),
									S_q_p = qpcr_coho_st %>% filter(pres==1) %>% mutate(st_conc=log(st_conc)) %>% pull(st_conc),
									# 
									# t1_idx = c(3,3,4,5,6,6),
									t1_idx = c(1,1,2,3,4,4),
									# t0_idx = c(2,3,3,4,5,6),
									# air_wat_idx = qpcr_coho_air_en %>% distinct(time,idx,t_idx) %>% pull(t_idx),
									# air_filt_idx = qpcr_coho_air_en %>% distinct(FilterType,time,idx,f_idx) %>% pull(f_idx),
									ft_idx = qpcr_coho_air_en %>% distinct(FilterType,time,idx,ft_idx) %>% pull(ft_idx),
									air_wat_idx2 = qpcr_coho_air_en %>% distinct(FilterType,time,f_idx,t_idx,ft_idx) %>% pull(t_idx),
									air_filt_idx2 = qpcr_coho_air_en %>% distinct(FilterType,time,f_idx,t_idx,ft_idx) %>% pull(f_idx),
									ft_idx2 = qpcr_coho_air_en %>% distinct(FilterType,time,f_idx,t_idx,ft_idx) %>% pull(ft_idx),
									# 
									# time = nrow(count_fish_effort_selected),
									time = 4,
									tau_p1 = 0,tau_p2=1,
									N = count_fish_effort_selected$Total,
									E = count_fish_effort_selected$time_days,
									logit_phi_mu = -3,
									logit_phi_sd = 1,
									label_filter = unique(qpcr_coho_air_en$FilterType),
									label_wat = qpcr_coho_wat_en %>% distinct(idx, Sample) %>% arrange(idx) %>% as.data.frame(),
									label_air = qpcr_coho_air_en %>% distinct(idx, Sample) %>% arrange(idx) %>% as.data.frame()
)

if(i==1){stan_data_a <- stan_data}
else if(i == 2){stan_data_b <- stan_data}
else if(i == 3){stan_data_c <- stan_data}
else if(i == 4){stan_data_d <- stan_data}
}



for (i in 1:4){
	if(i==1){stan_data <- stan_data_a}
	else if(i == 2){stan_data <- stan_data_b}
	else if(i == 3){stan_data <- stan_data_c}
	else if(i == 4){stan_data <- stan_data_d}

stanMod_count <- stan(
	file = here('Code','Count_model_5.stan'),
	data = stan_data,
	iter = 4000,
	control=list(max_treedepth=12),
	chains = 4)

if(i==1){stanMod_count_a <- stanMod_count}
else if(i == 2){stanMod_count_b <- stanMod_count}
else if(i == 3){stanMod_count_c <- stanMod_count}
else if(i == 4){stanMod_count_d <- stanMod_count}
}


table_1 <- 
join_ext_param(stanMod_count_a,'eta') %>% mutate(FT='gelatin') %>% rbind(.,
join_ext_param(stanMod_count_b,'eta') %>% mutate(FT='PTFE')) %>% rbind(.,
join_ext_param(stanMod_count_c,'eta') %>% mutate(FT='MCE Air')) %>% rbind(.,
join_ext_param(stanMod_count_d,'eta') %>% mutate(FT='MCE DI water')) %>% 
	rename(dilution='mean') %>% select(FT,dilution) %>% 
left_join(.,

join_ext_param(stanMod_count_a,'epsilon') %>% mutate(FT='gelatin') %>% rbind(.,
join_ext_param(stanMod_count_b,'epsilon') %>% mutate(FT='PTFE')) %>% rbind(.,
join_ext_param(stanMod_count_c,'epsilon') %>% mutate(FT='MCE Air')) %>% rbind(.,
join_ext_param(stanMod_count_d,'epsilon') %>% mutate(FT='MCE DI water')) %>% 
	rename(epsilon='mean') %>% select(FT,epsilon) %>% 
	group_by(FT) %>% 
	summarise(SEE=round(mean(abs(epsilon)),3)),
	# summarise(SEE=sqrt(sum(epsilon^2))),
by='FT') %>% 
left_join(.,

join_ext_param(stanMod_count_a,'delta') %>% mutate(FT='gelatin') %>% rbind(.,
join_ext_param(stanMod_count_b,'delta') %>% mutate(FT='PTFE')) %>% 
	rename(delta='mean') %>% select(FT,delta) %>% 
	group_by(FT) %>% 
	summarise(bio_rep=round(mean(abs(delta)),3)),
by='FT')

table_1


log_W <- join_ext_param(stanMod_count_a,'log_W') %>% 
	setNames(paste0('log_W_',colnames(.))) %>% 
left_join(.,
qpcr_coho_wat_en %>% select(Sample,idx,plate_idx,time),by=c('log_W_g_idx'='idx')) %>% 
	mutate(time=if_else(time=='2024-10-25','2024-10-24',time))
																																				 
log_A <- 
join_ext_param(stanMod_count_a,'log_A') %>% mutate(FT='Gelatin') %>% rbind(.,
join_ext_param(stanMod_count_b,'log_A') %>% mutate(FT='PTFE')) %>% rbind(.,
join_ext_param(stanMod_count_c,'log_A') %>% mutate(FT='MCE Air')) %>% rbind(.,
join_ext_param(stanMod_count_d,'log_A') %>% mutate(FT='MCE DI water')) %>% 
	setNames(paste0('log_A_',colnames(.))) %>% 
left_join(.,
rbind(qpcr_coho_air_en_a,qpcr_coho_air_en_b,qpcr_coho_air_en_c,qpcr_coho_air_en_d) %>% 
	distinct(Sample,j_idx,t_idx,idx,plate_idx,time,FilterType),
by=c("log_A_g_idx" = "j_idx","log_A_FT" = "FilterType"))
																																				 
mean_log_A <-
join_ext_param(stanMod_count_a,'mean_log_A') %>% mutate(FT='Gelatin') %>% rbind(.,
join_ext_param(stanMod_count_b,'mean_log_A') %>% mutate(FT='PTFE')) %>% rbind(.,
join_ext_param(stanMod_count_c,'mean_log_A') %>% mutate(FT='MCE Air')) %>% rbind(.,
join_ext_param(stanMod_count_d,'mean_log_A') %>% mutate(FT='MCE DI water')) %>% 
	setNames(paste0('mean_log_A_',colnames(.))) %>% 
left_join(.,
rbind(qpcr_coho_air_en_a,qpcr_coho_air_en_b,qpcr_coho_air_en_c,qpcr_coho_air_en_d) %>% 
	distinct(t_idx,time,FilterType),
by=c("mean_log_A_g_idx" = "t_idx","mean_log_A_FT" = "FilterType"))

N <-
count_fish_effort %>%
	mutate(idx=row_number()-2) %>% 
left_join(.,
join_ext_param(stanMod_count_a,'lambda') %>% 
	setNames(paste0('lambda_',colnames(.))),
by=c('idx'='lambda_g_idx')) %>% 
	mutate(lambda_mean=if_else(idx==0,Total/time_days,lambda_mean)) %>% 
filter(time>=as.Date('2024-10-10')&idx<=4)

omega <- extract_param(stanMod_count_a,'omega') %>% pull(mean)

p1 <-
ggplot()+
	geom_smooth(data=N,aes(x=time,y=log(lambda_mean)-omega),col='black',span=0.5)+
	geom_point(data=N,aes(x=time,y=log(lambda_mean)-omega),col='black',size=3)+
	geom_smooth(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',se=F,span=0.5)+
	geom_point(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',size=3)+
	scale_y_continuous(
		limits = c(8, 14),
		bquote(log[e]('eDNA Conc') ~ "from water samples"),
		sec.axis = sec_axis(~ . + omega, name = bquote(log[e]('Fish Density'))))+
	theme_bw() +
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"),  # Custom breaks every 2 weeks
		labels = date_format("%b %d")
	) +
	theme(
		axis.text = element_text(size = 14),
		axis.title.x = element_blank(),
		legend.text = element_text(size = 14),  
		legend.title = element_text(size = 15, face = "bold"),
		
		# Left axis
		axis.ticks.y.left = element_line(color = "#61BEA4", size = 1.2),
		axis.title.y.left = element_text(color = "#61BEA4", size = 14),
		axis.text.y.left = element_text(color = "#61BEA4", size = 13),
		
		# Right axis
		axis.ticks.y.right = element_line(color = "black", size = 1.2),
		axis.title.y.right = element_text(color = "black", size = 14),
		axis.text.y.right = element_text(color = "black", size = 13)
	)


p2 <-
	ggplot()+
	geom_smooth(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',se=F,span=0.5)+
	geom_point(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',size=3)+
	geom_point(data=log_A %>% filter(log_A_FT=='Gelatin'),
							aes(x=as.Date(time),y=table_1$dilution[1]+log_A_mean),col='#AB5971',size=6,pch=3)+
	geom_smooth(data=mean_log_A %>% filter(mean_log_A_FT=='Gelatin'),
							aes(x=as.Date(time),y=table_1$dilution[1]+mean_log_A_mean),col='#AB5971')+
	scale_y_continuous(
		limits = c(8, 14),
		bquote(log[e]('eDNA Conc') ~ " water samples"),
		sec.axis = sec_axis(~ . - table_1$dilution[1], name = bquote(log[e]('eDNA Conc') ~ "air samples - Gelatin filters")))+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"),  # Custom breaks every 2 weeks
		labels = date_format("%b %d")
	) +	
	theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title.x = element_blank(),
			legend.text = element_text(size = 14),  
			legend.title = element_text(size = 15, face = "bold"),
			
			# Left axis
			axis.ticks.y.left = element_line(color = "#61BEA4", size = 1.2),
			axis.title.y.left = element_text(color = "#61BEA4", size = 14),
			axis.text.y.left = element_text(color = "#61BEA4", size = 13),
			
			# Right axis
			axis.ticks.y.right = element_line(color = "#AB5971", size = 1.2),
			axis.title.y.right = element_text(color = "#AB5971", size = 14),
			axis.text.y.right = element_text(color = "#AB5971", size = 13)
		)

p3 <-
	ggplot()+
	geom_smooth(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',se=F,span=0.5)+
	geom_point(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',size=3)+
	geom_point(data=log_A %>% filter(log_A_FT=='PTFE'),
							aes(x=as.Date(time),y=table_1$dilution[2]+log_A_mean),col='#F49D4D',size=6,pch=3)+
	geom_smooth(data=mean_log_A %>% filter(mean_log_A_FT=='PTFE'),
							aes(x=as.Date(time),y=table_1$dilution[2]+mean_log_A_mean),col='#F49D4D')+
	scale_y_continuous(
		limits = c(8, 14),
		bquote(log[e]('eDNA Conc') ~ " water samples"),
		sec.axis = sec_axis(~ . - table_1$dilution[2], name = bquote(log[e]('eDNA Conc') ~ "air samples - PTFE filters")))+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"),  # Custom breaks every 2 weeks
		labels = date_format("%b %d")
	) +		
	theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title.x = element_blank(),
			legend.text = element_text(size = 14),  
			legend.title = element_text(size = 15, face = "bold"),
			
			# Left axis
			axis.ticks.y.left = element_line(color = "#61BEA4", size = 1.2),
			axis.title.y.left = element_text(color = "#61BEA4", size = 14),
			axis.text.y.left = element_text(color = "#61BEA4", size = 13),
			
			# Right axis
			axis.ticks.y.right = element_line(color = "#F49D4D", size = 1.2),
			axis.title.y.right = element_text(color = "#F49D4D", size = 14),
			axis.text.y.right = element_text(color = "#F49D4D", size = 13))

p4 <-
	ggplot()+
	geom_smooth(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',se=F,span=0.5)+
	geom_point(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',size=3)+
	geom_point(data=log_A %>% filter(log_A_FT=='MCE Air'),
							aes(x=as.Date(time),y=table_1$dilution[3]+log_A_mean),col='#D79FA7',size=6,pch=3)+
	geom_smooth(data=mean_log_A %>% filter(mean_log_A_FT=='MCE Air'),
							aes(x=as.Date(time),y=table_1$dilution[3]+mean_log_A_mean),col='#D79FA7',span=0.6)+
	scale_y_continuous(
		limits = c(8, 14),
		bquote(log[e]('eDNA Conc') ~ " water samples"),
		sec.axis = sec_axis(~ . - table_1$dilution[3], name = bquote(log[e]('eDNA Conc') ~ "air samples - MCE filters")))+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"),  # Custom breaks every 2 weeks
		labels = date_format("%b %d")
	) +		
	theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title.x = element_blank(),
			legend.text = element_text(size = 14),  
			legend.title = element_text(size = 15, face = "bold"),
			
			# Left axis
			axis.ticks.y.left = element_line(color = "#61BEA4", size = 1.2),
			axis.title.y.left = element_text(color = "#61BEA4", size = 14),
			axis.text.y.left = element_text(color = "#61BEA4", size = 13),
			
			# Right axis
			axis.ticks.y.right = element_line(color = "#D79FA7", size = 1.2),
			axis.title.y.right = element_text(color = "#D79FA7", size = 14),
			axis.text.y.right = element_text(color = "#D79FA7", size = 13)
		)

p5 <-
	ggplot()+
	geom_smooth(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',se=F,span=0.5)+
	geom_point(data=log_W,aes(x=as.Date(time),y=(log_W_mean)),col='#61BEA4',size=3)+
	geom_point(data=log_A %>% filter(log_A_FT=='MCE DI water'),
							aes(x=as.Date(time),y=table_1$dilution[4]+log_A_mean),col='#D85A44',size=6,pch=3)+
	geom_smooth(data=mean_log_A %>% filter(mean_log_A_FT=='MCE DI water'),
							aes(x=as.Date(time),y=table_1$dilution[4]+mean_log_A_mean),col='#D85A44',span=0.6)+
	scale_y_continuous(
		limits = c(8, 14),
		bquote(log[e]('eDNA Conc') ~ " water samples"),
		sec.axis = sec_axis(~ . - table_1$dilution[4], name = bquote(log[e]('eDNA Conc') ~ "air samples - DI water")))+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"),  # Custom breaks every 2 weeks
		labels = date_format("%b %d")
	) +		
	theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title.x = element_blank(),
			legend.text = element_text(size = 14),  
			legend.title = element_text(size = 15, face = "bold"),
			
			# Left axis
			axis.ticks.y.left = element_line(color = "#61BEA4", size = 1.2),
			axis.title.y.left = element_text(color = "#61BEA4", size = 14),
			axis.text.y.left = element_text(color = "#61BEA4", size = 13),
			
			# Right axis
			axis.ticks.y.right = element_line(color = "#D85A44", size = 1.2),
			axis.title.y.right = element_text(color = "#D85A44", size = 14),
			axis.text.y.right = element_text(color = "#D85A44", size = 13)
		)


p_leg <- cowplot::get_legend(
	data.frame(
		x = 1, 
		y = c(1, 2, 3, 4, 5, 6),
		group = c("1.Visual counts", "2.Water eDNA", "3.Air eDNA (Gelatin)", "4.Air eDNA (PTFE)", "5.Air eDNA (MCE)", "6.Air eDNA (water bucket)")) %>% 
		ggplot( aes(x = x, y = y, color = group)) +
		geom_line(size = 2) +
		scale_color_manual(values = c('black','#61BEA4', '#AB5971', '#F49D4D','#D79FA7','#D85A44'), name = "Observation method") +
		theme_minimal() +
		theme(legend.title = element_text(face = "bold",size = 16),
					legend.text = element_text(size = 14)))

pp1 <- cowplot::plot_grid(p2,p3,p4,p5,ncol = 2,align = 'v')
pp2 <- cowplot::plot_grid(p1,p_leg,ncol = 1,rel_heights = c(3,2))
fig_1 <- cowplot::plot_grid(pp2,pp1,ncol = 2,rel_widths = c(4,7))
fig_1
ggsave(fig_1,filename = here('Plots','Figure_1_new_optional.jpg'),width = 20,height = 10,dpi =300)

