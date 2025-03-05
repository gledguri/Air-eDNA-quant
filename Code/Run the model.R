# Libraries -----------------------------------------------------------------------------------

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(QM);#devtools::install_github("gledguri/QM")
library(rstan);options(mc.cores = parallel::detectCores())
library(tibble)
library(stringr)
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

# scientific_10 <- function(x) {
# 	c <- scales::scientific_format()(x)
# 	t <- gsub("e", "%*%10^", c)
# 	t2 <- gsub("10\\^\\+", "10\\^", t)
# 	str2expression(t2)
# }

scientific_10 <- function(x) {
	c <- log10(abs(x))
	expo <- floor(c)
	base <- round(x / 10^expo,3)
	formatted <- paste0(base,'%*%',"10^", expo)
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
	mutate(Date=if_else(Date==20241025,20241024,Date)) %>% 
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

qpcr_coho_air_en <-
	qpcr_coho_en_prep %>% 
	# filter(FilterType!='Gelatin') %>%
	# filter(FilterType!='PTFE') %>%
	# filter(FilterType!='MCE Air') %>%
	arrange(FilterType,Date) %>% 
	make_index("FilterType", "f_idx") %>%
	make_index("Date", "t_idx") %>% 
	make_index("Sample", "j_idx") %>% 
	combine_index(c('f_idx','t_idx','j_idx'),'idx') %>%
	combine_index(c('f_idx','t_idx'),'ft_idx') %>% 
	# arrange(j_idx) %>% 
	as.data.frame()

count_fish_effort_selected <-
	count_fish_effort %>% filter(time%in%unique(qpcr_coho_en_prep$time)) %>% 
	mutate(time=as.character(time)) %>% 
	left_join(.,
						qpcr_coho_wat_en %>% distinct(time,t_idx),
						by=c('time'='time'))

samp_tab_ijb <- qpcr_coho_air_en %>% distinct(FilterType,time,f_idx,t_idx,ft_idx,idx)
samp_tab_ij <- samp_tab_ijb %>% distinct(FilterType,time,f_idx,t_idx,ft_idx)

bio_rep_idx =
	samp_tab_ijb %>% mutate(pres=1) %>% 
	group_by(FilterType,time,f_idx,t_idx,ft_idx) %>% 
	summarise(bio_rep_idx=sum(pres)) %>% ungroup() %>% 
	pull(bio_rep_idx)
N_bio_rep_param = sum(bio_rep_idx-1)
N_bio_rep_idx <- length(bio_rep_idx)
N_bio_rep_RE <- sum(bio_rep_idx)

# This is supposed to be TRUE!
identical(qpcr_coho_air_en %>% distinct(time,t_idx), qpcr_coho_wat_en %>% distinct(time,t_idx))

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
									# t_lambda_idx = c(1,1,2,3,4,4),
									l_i = count_fish_effort_selected$t_idx,
									a_j = samp_tab_ijb %>% pull(f_idx),
									a_i = samp_tab_ijb %>% pull(t_idx),
									a_ij = samp_tab_ijb %>% pull(ft_idx),
									a_ijb = samp_tab_ijb %>% pull(idx),
									bio_rep = samp_tab_ijb %>% mutate(bio_rep=if_else(FilterType=='Gelatin',1,0)) %>% pull(bio_rep),
									N_bio_rep_param = N_bio_rep_param,
									N_bio_rep_RE = N_bio_rep_RE,
									N_bio_rep_idx = N_bio_rep_idx,
									bio_rep_idx = bio_rep_idx,
									tau_bio_rep_idx = samp_tab_ij$f_idx,
									# Mm = model.matrix(~ factor(t_idx) - 1, data = samp_tab_ijb) %>% t(),
									# 
									time = 4,
									tau_p1 = -1,tau_p2=1,
									N = count_fish_effort_selected$Total,
									E = count_fish_effort_selected$time_days,
									logit_phi_mu = -3,
									logit_phi_sd = 1,
									label_time = qpcr_coho_air_en %>% distinct(time,t_idx),
									label_filter = qpcr_coho_air_en %>% distinct(FilterType,f_idx),
									label_wat = qpcr_coho_wat_en %>% distinct(idx, Sample) %>% arrange(idx) %>% as.data.frame(),
									label_air = qpcr_coho_air_en %>% distinct(idx, Sample) %>% arrange(idx) %>% as.data.frame()
)


stanMod_count <- stan(
	file = here('Code','Count_model.stan'),
	data = stan_data,
	iter = 10000,
	warmup = 5000,
	chains = 4)

# saveRDS(stanMod_count,here('Output','stanMod_output.rds'))
# saveRDS(stan_data,here('Output','stan_data_input.rds'))

post_table_raw <- 
	bind_cols(stan_data$a_ij %>% as.data.frame() %>% setNames('a_ij') %>% mutate(a_ij=as.numeric(a_ij)),
						stan_data$a_i %>% as.data.frame() %>% setNames('a_i'),
						stan_data$a_j %>% as.data.frame() %>% setNames('a_j'),
						stan_data$a_ijb %>% as.data.frame() %>% setNames('a_ijb')) %>% 
	left_join(.,
						join_ext_param(stanMod_count,'log_W') %>% 
							select(mean,g_idx) %>% 
							rename(log_W='mean'),
						by=c('a_i'='g_idx')) %>% 
	left_join(.,
						cbind(extract_param(stanMod_count,'lambda'),
									as.data.frame(stan_data$E) %>% setNames('E'),
									as.data.frame(count_fish_effort_selected$t_idx) %>% setNames('g_idx')) %>% 
							select(mean,E,g_idx) %>%
							rename(lambda='mean'),
						by=c('a_i'='g_idx')) %>%
	left_join(.,
						join_ext_param(stanMod_count,'log_A') %>% 
							select(mean,g_idx) %>% 
							rename(log_A='mean'),
						by=c('a_ijb'='g_idx')) %>% 
	left_join(.,
						join_ext_param(stanMod_count,'bio_rep_RE') %>% 
							select(mean,g_idx) %>% 
							rename(bio_rep_RE='mean'),
						by=c('a_ijb'='g_idx')) %>% 
	left_join(.,
						join_ext_param(stanMod_count,'X_STATE') %>%
							select(mean,g_idx) %>%
							rename(X_STATE='mean'),
						by=c('a_i'='g_idx')) %>%
	left_join(.,
						join_ext_param(stanMod_count,'eta') %>% 
							select(mean,g_idx) %>% 
							rename(eta='mean'),
						by=c('a_j'='g_idx')) %>% 
	left_join(.,
						join_ext_param(stanMod_count,'epsilon') %>% 
							select(mean,g_idx) %>% 
							rename(epsilon='mean'),
						by=c('a_ij'='g_idx')) %>% 
	mutate(omega=
				 	extract_param(stanMod_count,'omega') %>% 
				 	pull(mean)) %>% 
	left_join(.,
						bind_cols(
							join_ext_param(stanMod_count,'mu_en_air') %>% select(mean) %>% 
								setNames(c('mu_A')),
							join_ext_param(stanMod_count,'sigma_en_air') %>% select(mean,g_idx) %>% 
								setNames(c('sigma_A','g_idx')),
							stan_data$j_qen_air_p_idx %>% as.data.frame() %>% setNames('j_qen_air_p_idx')
						) %>% 
							group_by(j_qen_air_p_idx) %>% 
							summarise(mu_A=mean(mu_A),
												sigma_A=mean(sigma_A)) %>% as.data.frame(),
						by=c('a_ijb'='j_qen_air_p_idx')) %>% 
	left_join(.,
						bind_cols(
							join_ext_param(stanMod_count,'psi_un_air') %>% select(mean),
							stan_data$j_qen_air_idx %>% as.data.frame() %>% setNames('j_qen_air_idx')) %>% 
							group_by(j_qen_air_idx) %>% 
							summarise(psi_un_air=mean(mean)) %>% as.data.frame(),
						by=c('a_ijb'='j_qen_air_idx'))

post_table <- post_table_raw %>%
	left_join(stan_data$label_time,
						by=c('a_i'='t_idx')) %>% 
	left_join(stan_data$label_filter,
						by=c('a_j'='f_idx'))

table_1 <-
	post_table %>% 
	group_by(FilterType) %>% 
	summarise(dilution=first(eta),
						SEE=mean(abs(epsilon)),
						bio_rep_delta=mean(abs(bio_rep_RE)))


omega <- post_table$omega %>% unique()

library(ggtext)
library(scales)
library(ggnewscale)

p1 <-
	post_table %>% filter(!duplicated(a_i)) %>% 
	ggplot()+
	geom_smooth(aes(x=as.Date(time),y=(X_STATE)),color='black',span=0.5)+
	geom_point(aes(x=as.Date(time),y=(lambda/E)),col='#AB5971',size=5)+
	labs(
		y = "<span style='color:black;'>X - Fish density (count &times; day<sup>-1</sup>)</span><br><br><span style='color:#AB5971;'> (λ &times; E<sup>-1</sup>)</span>"
	) +
	scale_y_log10()+
	scale_x_date(breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week")) +	
	theme_bw()+
	theme(
		strip.text = element_blank(),
		axis.title.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.y = element_markdown( size = 16),
		# axis.title.y = element_text(color = "black", size = 16),
		axis.text.y = element_text(color = "black", size = 14),
		legend.position = 'none',
	)


p2 <-
	post_table %>% filter(!duplicated(a_i)) %>%
	ggplot()+
	geom_smooth(aes(x=as.Date(time),y=(X_STATE/exp(omega))),color='black',span=0.5)+
	geom_point(aes(x=as.Date(time),y=exp(log_W)),col='deepskyblue2',size=5)+
	labs(y=bquote('eDNA concentration (water)'))+
	scale_y_log10(labels=scientific_10,breaks=c(50000,100000,200000,400000))+
	theme_bw()+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"), 
		labels = date_format("%b %d")) +	
	theme(strip.text = element_blank(),
				axis.title.x = element_blank(),
				axis.text.x = element_text(size=14),
				
				axis.title.y = element_text(color = "black", size = 16),
				axis.text.y = element_text(color = "black", size = 14),
				
				legend.position = 'none')

p3 <-
	post_table %>% 
	mutate(log_A=if_else(psi_un_air < -1.3,(-1/0),log_A)) %>% 
	ggplot()+
	geom_smooth(aes(x=as.Date(time),y=exp(log(X_STATE)-omega+eta)),color='black',span=0.5)+
	geom_point(aes(x=as.Date(time),y=exp(log_A),col=FilterType),pch=19, size=5)+
	scale_y_log10(labels=scientific_10,breaks=c(5,10,20,40,80,200,500,1000,2000))+
	labs(y=bquote('eDNA concentration (air)'))+
	scale_x_date(
		breaks = seq(as.Date("2024-10-17"), as.Date("2024-11-21"), by = "1 week"), 
		labels = date_format("%b %d")) +	
	scale_color_manual(values=c('#61BEA4','#D79FA7','#F49D4D','#D85A44'))+	
	geom_abline(intercept = -1.5,slope=0,lty=2)+
	facet_wrap(~FilterType,scales='free_y')+
	theme_bw()+
	theme(axis.title.x = element_blank(),
				strip.text = element_blank(),
				axis.title.y = element_text(size=16),
				axis.text.x = element_text(size=14),
				axis.text.y = element_text(size=14),
				legend.text = element_text(size = 14),
				legend.title = element_text(size = 15),
				legend.position = 'none')

leg_1 <- data.frame(x = 1:2,y1 = rep(10,2),y2 = rep(10,2),
										group = factor(c(1:2))) %>% 
	ggplot() +
	geom_line(aes(x = x, y = y1, linetype = 'X (fish/day)'), color = 'black', size = 1.2) +
	scale_linetype_manual(
		name = 'Unobserved state',
		values = 'solid') +
	guides(linetype = guide_legend(order = 1)) +
	
	geom_point(data = . %>% filter(group == 1), aes(x = x, y = y1, color = group), size = 4) +
	scale_color_manual(
		name = 'Observed state - visually',
		values = c('#AB5971'),
		labels = c('Visual counts')) +
	guides(color = guide_legend(order = 2)) +
	
	new_scale_color() +
	
	geom_point(data = . %>% filter(group == 2), aes(x = x, y = y1, color = group), size = 4) +
	scale_color_manual(
		name = 'Observed state - water eDNA',
		values = c('deepskyblue2'),
		labels = c('MCE (water samples)')) +
	guides(color = guide_legend(order = 3)) +
	
	theme_minimal() +
	theme(
		legend.position = 'right',
		legend.spacing.y = unit(0.1, "cm"),
		legend.key.size = unit(0.4, "cm"),
		legend.title = element_text(size = 14),
		legend.text = element_text(size = 14))

leg_2 <- data.frame(x = 1:4,y1 = rep(10,4),y2 = rep(10,4),
										group = factor(c(1:4))) %>% 
	ggplot() +
	geom_point(aes(x = x, y = y2, color = group), size = 4) +
	scale_color_manual(
		name = 'Observed state - air eDNA',
		values=c('#61BEA4','#D79FA7','#F49D4D','#D85A44'),
		labels = c('Gelatin', 'MCE (air eDNA)', 'MCE (Mili-q water)', 'PTFE')) +
	
	theme_minimal() +
	theme(
		legend.position = 'right',
		# legend.spacing.y = unit(0.4, "cm"),
		legend.key.size = unit(0.8, "cm"),
		legend.title = element_text(size = 14),
		legend.text = element_text(size = 14))

legend_1 <- cowplot::get_legend(leg_1+theme(legend.justification = c(0.5, 0.2)))
legend_2 <- cowplot::get_legend(leg_2+theme(legend.justification = c(-0.3, 0.2)))
legend <- cowplot::plot_grid(legend_1,legend_2)
legend

fig_1 <-
	cowplot::plot_grid(
		cowplot::plot_grid(
			cowplot::plot_grid(p1,p2,labels = c('A','B'),ncol = 1,align = 'v'),
			p3,rel_widths = c(1.1,2),labels=c('','C')),
		legend,ncol = 1,rel_heights = c(5,1.2))

ggsave(here('Plots','Figure_1.jpg'),fig_1,height = 10,width = 17,dpi =300)

source(here('Code','Dianostic_plots.R'))
