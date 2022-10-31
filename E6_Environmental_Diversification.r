# Setup Program (last modified 2022-04-28) ####
common_source_folder <- file.choose();	# grab Chronos.r;
print("Grab 'Chronos.r' from your common source folder)");
common_source_folder <- gsub("Chronos.r","",common_source_folder);
source(paste(common_source_folder,"Chronos.r",sep=""));
source(paste(common_source_folder,"Data_Downloading_v4.r",sep=""));
source(paste(common_source_folder,"Disparity.r",sep=""));
source(paste(common_source_folder,"General_Plot_Templates.r",sep=""));
source(paste(common_source_folder,"Historical_Diversity_Metrics.r",sep=""));
source(paste(common_source_folder,"Nexus_File_Routines.r",sep=""));
source(paste(common_source_folder,"Sampling_and_Occupancy_Distributions.r",sep=""));
source(paste(common_source_folder,"Wagner_Kluges.r",sep=""));
source(paste(common_source_folder,"Wagner_Stats_and_Probability_101.r",sep=""));

print("Grab 'PBDB_Data_Smøl.RData' from your data folder)");
data_for_R_folder <- file.choose();	# Grab PBDB_Data_Smøl.RData
data_for_R_folder <- gsub("PBDB_Data_Smøl.RData","",data_for_R_folder);
load(paste(data_for_R_folder,"PBDB_Data_Smøl.RData",sep="")); # load PBDB data
load(paste(data_for_R_folder,"Gradstein_2020_Augmented.RData",sep="")); # refined Gradstein 2020 timescale & biozonations

# Load Data Sets ####
time_scale <- gradstein_2020_emended$time_scale;
int_stages <- time_scale[time_scale$chronostratigraphic_rank=="Stage" & time_scale$scale=="International" & time_scale$interval_sr=="",];
int_stages <- int_stages[!is.na(int_stages$interval),];
temporal_precision <- 0.1;
pbdb_taxonomy <- pbdb_data_list_for_class$pbdb_taxonomy;
pbdb_finds <- pbdb_data_list_for_class$pbdb_finds;
pbdb_environments <- pbdb_data_list_for_class$pbdb_environments;
nfinds <- nrow(pbdb_finds);
pbdb_finds$test_group <- rep("other",nfinds);
pbdb_finds <- occidere_indeterminate_species(pbdb_finds);
pbdb_finds <- remove_duplicate_occurrences_paleodb(pbdb_finds);
pbdb_finds <- pbdb_finds[pbdb_finds$accepted_rank %in% taxonomic_rank[1:4],];
pbdb_finds <- pbdb_finds[!pbdb_finds$genus %in% "NO_GENUS_SPECIFIED",];
pbdb_sites <- pbdb_data_list_for_class$pbdb_sites;
finest_chronostrat <- pbdb_data_list_for_class$time_scale;
nslices <- nrow(finest_chronostrat);
age <- pbdb_sites$ma_lb;
pbdb_sites$bin_lb <- match(pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat),finest_chronostrat$interval);
age <- pbdb_sites$ma_ub;
pbdb_sites$bin_ub <- match(pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat),finest_chronostrat$interval);

# Choose your adventure! (THIS IS WHERE YOU TELL THE PROGRAM WHAT YOU WANT TO DO) ####
# Bivalve environmental breadth
test_taxa <- c("Bivalvia");
contrast_intervals <- c("Phanerozoic");
unique(pbdb_taxonomy$composition)
sampling_control_group <- c("Gastropoda","Hydrozoa","Hyolithomorpha","Orthothecimorpha","Nautiloidea","Polyplacophora","Rostroconchia","Scaphopoda","Tergomya","Helcionelloida");

# Let's start organizing the data given your request ####
test_taxa <- sort(test_taxa);
nttaxa <- length(test_taxa);
ma_start <- max(time_scale$ma_lb[match(contrast_intervals,time_scale$interval)]);
ma_end <- min(time_scale$ma_ub[match(contrast_intervals,time_scale$interval)]);

# first, reduce PBDB sites & finds to those within the time range.
pbdb_sites <- pbdb_sites[pbdb_sites$ma_lb>ma_end & pbdb_sites$ma_ub<ma_start,];
pbdb_finds <- pbdb_finds[pbdb_finds$collection_no %in% pbdb_sites$collection_no,];
nfinds <- nrow(pbdb_finds);
nsites <- nrow(pbdb_sites);

#### Identify the sampling control finds ####
control_ranks <- pbdb_taxonomy$accepted_rank[match(sampling_control_group,pbdb_taxonomy$accepted_name)];
for (ct in 1:length(control_ranks))	{
	if (control_ranks[ct] %in% colnames(pbdb_finds))	{
		tr <- match(control_ranks[ct],colnames(pbdb_finds));
		pbdb_finds$test_group[pbdb_finds[,tr]==sampling_control_group[ct]] <- "control";
		} else	{
		daughters <- accersi_daughter_taxa(parent_taxon=sampling_control_group[ct],pbdb_taxonomy = pbdb_taxonomy);
		daughters$daughter_rank[daughters$daughter_rank=="subgenus"] <- "genus";
		for (dd in 1:nrow(daughters))	{
			tr <- match(daughters$daughter_rank[dd],colnames(pbdb_finds))
			pbdb_finds$test_group[pbdb_finds[,tr] %in% daughters$daughter[dd]] <- "control";
			}
		}
	}
cfinds <- sum(pbdb_finds$test_group=="control");

#### Identify the test group finds ####
taxon_ranks <- pbdb_taxonomy$accepted_rank[match(test_taxa,pbdb_taxonomy$taxon_name)];
for (tx in 1:nttaxa)	{
	if (taxon_ranks[tx] %in% colnames(pbdb_finds))	{
		tr <- match(taxon_ranks[tx],colnames(pbdb_finds));
		pbdb_finds$test_group[pbdb_finds[,tr]==test_taxa[tx]] <- test_taxa[tx];
		} else	{
		daughters <- accersi_daughter_taxa(parent_taxon=test_taxa[tx],pbdb_taxonomy = pbdb_taxonomy,returnable_ranks=c("genus","sugenus"));
		daughters$daughter_rank[daughters$daughter_rank=="subgenus"] <- "genus";
		for (dd in 1:nrow(daughters))	{
			tr <- match(daughters$daughter_rank[dd],colnames(pbdb_finds))
			pbdb_finds$test_group[pbdb_finds[,tr] %in% daughters$daughter[dd]] <- test_taxa[tx];
			}
		}
#		print(unique(pbdb_finds$test_group))
	}

# Put together different counts of ecological distribution breadth ####
good_bins <- sort(unique(pbdb_sites$bin_lb[pbdb_sites$bin_lb==pbdb_sites$bin_ub]));
good_bins <- good_bins[good_bins>=min(pbdb_sites$bin_lb[match(pbdb_finds$collection_no[pbdb_finds$test_group %in% test_taxa],pbdb_sites$collection_no)])];

g_b <- length(good_bins);
test_group_breadth <- data.frame(interval=as.character(),
								 test_taxon_rocks=as.numeric(),cntrl_taxa_rocks=as.numeric(),all_rocks=as.numeric(),
								 test_taxon_sites=as.numeric(),cntrl_taxa_sites=as.numeric(),all_sites=as.numeric(),
								 test_taxon_finds=as.numeric(),cntrl_taxa_finds=as.numeric(),all_finds=as.numeric());

for (gb in 1:g_b)	{
	bn <- good_bins[gb];
	bin_sites <- pbdb_sites[(pbdb_sites$bin_lb==bn & pbdb_sites$bin_ub==bn),];
	bin_rocks <- unique(bin_sites$pbdb_formation_no[bin_sites$pbdb_formation_no>0]);
	bin_finds <- pbdb_finds[pbdb_finds$collection_no %in% bin_sites$collection_no,];
	cntrl_finds <- bin_finds[bin_finds$test_group %in% c(test_taxa,"control"),];
	cntrl_sites <- bin_sites[bin_sites$collection_no %in% cntrl_finds$collection_no,];
	cntrl_rocks <- unique(cntrl_sites$pbdb_formation_no[cntrl_sites$pbdb_formation_no>0]);
#	unique(cntrl_sites$environment)
	if (length(cntrl_rocks)>10)	{
		test_finds <- bin_finds[bin_finds$test_group %in% test_taxa,];
		test_sites <- bin_sites[bin_sites$collection_no %in% test_finds$collection_no,];
		test_rocks <- unique(test_sites$pbdb_formation_no[test_sites$pbdb_formation_no>0]);
		relv_bin_finds <- bin_finds[bin_finds$collection_no %in% cntrl_sites$collection_no,];

		dummy_frame <- data.frame(interval=as.character(finest_chronostrat$interval[bn]),
								  test_taxon_rocks=length(test_rocks),
								  cntrl_taxa_rocks=length(cntrl_rocks),
								  all_rocks=length(bin_rocks),
								  test_taxon_sites=nrow(test_sites),
								  cntrl_taxa_sites=nrow(cntrl_sites),
								  all_sites=nrow(bin_sites),
								  test_taxon_finds=sum(relv_bin_finds$test_group %in% test_taxa),
								  cntrl_taxa_finds=nrow(relv_bin_finds),
								  all_finds=nrow(bin_finds)
								  );
		test_group_breadth <- rbind(test_group_breadth,dummy_frame);
		}
	}
test_group_breadth$prop_cntrl_rocks <- test_group_breadth$test_taxon_rocks/test_group_breadth$cntrl_taxa_rocks;
test_group_breadth$prop_all_rocks <- test_group_breadth$test_taxon_rocks/test_group_breadth$all_rocks;
test_group_breadth$prop_cntrl_sites <- test_group_breadth$test_taxon_sites/test_group_breadth$cntrl_taxa_sites;
test_group_breadth$prop_all_sites <- test_group_breadth$test_taxon_sites/test_group_breadth$all_sites;
test_group_breadth$prop_cntrl_sites <- test_group_breadth$test_taxon_sites/test_group_breadth$cntrl_taxa_sites;
test_group_breadth$prop_all_sites <- test_group_breadth$test_taxon_sites/test_group_breadth$all_sites;
test_group_breadth$prop_cntrl_finds <- test_group_breadth$test_taxon_finds/test_group_breadth$cntrl_taxa_finds;
test_group_breadth$prop_all_finds <- test_group_breadth$test_taxon_finds/test_group_breadth$all_finds;

# Now, let's plot things out over time ####
par(mfrow=c(1,1));
ngroups <- length(test_taxa);
standard_time_scale <- accersi_appropriate_stratigraphic_intervals_for_plots(start=max(finest_chronostrat$ma_lb[good_bins]),finish=min(finest_chronostrat$ma_ub[good_bins]),time_scale=gradstein_2020_emended$time_scale);
time_scale_to_plot <- c(standard_time_scale$ma_lb,standard_time_scale$ma_ub[nrow(standard_time_scale)]);
yearbreaks <- sort(as.numeric(set_axis_breaks(max(abs(time_scale_to_plot)),min(abs(time_scale_to_plot)))));
strat_names <- standard_time_scale$st;
strat_names[strat_names=="Q"] <- "";
# set up plot details:
use_strat_labels <- T;						# if T, then strat_names will be plotted on X-axis inside boxes
alt_back <- F;								# if T, then the background will alternat shades between major intervals
hues <- "T";								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
ysize <- 4*4.285714285/6;
xsize <- 4;
strat_label_size <- 12.5;					# size of the labels for chronostratigraphic units

mxy <- ceiling(max(test_group_breadth$prop_cntrl_rocks)/0.1)/10;
ordinate <- "Prop. of Rocks w. Aragonitic Fossils"
Phanerozoic_Timescale_Plot_Flexible(onset=min(time_scale_to_plot),end=max(time_scale_to_plot),time_scale_to_plot,mxy=mxy,mny=-(mxy)/20,use_strat_labels,strat_names=strat_names,strat_colors=standard_time_scale$color,plot_title="Bivalve Formation Occupancy",ordinate=ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size);
specified_axis(axe=2,max_val=mxy,min_val=0,maj_break=0.1,med_break=0.05,min_break=0.01,linewd=4/3,orient=2,print_label=TRUE);

for (bn in 1:nrow(test_group_breadth))	{
	bnn <- match(test_group_breadth$interval[bn],finest_chronostrat$interval);
	myr <- -abs((finest_chronostrat$ma_lb[bnn]+finest_chronostrat$ma_ub[bnn])/2);
	bnc <- finest_chronostrat$color[bnn];
	points(myr,test_group_breadth$prop_cntrl_rocks[bn],pch=21,bg=bnc,cex=0.75);
	}

mxy <- ceiling(max(test_group_breadth$prop_cntrl_sites)/0.1)/10;
ordinate <- "Prop. of Sites w. Aragonitic Fossils"
Phanerozoic_Timescale_Plot_Flexible(onset=min(time_scale_to_plot),end=max(time_scale_to_plot),time_scale_to_plot,mxy=mxy,mny=-(mxy)/20,use_strat_labels,strat_names=strat_names,strat_colors=standard_time_scale$color,plot_title="Bivalve Site Occupancy",ordinate=ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size);
specified_axis(axe=2,max_val=mxy,min_val=0,maj_break=0.1,med_break=0.05,min_break=0.01,linewd=4/3,orient=2,print_label=TRUE);

for (bn in 1:nrow(test_group_breadth))	{
	bnn <- match(test_group_breadth$interval[bn],finest_chronostrat$interval);
	myr <- -abs((finest_chronostrat$ma_lb[bnn]+finest_chronostrat$ma_ub[bnn])/2);
	bnc <- finest_chronostrat$color[bnn];
	points(myr,test_group_breadth$prop_cntrl_sites[bn],pch=21,bg=bnc,cex=0.75);
	}

mxy <- ceiling(max(test_group_breadth$prop_all_finds)/0.1)/10;
ordinate <- "Prop. of Occurrences from Sites w. Aragonitic Fossils"
Phanerozoic_Timescale_Plot_Flexible(onset=min(time_scale_to_plot),end=max(time_scale_to_plot),time_scale_to_plot,mxy=mxy,mny=-(mxy)/20,use_strat_labels,strat_names=strat_names,strat_colors=standard_time_scale$color,plot_title="Bivalve Find Proportions",ordinate=ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size);
specified_axis(axe=2,max_val=mxy,min_val=0,maj_break=0.1,med_break=0.05,min_break=0.01,linewd=4/3,orient=2,print_label=TRUE);

for (bn in 1:nrow(test_group_breadth))	{
	bnn <- match(test_group_breadth$interval[bn],finest_chronostrat$interval);
	myr <- -abs((finest_chronostrat$ma_lb[bnn]+finest_chronostrat$ma_ub[bnn])/2);
	bnc <- finest_chronostrat$color[bnn];
	points(myr,test_group_breadth$prop_all_finds[bn],pch=21,bg=bnc,cex=0.75);
	}

# End ####
{}

list_names <- c(names(pbdb_data_list_for_class),"pbdb_environments");
pbdb_environments <- as.data.frame(readxl::read_xlsx(file.choose()));
pbdb_data_list_for_class <- rlist::list.append(pbdb_data_list_for_class,pbdb_environments);
names(pbdb_data_list_for_class) <- list_names;
load(paste(data_for_R_folder,"Paleobiology_Database.RData",sep=""));
list_names <- c(names(pbdb_data_list),"pbdb_environments");
pbdb_data_list <- rlist::list.append(pbdb_data_list,pbdb_environments);
names(pbdb_data_list) <- list_names;
save(pbdb_data_list,file=paste(data_for_R_folder,"Paleobiology_Database.RData",sep=""));
save(pbdb_data_list_for_class,file=paste(data_for_R_folder,"PBDB_Data_Smøl.RData",sep=""));

{}