ecosystem_engineers <- "~/Documents/R_Projects/Ecosystem_Engineers/";
setwd(ecosystem_engineers);
library(openxlsx);	#install.packages(openxlsx);
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');  #
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');  #
load("~/Documents/R_Projects/Data_for_R/Gradstein_2020_Augmented.RData"); # refined Gradstein 2020 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/PBDB_Data_Smøl.RData"); # refined Gradstein 2020 timescale & biozonations
land_plants <- c("Angiospermae","Coniferophyta","Ginkgophyta","Haptophyta","Pinophyta","Psilophytophyta","Pteridophyta","Spermatophyta","Tracheophyta");
#source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');  #

count_taxa_from_pbdb_finds <- function(taxon,pbdb_find_taxa)	{
return(sum(pbdb_find_taxa %in% taxon));
}

stage_slices <- gradstein_2020_emended$time_scale[gradstein_2020_emended$time_scale$scale %in% "Stage Slice",];
stage_slices$interval <- stage_slices$st;
stage_slices <- stage_slices[order(-stage_slices$ma_lb),];
pbdb_sites <- pbdb_data_list$pbdb_sites_refined;
pbdb_sites$bin_lb <- match(pbdb_sites$interval_lb,stage_slices$interval);
pbdb_sites$bin_ub <- match(pbdb_sites$interval_ub,stage_slices$interval);
pbdb_reef_sites <- pbdb_sites[pbdb_sites$environment %in% reef_environments,];
reef_phyla <- sort(unique(pbdb_data_list$pbdb_finds$phylum[pbdb_data_list$pbdb_finds$collection_no %in% pbdb_reef_sites$collection_no]));
reef_phyla <- reef_phyla[!reef_phyla %in% c("","NO_PHYLUM_SPECIFIED","Problematica",land_plants)];
rphyla <- length(reef_phyla);
relv_finds <- pbdb_data_list$pbdb_finds[!pbdb_data_list$pbdb_finds$phylum %in% land_plants,]
reef_classes <- sort(unique(relv_finds$class));
reef_classes <- reef_classes[!reef_classes %in% c("","NO_CLASS_SPECIFIED")];
rclass <- length(reef_classes);
nslice <- nrow(stage_slices);

archaeocyath_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$class=="Archaeocyatha",];
archaeocyath_sites <- pbdb_sites[pbdb_sites$collection_no %in% archaeocyath_finds$collection_no,];
max(archaeocyath_sites$ma_lb)
max(archaeocyath_sites$ma_ub)
min(archaeocyath_sites$ma_lb)
min(archaeocyath_sites$ma_ub)

stromatoporoid_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$class=="Stromatoporoidea",]
stromatoporoid_sites <- pbdb_reef_sites[pbdb_reef_sites$collection_no %in% stromatoporoid_finds$collection_no,];
max(stromatoporoid_sites$ma_lb)
max(stromatoporoid_sites$ma_ub)
min(stromatoporoid_sites$ma_lb)
min(stromatoporoid_sites$ma_ub)

rudists <- c("Caprinulidae","Caprotinidae","Diceratidae","Hippuritidae","Monopleuridae","Plagioptychidae","Polyconitidae","Radiolitidae","Trechmannellidae");
rudist_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$family %in% rudists,]
rudist_sites <- pbdb_sites[pbdb_sites$collection_no %in% rudist_finds$collection_no,];
max(rudist_sites$ma_lb)
max(rudist_sites$ma_ub)

bin_rocks <- bin_reef_rocks <- vector(length=nslice);
names(bin_rocks) <- stage_slices$st;
phylum_reef_counts <- data.frame(array(0,dim=c(nslice,rphyla)));
class_reef_counts <- data.frame(array(0,dim=c(nslice,rclass)));

rownames(class_reef_counts) <- rownames(phylum_reef_counts) <- stage_slices$st;
colnames(phylum_reef_counts) <- reef_phyla;
colnames(class_reef_counts) <- reef_classes;
phylum_reef_dominance <- phylum_reef_counts;
class_reef_dominance <- class_reef_counts;
for (ns in 1:nslice)	{
	slice_sites <- pbdb_sites[pbdb_sites$bin_lb==ns & pbdb_sites$bin_ub==ns,];
	slice_sites_reefs <- slice_sites[slice_sites$environment %in% reef_environments,];
	slice_sites_reefs <- slice_sites_reefs[!slice_sites_reefs$environment %in% "perireef or subreef",];
	if (nrow(slice_sites_reefs)>0)	{
		bin_reef_rocks[ns] <- length(bin_rock_nos);
		slice_reef_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$collection_no %in% slice_sites_reefs$collection_no,];
		bin_rock_nos <- unique(slice_sites_reefs$pbdb_rock_no_sr[slice_sites_reefs$pbdb_rock_no_sr>0]);
		taxon <- reef_phyla;
		phylum_reef_counts[ns,] <- sapply(taxon,count_taxa_from_pbdb_finds,pbdb_find_taxa=slice_reef_finds$phylum);
		taxon <- reef_classes;
		class_reef_counts[ns,] <- sapply(taxon,count_taxa_from_pbdb_finds,pbdb_find_taxa=slice_reef_finds$class);
		rr <- 0;
		while (rr < bin_reef_rocks[ns])	{
			rr <- rr+1;
			rock_coll_nos <- slice_sites_reefs$collection_no[slice_sites_reefs$pbdb_rock_no_sr == bin_rock_nos[rr]]
			this_reef_finds <- slice_reef_finds[slice_reef_finds$collection_no %in% rock_coll_nos,];
#			rock_coll_nos
			taxon <- reef_phyla;
			phylum_counts <- sapply(taxon,count_taxa_from_pbdb_finds,pbdb_find_taxa=this_reef_finds$phylum);
			if (max(phylum_counts)>1)	{
				dom_phylum <- (1:rphyla)[phylum_counts %in% max(phylum_counts)];
				pc <- bin_reef_rocks[ns]/sum(phylum_counts %in% max(phylum_counts));
				phylum_reef_dominance[ns,dom_phylum] <- phylum_reef_dominance[ns,dom_phylum]+pc;
				}
			taxon <- reef_classes;
			class_counts <- sapply(taxon,count_taxa_from_pbdb_finds,pbdb_find_taxa=this_reef_finds$class);
			if (max(class_counts)>1)	{
				dom_class <- (1:rclass)[class_counts %in% max(class_counts)];
				cc <- bin_reef_rocks[ns]/sum(class_counts %in% max(class_counts));
				class_reef_dominance[ns,dom_class] <- class_reef_dominance[ns,dom_class]+cc;
				}
			}
		}
	}
write.csv(phylum_reef_counts,"Reef_Phyla.csv",row.names = T);
write.csv(class_reef_counts,"Reef_Classes.csv",row.names = T);
write.csv(phylum_reef_dominance,"Reef_Phylum_Dominance.csv",row.names = T);
write.csv(class_reef_dominance,"Reef_Class_Dominance.csv",row.names = T);
reef_props <- bin_reef_rocks/bin_rocks;

strat_unit <- phanerozoic;
# Get those Strat Units from the Gradstein scale file
standard_time_scale <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$interval %in% strat_unit);
# Make sure that ages on time scale are negatives
standard_time_scale$ma_lb <- -abs(standard_time_scale$ma_lb);
standard_time_scale$ma_ub <- -abs(standard_time_scale$ma_ub);
# Make sure that time scale is correct order, with oldest first
standard_time_scale <- standard_time_scale[order(standard_time_scale$ma_lb),];
# Add names/symbols for plotting
strat_names <- standard_time_scale$st;
# Get the colors for time units
strat_colors <- standard_time_scale$color;
# Set the oldest and youngest intervals by names on whatever chronostratigraphic scale you used in the analysis
oldest_interval <- "Ediacaran";
youngest_interval <- "Quaternary";
# get the time scale that will be plotted
time_scale_to_plot <- unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub));
# Now, set the onset and end of the x-axis
onset <- min(time_scale_to_plot);
end <- max(time_scale_to_plot);
# Finally, set up breaks for x-axis
yearbreaks <- c(5,25,50);					# set breaks for x-axis (minor, medium & major)

# now, set up the y-axis: this will reflect your data
mxy <- 0.60;
mny <- 0;									# set maximum y value

# set up plot details:
use_strat_labels <- T;						# if T, then strat_names will be plotted on X-axis inside boxes
strat_label_size <- 10;					# size of the labels for chronostratigraphic units
alt_back <- F;								# if T, then the background will alternat shades between major intervals
ordinate <- "   Prop. PBDB Marine Rock Units";								# Label of Y-axis
hues <- "T";								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
ysize <- 4*4.285714285/6;
xsize <- 4; 

strat_names[strat_names=="Q"] <- "";
strat_names[strat_names=="Ng"] <- "N";
plot_title <- "Rocks with Sponge-Heavy Reefs*"			# Name of the plot; enter "" for nothing
Phanerozoic_Timescale_Plot_Flexible(onset,end,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels=T,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=15);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.1,med_break=0.05,min_break=0.01,linewd=4/3,orient=2,print_label=TRUE);

for (ns in 1:nslice)	if (reef_props[ns]>0)	rect(-stage_slices$ma_lb[ns],0,-stage_slices$ma_ub[ns],min(mxy,reef_props[ns]),col=stage_slices$color[ns],lwd=0.25);

plot_title <- "Rocks with Coral-Heavy Reefs*"			# Name of the plot; enter "" for nothing
Phanerozoic_Timescale_Plot_Flexible(onset,end,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels=T,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=15);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.1,med_break=0.05,min_break=0.01,linewd=4/3,orient=2,print_label=TRUE);


{}
