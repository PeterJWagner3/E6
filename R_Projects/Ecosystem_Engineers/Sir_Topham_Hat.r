ecosystem_engineers <- "~/Documents/R_Projects/Ecosystem_Engineers/";
setwd(ecosystem_engineers);
library(openxlsx);	#install.packages(openxlsx);
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');  #
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');  #
load("~/Documents/R_Projects/Data_for_R/Gradstein_2020_Augmented.RData"); # refined Gradstein 2020 timescale & biozonations
#source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');  #

thomas <- as.data.frame(read.xlsx("Ecosystem_Earth_Systems_Engineers_Table_1.xlsx"));

strat_unit <- periods;
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
oldest_interval <- "EA";
youngest_interval <- "Q";
# get the time scale that will be plotted
time_scale_to_plot <- unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub));
# Now, set the onset and end of the x-axis
onset <- min(time_scale_to_plot);
end <- max(time_scale_to_plot);
# Finally, set up breaks for x-axis
yearbreaks <- c(50,250,500);					# set breaks for x-axis (minor, medium & major)

# now, set up the y-axis: this will reflect your data
mxy <- nrow(thomas)+1;	# set maximum y value
mny <- 0;									# set maximum y value

# set up plot details:
use_strat_labels <- F;						# if T, then strat_names will be plotted on X-axis inside boxes
strat_label_size <- 2.5;					# size of the labels for chronostratigraphic units
alt_back <- F;								# if T, then the background will alternat shades between major intervals
plot_title <- "Ecosystem Engineers over Time"			# Name of the plot; enter "" for nothing
ordinate <- "";								# Label of Y-axis
hues <- "T";								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
ysize <- 4*4.285714285/6;
xsize <- 4; 

Phanerozoic_Timescale_Plot_Flexible(onset,end,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size);

thomas$Ma_End_lb

{}