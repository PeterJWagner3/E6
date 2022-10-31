# setup functions ####
update_opinions <- function(taxon,ops_modified_after,inc_children=T)	{
#	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod")
options(warn=-1);
if (inc_children)	{
	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}

pbdb_data <- read.csv(httpTO,header=T,stringsAsFactors = F);
new_opinions <- put_pbdb_dataframes_into_proper_type(pbdb_data);
options(warn=1);
return(new_opinions)
}

update_taxonomy <- function(taxon,taxa_modified_after,inc_children=T)	{
#httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,",&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,",&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="")
options(warn=-1);
if (inc_children)	{
	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&variant=all&rel=all_children&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
#	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#	httpTt <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",taxa_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
#	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#	httpTt <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&ops_modified_after=",taxa_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}
pbdb_data <- read.csv(httpTt,header=T,stringsAsFactors = F);
new_taxa <- put_pbdb_dataframes_into_proper_type(pbdb_data);
options(warn=1);
return(new_taxa)
}

clean_pbdb_fields <- function(pbdb_data)	{
dirty_fields <- c("collection_name","collection_aka","state","county","geogcomments","formation","member","stratgroup","county","localsection","regionalsection","stratcomments","lithdescript","geology_comments","author","ref_author","primary_reference","collection_comments","taxonomy_comments","primary_reference","authorizer","enterer","modifier");
d_f <- (1:ncol(pbdb_data))[colnames(pbdb_data) %in% dirty_fields];
for (df in 1:length(d_f))	{
	web_text <- pbdb_data[,d_f[df]];
	if (sum(web_text!="")>0)	pbdb_data[,d_f[df]][web_text!=""] <- pbapply::pbsapply(web_text[web_text!=""],mundify_web_text_dull);
	}
return(pbdb_data);
}

# main function ####
update_pbdb_RData <- function(pbdb_data_list,rock_unit_data,gradstein_2020_emended,paleodb_fixes,do_not_forget=c(),do_not_forget_start="Proterozoic",do_not_forget_end="Phanerozoic")	{
# reload current RData & setup relevant databases ####
pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_sites <- pbdb_data_list$pbdb_sites;
pbdb_sites_refined <- pbdb_data_list$pbdb_sites_refined;
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_opinions <- pbdb_data_list$pbdb_opinions;
pbdb_rocks <- pbdb_data_list$pbdb_rocks;
pbdb_rocks_sites <- pbdb_data_list$pbdb_rocks_sites;

time_scale <- gradstein_2020_emended$time_scale;
finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
zone_database <- gradstein_2020_emended$zones;

fossilworks_collections <- paleodb_fixes$fossilworks_collections;
pbdb_taxonomy_fixes <- paleodb_fixes$paleodb_taxonomy_edits;

pbdb_taxonomy_fixes <- clear_na_from_matrix(pbdb_taxonomy_fixes,"");
pbdb_taxonomy_fixes <- clear_na_from_matrix(pbdb_taxonomy_fixes,0);
pbdb_taxonomy_species_fixes <- pbdb_taxonomy_fixes[pbdb_taxonomy_fixes$taxon_rank %in% c("species","subspecies"),];
pbdb_taxonomy_fixes$created <- gsub(" UTC","",pbdb_taxonomy_fixes$created);
pbdb_taxonomy_fixes$modified <- gsub(" UTC","",pbdb_taxonomy_fixes$modified);

# I hate that I have to do this!!! ####
occurrence_no <- pbdb_finds$occurrence_no[is.na(pbdb_finds$modified)];
fixed_finds <- data.frame(base::t(pbapply::pbsapply(occurrence_no,update_occurrence_from_occurrence_no)));
new_modified <- unlist(fixed_finds$modified);
if (sum(new_modified=="0000-00-00 00:00:00")>0)	{
	new_created <- unlist(fixed_finds$created);
	new_modified[new_modified=="0000-00-00 00:00:00"] <- new_created[new_modified=="0000-00-00 00:00:00"];
	}

pbdb_finds$modified[is.na(pbdb_finds$modified)] <- as.character.Date(new_modified);
#pbdb_finds$modified[pbdb_finds$occurrence_no %in% occurrence_no];

# get new/revised data ####
print("Getting occurrence & locality data entered/modified since the last download");
occs_modified_after <- as.Date(max(pbdb_finds$modified[!is.na(pbdb_finds$modified)]))-1;
print(paste("Looking for data modified after",occs_modified_after));
new_finds <- update_occurrence_data(occs_modified_after = occs_modified_after,basic_environments = c("marine","unknown"),save_files = F);
print(paste("Adding data modified as recently as",max(new_finds$modified)));
if (length(do_not_forget)>0 && sum(do_not_forget=="Metazoa")!=length(do_not_forget))	{
	other_finds <- accersi_occurrence_data(taxa = do_not_forget,onset=do_not_forget_start,end=do_not_forget_end,save_files = F);
	new_finds <- rbind(new_finds,other_finds);
	}
new_finds <- clean_pbdb_fields(new_finds);
new_finds <- new_finds[match(unique(new_finds$occurrence_no),new_finds$occurrence_no),];

colls_modified_after <- as.Date(min(occs_modified_after,max(pbdb_sites$modified[!is.na(pbdb_sites$modified)])));
options(warn=-1);
new_sites <- update_collection_data(colls_modified_after = occs_modified_after,basic_environments = c("marine","unknown"),save_files = F);
coll_id <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites$collection_no]);
paste(length(coll_id),"lost sites");
ns <- 1;
if (length(coll_id)>0)	{
	lost_sites <- accersi_data_for_one_collection(coll_id[ns]);
	} else	{
	lost_sites <- accersi_data_for_one_collection(1);
	lost_sites <- lost_sites[lost_sites$collection_no==0,];
	}
while (ns < length(coll_id))	{
	ns <- ns+1;
	lost_sites <- rbind(lost_sites,accersi_data_for_one_collection(coll_id[ns]));
	}
#other_sites <- accersi_collection_data(taxa = "Ichthyosauromorpha,Plesiosauria,Mosasauria",onset="Permian",end="Cenozoic",save_files = F);
new_sites <- rbind(new_sites,lost_sites);
new_sites <- clean_pbdb_fields(new_sites);
cd <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites$collection_no]);
options(warn=1);

# get taxonomic updates ####
print("Getting taxonomic data entered/modified since the last download");
new_opinions <- update_opinions("Life",ops_modified_after = occs_modified_after);
new_opinions <- clean_pbdb_fields(new_opinions);
taxa_modified_after <- occs_modified_after-5;
#taxa_modified_after <- "2002-04-30";
new_taxa <- update_taxonomy(taxon = "Life",taxa_modified_after = taxa_modified_after);
new_taxa <- clean_pbdb_fields(new_taxa);
#new_taxa[new_taxa$taxon_name=="Isabelinia glabrata",]
new_opinions_spc <- new_opinions[new_opinions$taxon_rank %in% c("species","subspecies"),];
new_opinions_spc <- new_opinions_spc[new_opinions_spc$opinion_type=="class",];
taxon <- new_opinions_spc$taxon_name;
taxon_id <- unique(new_opinions_spc$orig_no);
#pbdb_list <- data.frame(base::t(pbapply::pbsapply(taxon_id,accersi_occurrences_for_one_taxon_no)));
#updated_finds <- list_to_dataframe_for_pbdb_data(pbdb_list);
updated_finds <- accersi_occurrences_for_one_taxon_no(taxon_id[1]);
for (td in 2:length(taxon_id))	updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(taxon_id[td]));
updated_finds$modified[is.na(updated_finds$modified)] <- updated_finds$created[is.na(updated_finds$modified)];

new_finds <- rbind(new_finds,updated_finds[!updated_finds$occurrence_no %in% new_finds$occurrence_no,]);
new_finds <- new_finds[match(unique(new_finds$occurrence_no),new_finds$occurrence_no),];

taxon_name <- new_finds$accepted_name;
informals <- pbapply::pbsapply(taxon_name,revelare_informal_taxa,keep_author_specific=T);
new_finds_spc <- new_finds[!informals,];
dummy_finds <- unique(rbind(pbdb_data_list$pbdb_finds,new_finds_spc));
new_sites_spc <- unique(new_sites[new_sites$collection_no %in% unique(dummy_finds$collection_no),]);
cd <- unique(pbdb_finds$collection_no[!new_sites_spc$collection_no %in% pbdb_sites$collection_no]);
dummy_finds <- NULL;

# repair & redate new sites ####
print("Repairing & Redating the new & updated localities.");
new_sites_spc <- unique(new_sites_spc); # new_sites_spc <- unique(new_sites);
new_sites_spc <- new_sites_spc[order(new_sites_spc$collection_no),];
new_sites_spc <- reparo_unedittable_paleodb_collections(new_sites_spc,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);

if (sum(new_sites_spc$direct_ma>0)>0)	{
	new_sites_redated <- redate_paleodb_collections_with_direct_dates(paleodb_collections=new_sites_spc,finest_chronostrat);
	new_sites_redated <- redate_paleodb_collections_with_time_scale(paleodb_collections=new_sites_redated,time_scale=time_scale,zone_database = zone_database);
	} else	{
	new_sites_redated <- redate_paleodb_collections_with_time_scale(paleodb_collections=new_sites_spc,time_scale=time_scale,zone_database = zone_database);
	}
if (sum(new_sites_redated$zone!="")>0)	new_sites_redated <- redate_paleodb_collections_with_zones(paleodb_collections = new_sites_redated,zone_database = zone_database,time_scale = time_scale,emend_paleodb=T);

coll_no <- unique(sort(c(pbdb_sites$collection_no[is.na(pbdb_sites$paleolat)],pbdb_sites$collection_no[pbdb_sites$geoplate==0])));
if (length(coll_no)>0)	{
	revised_paleogeography <- data.frame(base::t(pbapply::pbsapply(coll_no,paleogeographic_orphanarium)));
	for (cc in 1:ncol(revised_paleogeography))	revised_paleogeography[,cc] <- as.vector(unlist(revised_paleogeography[,cc]));
	revised_paleogeography <- put_pbdb_dataframes_into_proper_type(revised_paleogeography);
	leelas <- match(revised_paleogeography$collection_no,pbdb_sites$collection_no);
	mutants <- match(colnames(revised_paleogeography),colnames(pbdb_sites));
	pbdb_sites[leelas,mutants] <- revised_paleogeography;
	}

if (is.list(pbdb_sites$zone))	pbdb_sites$zone <- as.character(unlist(pbdb_sites$zone));

# separate new finds & sites from recently edited finds & sites ####
print("Integrating modified finds & new finds.");
edited_finds <- new_finds_spc[(1:nrow(new_finds_spc))[new_finds_spc$occurrence_no %in% pbdb_finds$occurrence_no],];
added_finds <- new_finds_spc[(1:nrow(new_finds_spc))[!new_finds_spc$occurrence_no %in% pbdb_finds$occurrence_no],];
pbdb_finds[match(edited_finds$occurrence_no,pbdb_finds$occurrence_no),] <- edited_finds;
pbdb_finds <- rbind(pbdb_finds,added_finds);
pbdb_finds <- pbdb_finds[order(pbdb_finds$collection_no,pbdb_finds$occurrence_no),];

# update taxonomy that PBDB will not
pbdb_finds$accepted_name[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_name[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_name_orig[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_name[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_rank[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_rank[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$subgenus_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$subgenus_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$genus_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$genus_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$genus[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$genus[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$family_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$family_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$family[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$family[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];

age <- new_sites_redated$max_ma;
new_sites_redated$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
age <- new_sites_redated$min_ma;
new_sites_redated$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
new_sites_redated$created <- as.Date(new_sites_redated$created);
new_sites_redated$modified <- as.Date(new_sites_redated$modified);

if (!is.null(pbdb_sites$rock_no))		pbdb_sites$rock_no <- NULL;
if (!is.null(pbdb_sites$formation_no))	pbdb_sites$formation_no <- NULL;
edited_sites <- new_sites_redated[new_sites_redated$collection_no %in% pbdb_sites$collection_no,];
added_sites <- new_sites_redated[(1:nrow(new_sites_redated))[!new_sites_redated$collection_no %in% pbdb_sites$collection_no],];
nn <- match(edited_sites$collection_no,pbdb_sites$collection_no);

edited_sites$created <- as.Date(edited_sites$created);
edited_sites$modified <- as.Date(edited_sites$modified);
added_sites$created <- as.Date(added_sites$created);
added_sites$modified <- as.Date(added_sites$modified);
if (is.null(pbdb_sites$interval_lb))	{
	age <- pbdb_sites$max_ma;
	pbdb_sites$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
	age <- pbdb_sites$min_ma;
	pbdb_sites$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
	}
pbdb_sites[nn,] <- edited_sites;
pbdb_sites <- rbind(pbdb_sites,added_sites);
pbdb_sites <- pbdb_sites[order(pbdb_sites$collection_no),];
nsites <- nrow(pbdb_sites);
pbdb_sites <- reparo_unedittable_paleodb_collections(pbdb_sites,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);
time_scale <- gradstein_2020_emended$time_scale;
time_scale <- rbind(time_scale,finest_chronostrat[!finest_chronostrat$interval %in% time_scale$interval,]);
pbdb_sites <- redate_paleodb_collections_with_time_scale(paleodb_collections=pbdb_sites,time_scale,zone_database);
age <- pbdb_sites$max_ma;
pbdb_sites$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
age <- pbdb_sites$min_ma;
pbdb_sites$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
if (sum(is.na(pbdb_sites$modified))>0)	{
	fix_these1 <- (1:nsites)[is.na(pbdb_sites$created)];
	fix_these2 <- (1:nsites)[is.na(pbdb_sites$modified)];
	fix_these <- sort(unique(c(fix_these1,fix_these2)));
	collection_nos <- pbdb_sites$collection_no[fix_these];
	effed_up_sites <- accersi_collection_data_for_list_of_collection_nos(collection_nos);
	pbdb_sites$created[fix_these] <- effed_up_sites$created;
	pbdb_sites$modified[fix_these] <- effed_up_sites$modified;
	}

pbdb_taxonomy <- rbind(pbdb_taxonomy[!pbdb_taxonomy$taxon_no %in% new_taxa$taxon_no,],new_taxa);
pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% pbdb_taxonomy_fixes$taxon_no,] <- pbdb_taxonomy_fixes[pbdb_taxonomy_fixes$taxon_no %in% pbdb_taxonomy$taxon_no,];
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];

edited_opinions <- new_opinions[(1:nrow(new_opinions))[new_opinions$opinion_no %in% pbdb_opinions$opinion_no],];
added_opinions <- new_opinions[(1:nrow(new_opinions))[!new_opinions$opinion_no %in% pbdb_opinions$opinion_no],];
pbdb_opinions[match(edited_opinions$opinion_no,pbdb_opinions$opinion_no),] <- edited_opinions;
pbdb_opinions <- rbind(pbdb_opinions,added_opinions);
pbdb_opinions <- pbdb_opinions[order(pbdb_opinions$opinion_no),];

# generate PBDB rock unit database ####
print("Generating a database of rock units in the PBDB.");
print("Cleaning stupid rock names yet AGAIN...")
named_rock_unit <- pbdb_sites$formation[pbdb_sites$formation!=""];
pbdb_sites$formation[pbdb_sites$formation!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$member[pbdb_sites$member!=""];
pbdb_sites$member[pbdb_sites$member!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""];
pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$formation_alt[pbdb_sites$formation_alt!=""];
pbdb_sites$formation_alt[pbdb_sites$formation_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$member_alt[pbdb_sites$member_alt!=""];
pbdb_sites$member_alt[pbdb_sites$member_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$stratgroup_alt[pbdb_sites$stratgroup_alt!=""];
pbdb_sites$stratgroup_alt[pbdb_sites$stratgroup_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);

# some legitimate rock names get eradicated this way; put them back if you can!
pbdb_sites <- reparo_unedittable_paleodb_collections(pbdb_sites,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);

pbdb_sites_w_rocks <- unique(rbind(subset(pbdb_sites,pbdb_sites$formation!=""),
								   subset(pbdb_sites,pbdb_sites$stratgroup!=""),
								   subset(pbdb_sites,pbdb_sites$member!="")));
pbdb_sites_w_rocks <- pbdb_sites_w_rocks[order(pbdb_sites_w_rocks$collection_no),];
pbdb_rock_info <- organize_pbdb_rock_data(paleodb_collections = pbdb_sites_w_rocks,geosplit = T);
pbdb_rocks <- pbdb_rock_info$pbdb_rocks;
#write.csv(pbdb_rock_info$pbdb_rocks,"PBDB_Rocks2.csv",row.names = F)
rock_names <- pbdb_rocks$rock_unit_clean_no_rock;
nrocks <- nrow(pbdb_rocks);
group_only <- (1:nrocks)[pbdb_rocks$formation=="" & pbdb_rocks$member==""];
last_group <- (min((1:nrocks)[pbdb_rocks$formation!=""])-1);
formation_names <- c(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],
					 sort(pbdb_rocks$formation_clean_no_rock[pbdb_rocks$formation_clean_no_rock!=""]));
rock_no_orig <- rock_no <- 1:nrocks;
rock_no_sr <- match(pbdb_rocks$rock_unit_clean_no_rock,rock_names);
formation_no <- rock_no_sr[match(pbdb_rocks$formation_clean_no_rock,formation_names)];
formation_no[1:last_group] <- rock_no_sr[match(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],formation_names)];
member_only <- (1:nrocks)[pbdb_rocks$formation=="" & pbdb_rocks$member!=""];
formations_and_members <- member_only[pbdb_rocks$member_clean_no_rock[member_only] %in% formation_names]
formation_no[formations_and_members] <- formation_no[match(pbdb_rocks$member_clean_no_rock[formations_and_members],formation_names)]
pbdb_rocks <- tibble::add_column(pbdb_rocks, formation_no=as.numeric(formation_no), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no_sr=as.numeric(rock_no_sr), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no=as.numeric(rock_no), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no_orig=as.numeric(rock_no_orig), .before = 1);
formation_no_dummy <- match(pbdb_rocks$formation_clean_no_rock,rock_names);
#formation_no_dummy <- match(pbdb_rocks$formation_clean_no_rock,pbdb_rocks$rock_unit_clean_no_rock);
formation_no_dummy[1:last_group] <- match(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],formation_names);
#
add_here <- (1:nrow(pbdb_rocks))[is.na(formation_no_dummy)];
add_formations <- pbdb_rocks$formation_clean_no_rock[is.na(formation_no_dummy)];
#match("Ablakoskovolgy",pbdb_rocks$)
kluge_city <- cbind(add_here,add_formations);#add_formations <- unique(pbdb_rocks$formation_clean_no_rock[is.na(formation_no_dummy)])
kluge_city <- kluge_city[match(unique(add_formations),add_formations),];
add_here <- as.numeric(kluge_city[,1]);
a_f <- length(add_here);
pbdb_rocks_old <- pbdb_rocks;
for (af in a_f:1)	{
	if (af>1)
		while (af>1 & add_here[af]==1+add_here[af-1] & pbdb_rocks$formation_clean_basic[add_here[af]]==pbdb_rocks$formation_clean_basic[add_here[af-1]])
			af <- af-1;
	nn <- add_here[af];
	add_rock <- pbdb_rocks[nn,];
	add_rock$rock_no_orig <- -1;
	add_rock$rock_no <- add_rock$rock_no-0.5;
	add_rock$rock_no_sr <- add_rock$rock_no_sr-0.5;
	add_rock$member <- add_rock$member_clean_basic <- add_rock$member_clean_no_rock <- add_rock$member_clean_no_rock_formal <- "";
	add_rock$full_name <- add_rock$formation;
	add_rock$rock_unit_clean_basic <- add_rock$formation_clean_basic;
	add_rock$rock_unit_clean_no_rock <- add_rock$formation_clean_no_rock;
	add_rock$rock_unit_clean_no_rock_formal <- add_rock$formation_clean_no_rock_formal;
	yin <- 1:(nn-1);
	yang  <- nn:nrow(pbdb_rocks);
	pbdb_rocks <- rbind(pbdb_rocks[yin,],add_rock,pbdb_rocks[yang,]);
#	pbdb_rocks[(nn-5):(nn+2),];
	}
nrocks <- nrow(pbdb_rocks);
pbdb_rocks_redone <- pbdb_rocks;
#nn <- match("Ablakoskovolgy",pbdb_rocks$formation_clean_no_rock)
new_rock_no <- 1:nrocks;
new_rock_no_sr <- new_rock_no[match(pbdb_rocks$rock_no_sr,pbdb_rocks$rock_no_sr)];
new_formation_no <- new_rock_no_sr[match(pbdb_rocks$formation_no,pbdb_rocks$formation_no)];
pbdb_rocks_redone$rock_no <- new_rock_no;
pbdb_rocks_redone$rock_no_sr <- new_rock_no_sr;
pbdb_rocks_redone$formation_no <- new_formation_no;

pbdb_rocks_site_lists <- pbdb_rock_info$site_list;

# eliminate accidental duplicate sites 1####
print("We will now pause to remove replicates from the data...");
dup_colls <- hist(pbdb_sites$collection_no,breaks=0:max(pbdb_sites$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites[pbdb_sites$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
#	if (nrow(dup_sites)>1 && max(dup_sites$pbdb_rock_no)>0)
#		dup_sites <- dup_sites[dup_sites$pbdb_rock_no>0,];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites[pbdb_sites$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites <- unique(pbdb_sites);

dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[(dup_sites$ma_lb-dup_sites$ma_ub)==min(dup_sites$ma_lb-dup_sites$ma_ub),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

cd <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites$collection_no]);
# Refine Rocks by Period ####
print("Refining information about rock units in each geological period.");
effed <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
#sum(pbdb_sites_refined$interval_lb[effed]==pbdb_sites_refined$interval_ub[effed])
effed_int <- pbdb_sites_refined$interval_lb[effed];
pbdb_sites_refined$interval_lb[effed] <- pbdb_sites_refined$interval_ub[effed];
pbdb_sites_refined$interval_ub[effed] <- effed_int;
pbdb_sites_refined$ma_lb[effed] <- finest_chronostrat$ma_lb[match(pbdb_sites_refined$interval_lb[effed],finest_chronostrat$interval)];
pbdb_sites_refined$ma_ub[effed] <- finest_chronostrat$ma_ub[match(pbdb_sites_refined$interval_ub[effed],finest_chronostrat$interval)];

age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
zone_database <- gradstein_2020_emended$zones;
chronostrat_units <- unique(c(pbdb_sites$early_interval,pbdb_sites$late_interval));
hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=chronostrat_units,time_scale=gradstein_2020_emended$time_scale,regional_scale="Stage Slice");
myrs <- sort(unique(round(c(hierarchical_chronostrat$ma_lb,hierarchical_chronostrat$ma_ub),5)),decreasing=T)
hierarchical_chronostrat$bin_first <- match(hierarchical_chronostrat$ma_lb,myrs);
hierarchical_chronostrat$bin_last <- match(hierarchical_chronostrat$ma_ub,myrs)-1;
#write.csv(hierarchical_chronostrat,"Hierarchical_Chronostrat.csv",row.names = F);
#hierarchical_chronostrat$bin_last <- match(hierarchical_chronostrat$bin_last,sort(unique(c(hierarchical_chronostrat$bin_first,hierarchical_chronostrat$bin_last))));
opt_periods <- c("Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene");
#hierarchical_chronostrat$interval[hierarchical_chronostrat$interval!=hierarchical_chronostrat$st] <- hierarchical_chronostrat$st[hierarchical_chronostrat$interval!=hierarchical_chronostrat$st];

for (b1 in 1:(length(opt_periods)))	{
	print(paste("doing the",opt_periods[b1]));
	max_ma <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b1],gradstein_2020_emended$time_scale$interval)];
	if (b1<length(opt_periods))	{
		min_ma <- gradstein_2020_emended$time_scale$ma_ub[match(opt_periods[b1],gradstein_2020_emended$time_scale$interval)];
		} else	{
		min_ma <- 0;
		}
	interval_sites <- pbdb_sites[pbdb_sites$max_ma>=min_ma & pbdb_sites$min_ma<=max_ma,];

	rock_database <- rock_unit_data$rock_unit_database[rock_unit_data$rock_unit_database$ma_ub<(max_ma+25) & rock_unit_data$rock_unit_database$ma_lb>(min_ma-25),];
	rock_to_zone_database <- rock_unit_data$rock_to_zone_database[rock_unit_data$rock_to_zone_database$rock_no %in% rock_database$rock_no,];
	zone_database <- gradstein_2020_emended$zone[gradstein_2020_emended$zone$ma_ub<=(max_ma+25) & gradstein_2020_emended$zone$ma_lb>=(min_ma-25),];
	zone_database <- zone_database[!is.na(zone_database$ma_lb),];

	refined_info <- refine_pbdb_collections_w_external_databases(paleodb_collections=interval_sites,rock_database,zone_database,rock_to_zone_database,finest_chronostrat);
	#colnames(refined_info$refined_collections)[!colnames(refined_info$refined_collections) %in% colnames(pbdb_sites)]
	new_sites <- refined_info$refined_collections[!refined_info$refined_collections$collection_no %in% pbdb_sites_refined$collection_no,];
	old_sites <- refined_info$refined_collections[refined_info$refined_collections$collection_no %in% pbdb_sites_refined$collection_no,];
	old_sites <- old_sites[order(old_sites$collection_no),];
#		now_rocked <- old_sites$collection_no[old_sites$rock_unit_senior!=""];
#		was_rocked <- pbdb_sites_refined$collection_no[pbdb_sites_refined$rock_unit_senior!=""];
#		rock_change <- now_rocked[!now_rocked %in% was_rocked]
	now_rockless <- old_sites$collection_no[!old_sites$rock_unit_senior!=""];
	was_rockless <- pbdb_sites_refined$collection_no[!pbdb_sites_refined$rock_unit_senior!=""];
	change_rock <- now_rockless[!now_rockless %in% was_rockless]
	old_sites <- old_sites[!old_sites$collection_no %in% change_rock,];
#		pbdb_sites_refined$rock_unit_senior[pbdb_sites_refined$collection_no %in% change_rock]
		# update sites
	old_sites <- old_sites[,colnames(old_sites) %in% colnames(pbdb_sites_refined)];
	# only update if there is now rock info!!!!
#	newly_rocked_collections <- old_sites$collection_no[old_sites$rock_no_sr>0];
	old_sites_rocked <- subset(old_sites,old_sites$rock_no>0);
	pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
	pbdb_sites_refined[match(old_sites_rocked$collection_no,pbdb_sites_refined$collection_no),match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))] <- old_sites_rocked[,colnames(old_sites_rocked) %in% colnames(pbdb_sites_refined)];
#	pbdb_sites_refined[match(old_sites_rocked$collection_no,pbdb_sites_refined$collection_no),match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))] <- old_sites_rocked[,match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))];
#	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% old_sites$collection_no,] <- old_sites;

		# add new sites
	new_sites <- new_sites[,colnames(new_sites) %in% colnames(pbdb_sites_refined)];
	if (nrow(new_sites)>0)	{
		if (ncol(new_sites)<ncol(pbdb_sites_refined))	{
			lost_fields <- colnames(pbdb_sites_refined)[!colnames(pbdb_sites_refined) %in% colnames(new_sites)];
			dummy <- array(0,dim=c(nrow(new_sites),length(lost_fields)));
			colnames(dummy) <- lost_fields;
			new_sites <- cbind(new_sites,dummy);
			}
#		if (ncol(pbdb_sites_refined)<=ncol(new_sites))	{
		pbdb_sites_refined <- rbind(pbdb_sites_refined,new_sites[,match(colnames(new_sites),colnames(pbdb_sites_refined))]);
#			}
		}
	pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
	}

dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1 && max(dup_sites$rock_no)>0 && !is.infinite(abs(dup_sites$rock_no)))
		dup_sites <- dup_sites[dup_sites$rock_no>0,];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

effed <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
#sum(pbdb_sites_refined$interval_lb[effed]==pbdb_sites_refined$interval_ub[effed])
effed_int <- pbdb_sites_refined$interval_lb[effed];
pbdb_sites_refined$interval_lb[effed] <- pbdb_sites_refined$interval_ub[effed];
pbdb_sites_refined$interval_ub[effed] <- effed_int;
pbdb_sites_refined$ma_lb[effed] <- finest_chronostrat$ma_lb[match(pbdb_sites_refined$interval_lb[effed],finest_chronostrat$interval)];
pbdb_sites_refined$ma_ub[effed] <- finest_chronostrat$ma_ub[match(pbdb_sites_refined$interval_ub[effed],finest_chronostrat$interval)];

# update rock numbers if need be
rock_database <- rock_unit_data$rock_unit_database;
pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$rock_no>0] <- rock_database$rock_no_sr[match(pbdb_sites_refined$rock_no[pbdb_sites_refined$rock_no>0],rock_database$rock_no)];
beepr:beep("wilhelm");
# add PBDB Rock numbers ####
print("Add numbers to rock units from recently created database.");
p_r_s_l <- length(pbdb_rocks_site_lists);
nsites <- nrow(pbdb_sites_refined);
pbdb_sites_refined$pbdb_formation_no <- pbdb_sites_refined$pbdb_rock_no_sr <- pbdb_sites_refined$pbdb_rock_no <- rep(0,nsites);
#pbdb_rocks_refined$rock_no <- 1:nrow(pbdb_rocks);
for (pr in 1:p_r_s_l)	{
	rd <- match(pr,pbdb_rocks_redone$rock_no_orig);
	pbdb_rows <- match(pbdb_rocks_site_lists[[pr]],pbdb_sites_refined$collection_no);
	pbdb_sites_refined$pbdb_rock_no[pbdb_rows] <- pbdb_rocks_redone$rock_no[rd];
	pbdb_sites_refined$pbdb_rock_no_sr[pbdb_rows] <- pbdb_rocks_redone$rock_no_sr[rd];
	pbdb_sites_refined$pbdb_formation_no[pbdb_rows] <- pbdb_rocks_redone$formation_no[rd];
	}
beepr:beep("wilhelm");
# deal with duds ####
nsites <- nrow(pbdb_sites_refined);
dud_sites <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];

# Optimize by Period ####
print("Optimize localities using limited Unitary Association.");
single_binners <- pbdb_data_list$single_binners;
single_binners_a <- single_binners_z <- c();
age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
single_binners_a <- vector(length=length(opt_periods)-1);
for (b2 in 1:(length(opt_periods)-1))	{
	print(paste("Optimizing the",opt_periods[b2],"+",opt_periods[b2+1]));
	max_ma <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b2],gradstein_2020_emended$time_scale$interval)];
	min_ma <- gradstein_2020_emended$time_scale$ma_ub[match(opt_periods[b2+1],gradstein_2020_emended$time_scale$interval)];
#	max_ma_b <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b2-1],gradstein_2020_emended$time_scale$interval)];
	interval_sites <- pbdb_sites_refined[pbdb_sites_refined$ma_lb<=max_ma & pbdb_sites_refined$ma_ub>=min_ma,];
	single_binners_a[b2] <- sum(interval_sites$interval_lb==interval_sites$interval_ub);
#	interval_sites_b <- pbdb_sites_refined[pbdb_sites_refined$ma_lb<=max_ma_b & pbdb_sites_refined$ma_ub>=min_ma,];
#	interval_finds <- pbdb_finds[pbdb_finds$collection_no %in% interval_sites_b$collection_no,];
	interval_finds <- pbdb_finds[pbdb_finds$collection_no %in% interval_sites$collection_no,];
	interval_finds <- interval_finds[interval_finds$identified_rank %in% c("species","subspecies"),];
	if (is.null(interval_sites$ma_lb))	interval_sites$ma_lb <- interval_sites$max_ma;
	if (is.null(interval_sites$ma_ub))	interval_sites$ma_ub <- interval_sites$min_ma;
	interval_zones <- gradstein_2020_emended$zone[gradstein_2020_emended$zone$ma_ub<=(max_ma+25) & gradstein_2020_emended$zone$ma_lb>=(min_ma-25),];
	options(warn=-1);
	#write.csv(interval_sites,"Ediacaran-Cambrian_Sites.csv",row.names = F)
	interval_finds$early_interval <- interval_sites$early_interval[match(interval_finds$collection_no,interval_sites$collection_no)];
	interval_finds$late_interval <- interval_sites$late_interval[match(interval_finds$collection_no,interval_sites$collection_no)];
	stages <- sort(unique(c(interval_finds$early_interval,interval_finds$late_interval)));
	#match(stages,time_scale$interval)
	# paleodb_finds=interval_finds;paleodb_collections=interval_sites;hierarchical_chronostrat=hierarchical_chronostrat;zone_database=interval_zones;
	optimized_sites <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds=interval_finds,paleodb_collections=interval_sites,hierarchical_chronostrat=hierarchical_chronostrat,zone_database=interval_zones);
	#write.csv(optimized_sites,"Ediacaran-Cambrian_Sites_Refined.csv",row.names = F)
#	keep_these <- match(colnames(optimized_sites)[colnames(optimized_sites) %in% colnames(pbdb_sites)],colnames(optimized_sites));
#	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no,] <- optimized_sites[,keep_these];
	single_binners_z <- c(single_binners_z,sum(optimized_sites$interval_lb==optimized_sites$interval_ub));
#	pbdb_sites_refined$ma_lb[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no] <- optimized_sites$ma_lb;
	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no,] <- optimized_sites;
	options(warn=1);
	}
if (ncol(single_binners)-1 < 10)	{
	new_singles <- paste("single_binners_0",ncol(single_binners)-1,sep="");
	} else	{
	new_singles <- paste("single_binners_",ncol(single_binners)-1,sep="");
	}

age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);

single_binners <- cbind(single_binners,single_binners_z);
colnames(single_binners)[ncol(single_binners)] <- new_singles;
write.csv(single_binners,"Single_Binners.csv",row.names = F);
beepr:beep("wilhelm");

# deal with duds 2 ####
nsites <- nrow(pbdb_sites_refined);
dud_sites <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];

# eliminate accidental duplicate sites 2####
print("Eliminate any duplicate information.");
dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)	dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1 && max(dup_sites$rock_no)>0)	dup_sites <- dup_sites[dup_sites$rock_no>0,];
	if (nrow(dup_sites)>1)	dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
#names(dup_colls) <- (1:max(pbdb_sites_refined$collection_no));
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[(dup_sites$ma_lb-dup_sites$ma_ub)==min(dup_sites$ma_lb-dup_sites$ma_ub),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

# I hate that I have to do this 2!!! ####
print("Make sure that dates are dates & not numbers.");
occurrence_no <- pbdb_finds$occurrence_no[is.na(pbdb_finds$modified)];
fixed_finds <- data.frame(base::t(pbapply::pbsapply(occurrence_no,update_occurrence_from_occurrence_no)));
new_modified <- unlist(fixed_finds$modified);
if (sum(new_modified=="0000-00-00 00:00:00")>0)	{
	new_created <- unlist(fixed_finds$created);
	new_modified[new_modified=="0000-00-00 00:00:00"] <- new_created[new_modified=="0000-00-00 00:00:00"];
	}

pbdb_finds$modified[is.na(pbdb_finds$modified)] <- as.character.Date(new_modified);
#pbdb_finds$modified[pbdb_finds$occurrence_no %in% occurrence_no];

# put a bow on it and gift it to yourself ####
print("You deserve this....");
pbdb_data_list$pbdb_finds <- pbdb_finds;
pbdb_data_list$pbdb_sites <- pbdb_sites;
pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
pbdb_data_list$pbdb_opinions <- pbdb_opinions;
pbdb_data_list$time_scale <- finest_chronostrat;
pbdb_data_list$pbdb_rocks <- pbdb_rocks_redone;
pbdb_data_list$pbdb_rocks_sites <- pbdb_rocks_sites;
pbdb_data_list$single_binners <- single_binners;
save(pbdb_data_list,file=paste(getwd(),"/data/Paleobiology_Database.RData",sep=""));
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Marine_Tetrapod_Diversification/data/Paleobiology_Database.RData",sep=""));

pbdb_data_list_for_class <- list(pbdb_sites_refined,pbdb_finds,pbdb_taxonomy,finest_chronostrat);
names(pbdb_data_list_for_class) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
save(pbdb_data_list_for_class,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_for_Invert_Paleo.RData",sep=""));
}

{}
