#######################################################################################
###                                                                                 ###
###    Copyright (C) 2020,  Burair Alsaihati                                        ###
###                                                                                 ###
###    This program is free software: you can redistribute it and/or modify         ###
###    it under the terms of the GNU General Public License as published by         ###
###    the Free Software Foundation version 3                                       ###
###                                                                                 ###
###    This program is distributed in the hope that it will be useful,              ###
###    but WITHOUT ANY WARRANTY; without even the implied warranty of               ###
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                ###
###    GNU General Public License for more details.                                 ###
###                                                                                 ###
###    You should have received a copy of the GNU General Public License            ###
###    along with this program.  If not, see <https://www.gnu.org/licenses/>.       ###
###                                                                                 ###
###    Email: burair.alsaihati25@uga.edu, burair_99@yahoo.com, szhao@uga.edu        ###
###                                                                                 ###
#######################################################################################


################# Parameters related to the VAF file structure ######################## 
# Please don't change any of these parameters
residue_column_count <- 7;
meta_row_count <- 1;
################# End of parameters related to the VAF file structure ######################## 

depth_cutoff <- 10; # Cutoff for minimum coverage. Any variant with coverage < 10 will be assigned a VAF value of NA.

base_dir <- paste(lab_path, "pancancer\\breed_predict_publish\\sample_files\\", sep="");

################# Input files ######################## 
# make sure to modify the paths to the correct ones
VAF_file <- paste(base_dir, "germline_VAF_matrix.txt.gz", sep="");
depth_file <- paste(base_dir, "germline_depth_matrix.txt.gz", sep="");

################# Output files ########################
# make sure to modify the paths to the correct ones
VAF_output_file <- paste(base_dir, "germline_VAF_matrix.reset_low_coverage.txt.gz", sep="");


################ Main code ############################

VAF_data <- read.table(VAF_file, header=F, sep="\t", check.names=F, stringsAsFactors=F);
depth_data <- read.table(depth_file, header=F, sep="\t", check.names=F, stringsAsFactors=F);


# filtering VAF based on BED depth
for(i in (residue_column_count+1):ncol(depth_data)) {
	depths_as_strings <- as.vector(unlist(depth_data[-c(1:meta_row_count), i]));
	bed_depths_as_numbers <- as.numeric(sapply(depths_as_strings, function(x) {strsplit(x, ",", fixed=T)[[1]][1]},  USE.NAMES=F));
	vcf_depths_as_numbers <- as.numeric(sapply(depths_as_strings, function(x) {strsplit(x, ",", fixed=T)[[1]][2]},  USE.NAMES=F));
	na_depth <- which(is.na(bed_depths_as_numbers));
	vcf_na_depth <- which(is.na(vcf_depths_as_numbers));
	low_depth <- which(bed_depths_as_numbers < depth_cutoff);
	reset_rows <- union(c(na_depth, low_depth), vcf_na_depth) + meta_row_count;
	VAF_data[reset_rows, i] <- "NA";
	print(paste("col", (i-residue_column_count), "of", (ncol(depth_data)-residue_column_count)));
	flush.console();
}

# removing duplicated variants if any
# now make variant names
variant_names <- apply(VAF_data[-c(1:meta_row_count), c(1:5)], MARGIN=1, function(x) {paste(as.vector(unlist(x)), collapse="_")});
names(variant_names) <- NULL;
duplicated_indices <- which(duplicated(variant_names) == TRUE) + meta_row_count;
length(duplicated_indices);
if(length(duplicated_indices) > 0) {
	VAF_data <- VAF_data[-duplicated_indices,];
}

gz_file <- gzfile(VAF_output_file, "w");
write.table(VAF_data, file=gz_file, sep="\t", quote=F, row.names=F, col.names=F);
close(gz_file);
