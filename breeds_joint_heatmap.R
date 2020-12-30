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

library("ComplexHeatmap");

# In the germline VAF file, VAF values for low coverage bases are expected to have NA values
# Varaints with enough coverage that are not called by GATK HaplotypeCaller or fail VariantFilteration should have a VAF value of 0

# These parameters are related to the structure of the input file containing VAF values (don't change any of them)
residue_column_count <- 7; # number of columns describing each variant
meta_row_count <- 1; # number of rows dedicated to meta data (sample ids, and others if applicable) in the VAF input file

############ Script customization parameters ########################
# You may modify these parameters as desired
non_na_percentage_cutoff <- 0.8; # all samples must have known VAF values in at least 80% of the breed-specific variants
examined_breeds <- c("Yorkshire Terrier", "Shih Tzu", "Greyhound", "Golden Retriever");
breed_order <- 1:5; # This will define the order of which heatmap breed color legends will be displayed
cancer_types <- c("MT", "OM", "OSA");
disease_order <- c(1:3); # This will define the order of which heatmap disease color legends will be displayed
# See Glasbey palette "Polychrome: Creating and Assessing Qualitative Palettes With Many Colors" (https://www.biorxiv.org/content/10.1101/303883v1.full)
cancer_pallete <- c("blue", "red", "green", "#9A4D43", "#FF00B6", "#FFD300", "#000032", "#009EFF", "#00FEBE");
# See Kelly palette "Polychrome: Creating and Assessing Qualitative Palettes With Many Colors" (https://www.biorxiv.org/content/10.1101/303883v1.full)
breed_pallete <- c("yellow", "lightblue", "orange", "red", "purple", "gray", "green", "#E58FAC", "#964B00", "black");

############ code dependency paths ########################
# Code to build sample meta data
build_meta_data_code_path <- paste(lab_path, "pancancer\\r-code\\breed_prediction\\publish\\build_sample_meta_data.R", sep="");

############ Input and output paths ########################
# Please modify these file paths as needed
base_dir <- paste(lab_path, "pancancer\\WES\\germline\\germline_residues\\", sep="");
output_base <- paste(base_dir, "plots\\predict_breeds\\publish2\\", sep="");

# Input file containing VAF values for all samples for each germline variant: samples as columns and variants as rows
VAF_input_file <- paste(base_dir, "germline_VAF_matrix.reset_low_coverage.txt.gz", sep="");
# Input file containing all breed-specific variants
specific_variants_file <- paste(output_base, "all_breed_specific_variants.txt", sep="");
# Input file containing all samples meta data
meta_data_file <- paste(lab_path, "pancancer\\metadata\\breed_prediction_metadata.txt", sep="");

output_png1 <- paste(output_base, "breeds_heatmap_main_305_dpi.png", sep=""); # this heatmap won't contain samples with unknown breeds
output_png2 <- paste(output_base, "breeds_heatmap_assignment_305_dpi.png", sep=""); # this heatmap will contain samples with unknown breeds
output_png <- c(output_png1, output_png2);
# output_clusters files will have the list of samples ordered as they appear in the heatmaps (used for supplementary tables and breed validation/prediction results)
output_clusters <- c(paste(output_base, "main_clusters.txt", sep=""), 
			paste(output_base, "assignment_clusters.txt", sep=""));

################ Main code ############################

VAF_data <- read.table(VAF_input_file, header=F, sep="\t", check.names=F, stringsAsFactors=F);
specific_variants_data <- read.table(specific_variants_file, header=T, sep="\t", check.names=F, stringsAsFactors=F);

### building meta_data data frame
source(build_meta_data_code_path);
meta_data <- build_meta_data(meta_data_file);
meta_data <- add_tumor_normal_columns(meta_data, unlist(VAF_data[meta_row_count, ]));

# now make variant names
variant_names <- apply(VAF_data[-c(1:meta_row_count), c(1:6)], MARGIN=1, function(x) {paste(as.vector(unlist(x)), collapse="_")});
names(variant_names) <- NULL;

# now reading VAF values for normal samples
normal_VAF_data <- data.matrix(VAF_data[-c(1:meta_row_count), meta_data[, "NormalCol"]]);
colnames(normal_VAF_data) <- rownames(meta_data);
rownames(normal_VAF_data) <- variant_names;

# now getting the heatmap samples
heatmap_breed_samples <- rownames(meta_data)[which(meta_data[, "Breed"] %in% examined_breeds)];
na_breed_samples <- rownames(meta_data)[which(is.na(meta_data[, "Breed"]) == TRUE)];

heatmap_samples_1 <- c(heatmap_breed_samples); # Samples for first heatmap (no unknown breeds)
heatmap_samples_2 <- c(heatmap_samples_1, na_breed_samples); # samples for second heatmap (with unknown breeds)
sample_list <- list(heatmap_samples_1, heatmap_samples_2);

# now converting breed-specific variants to variant names for the heatmaps
text_tokens_to_variant_name <- function(text_tokens) {
	gene <- text_tokens[1];
	locus_split <- strsplit(text_tokens[2], fixed=TRUE, split=":")[[1]];
	chromosome <- locus_split[1];
	position <- locus_split[2];
	allele_split <- strsplit(text_tokens[3], fixed=TRUE, split=">")[[1]];
	ref_allele <- allele_split[1];
	alt_allele <- allele_split[2];
	protein <- text_tokens[4];
	return(paste(gene, chromosome, position, ref_allele, alt_allele, protein, sep="_"));
}
heatmap_variants <- apply(specific_variants_data, 1, function(x) {text_tokens_to_variant_name(x)});
heatmap_variants <- intersect(heatmap_variants, rownames(normal_VAF_data));

heatmap_data_list <- list();
for(heatmap_version in c(1, 2)) {
	heatmap_samples <- sample_list[[heatmap_version]];
	# building the heatmap data and removing bad samples
	heatmap_data <- normal_VAF_data[heatmap_variants, heatmap_samples];
	VAF_non_na_counts <- apply(heatmap_data, 2, function(x) {length(which(is.na(x) == FALSE))});
	sample_count_cutoff <- length(heatmap_variants) * non_na_percentage_cutoff;
	bad_sample_indices <- which(VAF_non_na_counts < sample_count_cutoff);
	if(length(bad_sample_indices) > 0) {
		heatmap_data <- heatmap_data[, -bad_sample_indices];
	}
	heatmap_data_list[[heatmap_version]] <- heatmap_data;
}

# rm(normal_VAF_data, VAF_data);

# This variable is for debugging only (it will store the sample order for each heatmap)
backup_samples <- list();

for(heatmap_version in c(1,2)) {
	heatmap_data <- heatmap_data_list[[heatmap_version]];
	heatmap_samples <- colnames(heatmap_data);

	# randomly assigning VAF values to low coverage samples (NA)
	for(variant in rownames(heatmap_data)) {
		na_columns <- which(is.na(heatmap_data[variant,]) == TRUE);
		if(length(na_columns) > 0) {
			known_variant_mut_rates <- heatmap_data[variant, -na_columns];
			random_assignments <- sample(known_variant_mut_rates, length(na_columns), replace=TRUE);
			heatmap_data[variant, na_columns] <- random_assignments;
		}
	}

	# defining heatmap annotations and legends
	heatmap_breeds <- c(examined_breeds, "Unknown");
	breed_colors <- breed_pallete[c(1:length(examined_breeds), length(breed_pallete))]; # Unknown is always assigned last "black" color
	names(breed_colors) <- heatmap_breeds;
	breed_colors <- breed_colors[breed_order];
	breed_info <- meta_data[heatmap_samples, "Breed"];
	
	disease_colors <- cancer_pallete[disease_order];
	names(disease_colors) <- cancer_types;
	disease_info <- meta_data[heatmap_samples, "DiseaseAcronym"];
	
	if(heatmap_version == 1) {
		# do nothing
	} else {
		breed_info[which(heatmap_samples %in% na_breed_samples)] <- "Unknown";
	}
	
	annotation_legend_param <- list(labels_gp=gpar(fontsize=20), title_gp=gpar(fontsize=20, fontface="plain"), grid_height=unit(6, "mm"));
	disease_legend_param <- annotation_legend_param;
	disease_legend_param[["ncol"]] <- 2;
	
	top_annotation <- HeatmapAnnotation(col=list(Breed = breed_colors, Disease=disease_colors), Breed=breed_info, Disease=disease_info,
			annotation_legend_param=list(Breed=annotation_legend_param, Disease=disease_legend_param));
	
	set.seed(1000);
	heatmap_object <- Heatmap(name="ht", heatmap_data, col=colorRampPalette(c("white","blue"))(256), top_annotation=top_annotation, 
			show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE, show_row_dend=FALSE, column_dend_height=unit(12, "mm")); # previous 22 mm
	
	png(output_png[heatmap_version], height=3.7, width=9.3, units="in", res=305);
		heatmap_object <- draw(heatmap_object);
		decorate_heatmap_body("ht", {grid.rect(gp = gpar(fill = "transparent", col = "black"))})
	dev.off();

	ordered_samples <- heatmap_samples[column_order(heatmap_object)[[1]]];
	backup_samples[[heatmap_version]] <- ordered_samples;
	temp_dataset <- data.frame("SampleName"=ordered_samples, "Breed"=meta_data[ordered_samples, "Breed"]);
	write.table(temp_dataset, file=output_clusters[heatmap_version], quote=FALSE, sep="\t", row.names = FALSE, col.names=TRUE);
	
}

