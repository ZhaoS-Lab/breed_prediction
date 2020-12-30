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

# In the germline VAF file, VAF values for low coverage bases are expected to have NA values
# Varaints with enough coverage that are not called by GATK HaplotypeCaller or fail VariantFilteration should have a VAF value of 0

# These parameters are related to the structure of the input file containing VAF values (don't change any of them)
residue_column_count <- 7; # number of columns describing each variant
meta_row_count <- 1; # number of rows dedicated to meta data (sample ids, and others if applicable) in the VAF input file

######### Parameters for both breed-unique and breed-enriched variants discovery
global_sufficient_cov_cutoff <- 0.8; # all variants must have sufficient coverage in at least 80% of the samples in the discovery dataset
# Breeds for which breed-specific variants, breed validation or prediction to be performed
examined_breeds <- c("Yorkshire Terrier", "Shih Tzu", "Greyhound", "Golden Retriever");
# Samples with Mixed breed will be excluded from all analyses (breed-specific variant discovery, breed validation and prediction)
excluded_breeds <- c("Mixed");
# All samples with known breeds except examined or Mixed are assigned to Other. They will be used to filter out some candidate breed-specific variants, but will not be used in breed validation or prediction.
breed_names <- c(examined_breeds, "Other");
# For breed-unique and breed-enriched variants discovery, any variant with VAF >= 0.2 is considered a non-reference variant (heterozygous or homozygoug alternative)
non_ref_VAF_cutoff <- 0.2;

######### Parameters for breed-unique variants discovery
# For breed-unique variants discovery, VAF interval for heterozygous alleles is [0.2, 0.8) and for homozygoug alternative allele is [0.8, 1]
hom_VAF_cutoff <- 0.8;
breed_sufficient_cov_cutoff <- 0.5; # every unique variant must have sufficient coverage in at least 50% of the samples in every examined breed
unique_variants_sample_cutoffs <- c(5, 0.4); # A breed-unique variant of breed A must be observed in at least 5 samples or 40% of samples in breed A, whichever more. So, this variable includes both cutoffs (count and fraction)  

######### Parameters for breed-enriched variants discovery
fisher_pvalue_cutoff <- 0.1; # A breed-enriched variant must be enriched in a breed against every other breed at Fisher's p-value of 0.1
# All breed-enriched variants must be enriched in a breed "A" against every other breed at fisher_pvalue_cutoff.
# If the below value is TRUE, variants enriched in any comparison between two breeds (neither of which is breed "A") will be filtered out
filter_non_specific_variants <- TRUE;

############ code dependency paths ########################
# Code to build sample meta data
build_meta_data_code_path <- paste(lab_path, "pancancer\\r-code\\breed_prediction\\publish\\build_sample_meta_data.R", sep="");

############ Input and output paths ########################
# Please modify these file paths as needed
base_dir <- paste(lab_path, "pancancer\\WES\\germline\\germline_residues\\", sep="");
# Input file containing VAF values for all samples for each germline variant: samples as columns and variants as rows
VAF_input_file <- paste(base_dir, "germline_VAF_matrix.reset_low_coverage.txt.gz", sep="");
# Input file containing all samples meta data
meta_data_file <- paste(lab_path, "pancancer\\metadata\\breed_prediction_metadata.txt", sep="");

output_base <- paste(base_dir, "plots\\predict_breeds\\publish2\\", sep="");
unique_variants_output_file <- paste(output_base, "breed_unique_variants.txt", sep="");
enriched_variants_output_file <- paste(output_base, "breed_enriched_variants.txt", sep="");
specific_variants_output_file <- paste(output_base, "all_breed_specific_variants.txt", sep="");

############ End of input and output paths ########################


VAF_data <- read.table(VAF_input_file, header=F, sep="\t", check.names=F, stringsAsFactors=F);

### building meta_data data frame
source(build_meta_data_code_path);
meta_data <- build_meta_data(meta_data_file);
meta_data <- add_tumor_normal_columns(meta_data, unlist(VAF_data[meta_row_count, ]));

# now make residue variant names
variant_names <- apply(VAF_data[-c(1), c(1:6)], MARGIN=1, function(x) {paste(as.vector(unlist(x)), collapse="_")});
names(variant_names) <- NULL;

variant_name_to_output <- function(variant_name) {
	variant_tokens <- strsplit(variant_name, fixed=TRUE, split="_")[[1]];
	gene <- variant_tokens[1];
	locus <- paste(variant_tokens[2], variant_tokens[3], sep=":");
	mutation <- paste(variant_tokens[4], variant_tokens[5], sep=">");
	protein <- variant_tokens[6];
	return(paste(gene, locus, mutation, protein, sep="\t"));
}

# now reading VAF values for normal samples
normal_VAF_data <- data.matrix(VAF_data[-c(1:meta_row_count), meta_data[, "NormalCol"]]);
colnames(normal_VAF_data) <- rownames(meta_data);
rownames(normal_VAF_data) <- variant_names;

# now get samples for each breed
breed_info <- meta_data[, "Breed"];
names(breed_info) <- rownames(meta_data);
breed_info[which(breed_info %in% excluded_breeds)] <- NA;

breed_samples <- list();
for(breed in examined_breeds) {
	breed_samples[[breed]] <- names(which(breed_info == breed));
}


# Check first to see if there are samples with pure breed other than the examined breeds
other_breed_samples <- setdiff(names(breed_info), names(which(is.na(breed_info))));
other_breed_samples <- setdiff(breed_samples[["Other"]], names(breed_info)[which(breed_info %in% examined_breeds)]);
if(length(other_breed_samples) > 0) {
	breed_samples[["Other"]] <- setdiff(names(breed_info), names(which(is.na(breed_info))));
	breed_samples[["Other"]] <- setdiff(breed_samples[["Other"]], names(breed_info)[which(breed_info %in% examined_breeds)]);
} else {
	# There are no samples with a pure breed other than the examined breed. We must remove "Other" breed category
	breed_names <- setdiff(breed_names, "Other");
}


non_na_count_cutoff <- ncol(normal_VAF_data) * global_sufficient_cov_cutoff;
############ Finding breed-unique variants ########################
unique_variant_results <- paste(c("Gene", "Locus", "Mutation", "Protein", "Breed", "Variant_Type", "Sample_Count", "Sample_Percentage"), collapse="\t");
unique_variants_summary <- c();
unique_sample_counts <- list("Hom"=list(), "Het"=list(), "Het,Hom"=list(), "Any"=list());
for(variant_type_counts in names(unique_sample_counts)) {
	for(breed in examined_breeds) {
		unique_sample_counts[[variant_type_counts]][[breed]] <- list("Count"=c(), "Percentage"=c());
	}
}

for(i in 1:nrow(normal_VAF_data)) {
	variant <- variant_names[i];
	variant_VAF <- normal_VAF_data[variant,];
	non_na_samples <- names(which(is.na(variant_VAF) == FALSE));
	if(length(non_na_samples) >= non_na_count_cutoff) {
		variant_VAF <- variant_VAF[non_na_samples];
		ref_samples <- names(which(variant_VAF < non_ref_VAF_cutoff));
		homoz_samples <- names(which(variant_VAF >= hom_VAF_cutoff));
		heter_samples <- names(which(variant_VAF >= non_ref_VAF_cutoff));
		heter_samples <- setdiff(heter_samples, homoz_samples);
		
		samples_by_variant_type <- list("Ref"=ref_samples, "Het"=heter_samples, "Hom"=homoz_samples);
		breeds_by_variant_type <- list();
		for(variant_type in names(samples_by_variant_type)) {
			breeds_by_variant_type[[variant_type]] <- setdiff(unique(meta_data[samples_by_variant_type[[variant_type]], "Breed"]), NA);
		}
		
		## check uniqueness (exactly one breed is observed in Het and/or Hom)
		unique_variant_type <- c();
		if(length(breeds_by_variant_type[["Hom"]]) == 1 && length(breeds_by_variant_type[["Het"]]) == 1) {
			# is it the same unique breed that has both homozygoug and heterozygous variants
			if(breeds_by_variant_type[["Hom"]][1] == breeds_by_variant_type[["Het"]][1]) {
				unique_variant_type <- c("Het", "Hom");
			}
		} else if(length(breeds_by_variant_type[["Hom"]]) == 1 && length(breeds_by_variant_type[["Het"]]) == 0) {
			unique_variant_type <- c("Hom");
		} else if(length(breeds_by_variant_type[["Hom"]]) == 0 && length(breeds_by_variant_type[["Het"]]) == 1) {
			unique_variant_type <- c("Het");
		} # else unique_variant_type is empty, meaning this variant is not unique to any breed
		
		if(length(unique_variant_type) > 0) {
			## this variant is unique to one breed
			## Check condition 1
			## 1. if there are enough samples this breed with the variant and if it belongs to an examined breed
			breed_of_unique_variant <- breeds_by_variant_type[[unique_variant_type[1]]];
			if(breed_of_unique_variant %in% examined_breeds) {
				breed_samples_with_variant <- c();
				for(variant_type in unique_variant_type) {
					breed_samples_with_variant <- c(breed_samples_with_variant, samples_by_variant_type[[variant_type]]);
				}
				breed_samples_with_variant <- intersect(breed_samples_with_variant, breed_samples[[breed_of_unique_variant]]);
				condition1 <- length(breed_samples_with_variant) >= unique_variants_sample_cutoffs[1] && 
						length(breed_samples_with_variant)/length(breed_samples[[breed_of_unique_variant]]) >= unique_variants_sample_cutoffs[2];
				
				if(condition1 == TRUE) {
					## There are enough samples in this breed with Hom and/or Het variant(s).
					## 2. Check if every other breed has enough Ref variant samples
					condition2 <- TRUE;
					condition2_iteration <- 0;
					while(condition2 == TRUE && condition2_iteration < length(breed_samples)) {
						condition2_iteration <- condition2_iteration + 1;
						other_breed <- names(breed_samples)[condition2_iteration];
						if(other_breed != breed_of_unique_variant) {
							other_breed_ref_samples <- intersect(ref_samples, breed_samples[[other_breed]]);
							if(length(other_breed_ref_samples)/length(breed_samples[[other_breed]]) < breed_sufficient_cov_cutoff) {
								condition2 <- FALSE;
							}
						}
					}
					
					if(condition2 == TRUE) {
						## This is a breed-unique variant with all conditions satisfied.
						## Add it to the final results
						unique_variants_summary[variant] <- breed_of_unique_variant;
						summary_type <- paste(unique_variant_type, collapse=",");
						sample_count <- length(breed_samples_with_variant);
						sample_percentage <- round(100 * length(breed_samples_with_variant)/length(breed_samples[[breed_of_unique_variant]]), 2);
						summary_row <- c(variant_name_to_output(variant), breed_of_unique_variant, summary_type, sample_count, sample_percentage);
						unique_variant_results <- c(unique_variant_results, paste(summary_row, collapse="\t"));
					}
				}
			}
		}
	}
	
	if(i %% 1000 == 0) {
		print(paste("Testing for breed-unique variants:", i, "of", length(variant_names)));
		flush.console();
	}
}

# Saving results for breed-unique variants
sink(unique_variants_output_file);
	for(line in unique_variant_results) {
		cat(line, "\n", sep="");
	}
sink();


############ Finding breed-enriched variants ########################
enriched_variants_summary <- c();
enriched_variant_results <- paste("Gene", "Locus", "Mutation", "Protein", "Enriched in", sep="\t");

for(i in 1:nrow(normal_VAF_data)) {
	variant <- variant_names[i];
	variant_VAF <- normal_VAF_data[variant,];
	non_na_samples <- names(which(is.na(variant_VAF) == FALSE));
	if(length(non_na_samples) >= non_na_count_cutoff) {
		enrichment_matrix <- matrix(nrow=length(breed_names), ncol=length(breed_names));
		rownames(enrichment_matrix) <- breed_names;
		colnames(enrichment_matrix) <- breed_names;
		
		for(b1 in 1:(length(breed_names)-1)) {
			breed1 <- breed_names[b1];
			b1_samples <- intersect(breed_samples[[breed1]], non_na_samples);
			b1_samples_with_variant <- which(variant_VAF[b1_samples] >= non_ref_VAF_cutoff);
			b1_variant_ratio <- length(b1_samples_with_variant)/length(b1_samples);
			
			for(b2 in (b1+1):length(breed_names)) {
				breed2 <- breed_names[b2];
				b2_samples <- intersect(breed_samples[[breed2]], non_na_samples);
				# This is fisher test (we shouldn't use fisher test which is why it's commented out)
				b2_samples_with_variant <- which(variant_VAF[b2_samples] >= non_ref_VAF_cutoff);
				x <- cbind(c(length(b1_samples_with_variant), length(b1_samples) - length(b1_samples_with_variant)),
							c(length(b2_samples_with_variant), length(b2_samples) - length(b2_samples_with_variant)));
				fisher_test <- fisher.test(x);
				if(is.na(fisher_test[["p.value"]]) == FALSE && fisher_test[["p.value"]] < fisher_pvalue_cutoff) {
					b2_variant_ratio <- length(b2_samples_with_variant)/length(b2_samples);
					if(b1_variant_ratio > b2_variant_ratio) {
						enrichment_matrix[breed1, breed2] <- "Enriched";
					} else {
						enrichment_matrix[breed2, breed1] <- "Enriched";
					}
				}
			}
		}
		
		breed_enrichment_counts <- apply(enrichment_matrix, 1, function(x) {length(which(x == "Enriched"))});
		names(breed_enrichment_counts) <- breed_names;
		
		if(filter_non_specific_variants == TRUE) {
			# Only allow breed-specific variants
			if(length(which(breed_enrichment_counts == 0)) == (length(breed_names)-1)) {
				# Only one breed has been reported enriched in at least one fisher test
				enriched_breed <- names(which(breed_enrichment_counts > 0));
				if(breed_enrichment_counts[enriched_breed] == (length(breed_names)-1)) {
					# This breed is enriched in this variant against every other breed.
					enriched_variants_summary[variant] <- enriched_breed;
					enriched_variant_results <- c(enriched_variant_results, paste(variant_name_to_output(variant), enriched_breed, sep="\t"));
				}
			}
		} else {
			# Allow all breed-enriched variants
			enriched_breed <- which(breed_enrichment_counts == (length(breed_names)-1));
			if(length(enriched_breed) > 0) {
				enriched_breed <- names(breed_enrichment_counts)[enriched_breed];
				enriched_variants_summary[variant] <- enriched_breed;
				enriched_variant_results <- c(enriched_variant_results, paste(variant_name_to_output(variant), enriched_breed, sep="\t"));
			}
		}
		
	}
	
	if(i %% 250 == 0) {
		print(paste("Testing for breed-enriched variants:", i, "of", length(variant_names)));
		flush.console();
	}
}

# Saving results for breed-unique variants
sink(enriched_variants_output_file);
	for(line in enriched_variant_results) {
		cat(line, "\n", sep="");
	}
sink();

# Saving results for all breed-specific variants
breed_unique_only <- setdiff(names(unique_variants_summary), names(enriched_variants_summary));
all_specific_variants <- enriched_variants_summary;
if(length(breed_unique_only) > 0) {
	all_specific_variants <- c(all_specific_variants, breed_unique_only);
}
all_specific_variants <- all_specific_variants[order(all_specific_variants)];

sink(specific_variants_output_file);
	# write the header line first
	cat(enriched_variant_results[1], "\n", sep="");
	for(i in 1:length(all_specific_variants)) {
		line <- paste(variant_name_to_output(names(all_specific_variants[i])), all_specific_variants[i], sep="\t");
		cat(line, "\n", sep="");
	}
sink();
