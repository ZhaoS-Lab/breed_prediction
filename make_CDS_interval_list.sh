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

# Run this shell script to make an interval list file for coding sequence (CDS) regions from a GTF format reference genome annotation
# The interval list file is needed for GATK DepthOfCoverage to calculate the read coverage of bases within specific intervals of interest
# Required software: None
# GNU/Linux used: awk

######## expected input arguments ########
# path to GTF format reference genome annotation
gtf_file=Canis_familiaris.CanFam3.1.99.chr.gtf
# path ot output interval list file
output_file=Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS.interval_list

cat $gtf_file | awk -F'\t' '$3 == "CDS" {print $1":"$4"-"$5}' | sort -Vs > $output_file
