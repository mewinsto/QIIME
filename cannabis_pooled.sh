#######################################
## QIIME ANALYSIS FOR POOLED SAMPLES ##
#######################################

##Open-Reference OTU Picking (As advised by Sean)

#pick_subsampled_reference_otus_through_otu_table.py -i /home/mwinston/cannabis/cannabis_pooled/seqs_one.fna,/home/mwinston/cannabis/cannabis_pooled/seqs_two.fna -o /home/mwinston/cannabis/cannabis_pooled/ucrss/ -r /home/mwinston/greengenes/gg_12_10_otus/rep_set/97_otus.fasta -s 0.1

##Summarizes taxa through plots                                                                                                                                                                 
#summarize_taxa_through_plots.py -i ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o ucrss/wf_taxa_summary -m combined_data_mapping_corrected.txt                         

##Simple Beta Diversity Plots (No phylogenetic trees etc.), randomly subsampling OTU-table.biom at 1,700 sequences/sample                                                                   

#beta_diversity_through_plots.py -i ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o bdiv_even1700/ -m combined_data_mapping_corrected.txt -e 1700 -t /home/mwinston/cannabis/cannabis_pooled/ucrss/rep_set.tre

##Computes Alpha Diversity Metrics                                                                                                                                        

#alpha_rarefaction.py -i ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -m combined_data_mapping_corrected.txt -o ucrss/wf_arare/ -p alpha_params.txt -t /home/mwinston/cannabis/cannabis_pooled/ucrss/rep_set.tre        

##Computes OTU category significance for Strain, Soil, and Type -- Must use rarified OTU table or directory of rarified OTU tables                                                            

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Strain -s g_test -o single_g_test_Strain.txt

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Soil -s g_test -o single_g_test_Soil.txt

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Type -s g_test -o single_g_test_Type.txt

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Strain -s ANOVA -o single_ANOVA_Strain.txt

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Soil -s ANOVA -o single_ANOVA_Soil.txt

otu_category_significance.py -i ./bdiv_even1700/otu_table_mc2_w_tax_no_pynast_failures_even1700.biom -m combined_data_mapping_corrected.txt -c Type -s ANOVA -o single_ANOVA_Type.txt


#otu_category_significance.py -i./wf_arare/rarefaction/ -m cannabis_data_corrected.txt -c Strain -s ANOVA -o multiple_anova_Strain.txt -w                                                     

#otu_category_significance.py -i./wf_arare/rarefaction/ -m cannabis_data_corrected.txt -c Soil -s ANOVA -o multiple_anova_Soil.txt -w                                                         

#otu_category_significance.py -i./wf_arare/rarefaction/ -m cannabis_data_corrected.txt -c Type -s ANOVA -o multiple_anova_Type.txt -w

#####################################
##Tests for categorical differences##
#####################################

##Runs RDA for sample type, strain, soil type, and cannabinoids                                                                                                                                                                                      


#compare_categories.py --method dbrda -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o rda_type                                                                                                         

#compare_categories.py --method dbrda -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o rda_strain                                                                                                    

#compare_categories.py --method dbrda -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o rda_soil                                                                                                        

compare_categories.py --method dbrda -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o rda_type_weighted

compare_categories.py --method dbrda -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o rda_strain_weighted

compare_categories.py --method dbrda -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o rda_soil_weighted


#Tests significance of these groups using ADONIS and ANOSIM functions (using 999 reps)                                                                                                                                                                

#compare_categories.py --method anosim -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o anosim -n 999                                                                                                  

#compare_categories.py --method adonis -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o adonis -n 999                                                                                                   

#compare_categories.py --method anosim -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o anosim2 -n 999                                                                                               

#compare_categories.py --method adonis -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o adonis2 -n 999                                                                                               

#compare_categories.py --method adonis -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o anosim3 -n 999                                                                                                

#compare_categories.py --method adonis -i bdiv_even1700/unweighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o adonis3 -n 999                                                                                                

compare_categories.py --method anosim -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o anosim_weighted -n 999

compare_categories.py --method adonis -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Type -o adonis_weighted -n 999

compare_categories.py --method anosim -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o anosim2_weighted -n 999

compare_categories.py --method adonis -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Strain -o adonis2_weighted -n 999

compare_categories.py --method adonis -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o anosim3_weighted -n 999

compare_categories.py --method adonis -i bdiv_even1700/weighted_unifrac_dm.txt -m combined_data_mapping_corrected.txt -c Soil -o adonis3_weighted -n 999


