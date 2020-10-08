#Coordination script for generating data/figures for CEBP Interception Model paper
library(tidyverse)

#edit this file to label appropriate data versions
source("scripts/current_date_code.R")

#basic data preparation from raw data files
#no need to run this if already pre-processed
{
  source("scripts/prep_input_data/a0_prep_sensitivity_values.R")
  source("scripts/prep_input_data/a1_construct_imputed_SEER_draw.R")
  source("scripts/prep_input_data/a2_join_seer_to_sensitivity.R")
  #the output of this section is a data file containing 
  #per-cancer-type-per-stage-incidence,survival,non-decreasing sensitivity
}

#basic interception model scenarios
#no need to run this if already generated for making figures
{ 
  source("scripts/01_generate_basic_interception_scenarios.R")
}

#generate useful summaries for dwell time groups
#that will be used later for making tables
source("scripts/utility_describe_dwell_groups.R")

#generate tables for paper
source("scripts/02_generate_tables.R")

#generate figure 1 flow diagrams
#involves recomputing flows in detail
source("scripts/03_generate_figure_one_flow_diagram.R")

#generate figure 2 symbolic diagram
source("scripts/04_generate_figure_two_screening_diagram.R")

#generate figure 3 flow diagram
source("scripts/05_generate_figure_three_sankey_plot.R")

##generate supplementals
source("scripts/101_supplemental_SEER_and_sensitivity.R")
source("scripts/102_supplemental_detectable_explanation.R")
#generate modified table if USPSTF A/B approved screening cancer cases cannot be improved
source("scripts/103_supplemental_modified_by_USPSTF_screening.R")
source("scripts/104_supplemental_stage_shift_barplots.R")

#extra analysis
source("scripts/201_auxiliary_power_analysis.R")

#use at top directory to generate zip file if needed
#zip -r -X ExternalInterceptPaper.zip ExternalInterceptPaper -x "*.Rproj" -x "*/\.*