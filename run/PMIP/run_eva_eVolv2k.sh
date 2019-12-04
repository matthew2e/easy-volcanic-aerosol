#! /bin/bash
####################################################################################################
# EVA runscript

EXP=eVolv2k_v3

cat > eva_namelist << EOF
&eva_input
 eruption_list_filename = 'eVolv2k_volcanic_stratospheric_sulfur_injection_-500_1900_v3.nc',
 parameter_set_filename = 'EVAv1_parameter_set_piBG.nc',
 Lookuptable_filename   = 'eva_Mie_lookuptables.nc'
/
&sulfate_input
 sulfate_start_year     = -600,
 sulfate_end_year       = 1900
/
&aod_output_files
 aod_output_dir         = '.',
 aod_file_savename      = "eva_aod_${EXP}",
 aod_start_year         = -500,
 aod_end_year           = 1900
/
&forcing_output_files
 forcing_output_dir     = '.',
 forcing_file_savename  = "eva_forcing_${EXP}",
 forcing_start_year     = 1814,
 forcing_end_year       = 1820
/
EOF

./eva_build_sulfate_file
./eva_build_forcing_files_cmip6
./eva_build_aod_file


