#! /bin/bash
####################################################################################################
# EVA runscript

EXP=Pinatubo

cat > eva_namelist << EOF
&eva_input
 eruption_list_filename = 'Eruption_list_GMD_1960_2015.nc',
 parameter_set_filename = 'EVAv1_parameter_set_modernBG.nc',
 Lookuptable_filename   = 'eva_Mie_lookuptables.nc'
/
&sulfate_input
 sulfate_start_year     = 1960,
 sulfate_end_year       = 2015
/
&aod_output_files
 aod_output_dir         = '.',
 aod_file_savename      = "eva_aod_${EXP}",
 aod_start_year         = 1960,
 aod_end_year           = 2010
/
&forcing_output_files
 forcing_output_dir     = '.',
 forcing_file_savename  = "eva_forcing_${EXP}",
 forcing_start_year     = 1990,
 forcing_end_year       = 1995
/
EOF

./eva_build_sulfate_file
./eva_build_forcing_files_cmip6
./eva_build_aod_file


