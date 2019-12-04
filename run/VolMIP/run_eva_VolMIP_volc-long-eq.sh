#! /bin/bash
####################################################################################################
# EVA runscript

EXP=VolMIP_volc-long-eq
SAVEDIR="output"


cd ../..

cat > eva_namelist << EOF
&eva_input
 eruption_list_filename = 'VolMIP_volc-long-eq_volcanic_stratospheric_sulfur_injection.nc',
 parameter_set_filename = 'EVAv1_parameter_set_piBG.nc',
 Lookuptable_filename   = 'eva_Mie_lookuptables.nc'
/
&sulfate_input
 sulfate_start_year     = 1750,
 sulfate_end_year       = 1900
/
&aod_output_files
 aod_output_dir         = '${SAVEDIR}',
 aod_file_savename      = "eva_aod_${EXP}",
 grid_filename          = "./gridfiles/CMIP6_ECHAM6_gridfile.nc",
 aod_start_year         = 1800,
 aod_end_year           = 1850
/
&forcing_output_files
 forcing_output_dir     = '${SAVEDIR}',
 forcing_file_savename  = "eva_forcing_${EXP}",
 grid_filename          = "./gridfiles/CMIP6_ECHAM6_gridfile.nc",
 forcing_start_year     = 1815,
 forcing_end_year       = 1820
/
EOF

./eva_build_sulfate_file
./eva_build_forcing_files_cmip6
./eva_build_aod_file


