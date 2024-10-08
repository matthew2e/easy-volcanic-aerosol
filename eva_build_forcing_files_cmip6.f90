!! @brief Program to create volcanic aerosol forcing files with the EVA module 
!!
!! @author Matthew Toohey, MPI-M, Hamburg, GEOMAR, Kiel (2016-11-03)
!!
!! @par Copyright
!! 

PROGRAM eva_build_forcing_file_cmip6 

  USE mo_eva
  USE netcdf

  IMPLICIT NONE

  INTEGER :: &
       forcing_start_year = 1990, &
       forcing_end_year   = 1995

  CHARACTER(len=50) :: grid_filename         = "CMIP6_ECHAM6_gridfile.nc"
  CHARACTER(len=50) :: sulfate_filename      = "eva_sulfate_timeseries.nc"
  CHARACTER(len=100) :: forcing_output_dir    = "."
  CHARACTER(len=100) :: forcing_file_savename = "eva_forcing_echam_T63_sw"
  CHARACTER(len=5)  :: yearstr

  INTEGER, PARAMETER :: nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)

  INTEGER ::                &
       nplume             , & ! number of plume
       nSO4               , & ! number of timesteps in SO4 file
       nyear              , & ! number of years to save forcing files
       ntime              , & ! number of time steps (nyear*12)
       nz                 , & ! number of vertical levels
       nlat               , & ! number of latitudes
       nsw                , & ! number of shortwave wavelengths
       nlw                , & ! number of longwave wavelengths
       iret               , & ! netCDF reading return variable
       ncid               , & ! netCDF file ID
       sulfate_ncid       , & ! netcdf file ID for sulfate file
       VarID              , & ! pointer to generic dimension in netCDF file
       timeID             , & ! pointer to time dimension in netCDF file
       zID                , & ! pointer to z dimension in netCDF file
       latID              , & ! pointer to latitude dimension in netCDF file
       wlID               , & ! pointer to wavelength dimension in netCDF file
       lwID               , & ! pointer to longwave wavelength dimension in netCDF file
       var_t_ID           , & ! pointer to time variable in netCDF file
       var_dy_ID          , & ! pointer to decimal year variable in netCDF file
       var_z_ID           , & ! pointer to z variable in netCDF file
       var_lat_ID         , & ! pointer to lat variable in netCDF file
       var_sw_ID          , & ! pointer to wavelenth variable in netCDF file
       var_wl1_sw_ID      , & ! pointer to wavelenth variable in netCDF file
       var_wl2_sw_ID      , & ! pointer to wavelenth variable in netCDF file
       var_ext_sw_ID      , & ! pointer to extinction variable in netCDF file
       var_ssa_sw_ID      , & ! pointer to SSA variable in netCDF file
       var_asy_sw_ID      , & ! pointer to ASY variable in netCDF file
       var_wl1_lw_ID      , & ! pointer to wavelenth variable in netCDF file
       var_wl2_lw_ID      , & ! pointer to wavelenth variable in netCDF file
       var_ext_lw_ID      , & ! pointer to extinction variable in netCDF file
       var_ssa_lw_ID      , & ! pointer to SSA variable in netCDF file
       var_asy_lw_ID      , & ! pointer to ASY variable in netCDF file
       var_aod_sw_ID      , & ! pointer to aod variable in netCDF file
       var_aod_lw_ID      , & ! pointer to aod variable in netCDF file
       var_reff_ID        , & ! pointer to reff variable in netCDF file
       SO4_i0             , & ! index of time in SO4 file corresponding to first time of forcing file
       ind(12)            , & ! indices referencing years in full timeseries
       this_year          , & ! a single year, used in num2str 
       mons_since_1850(12), & ! for time variable consistent with CMIP6 historical
       iyear              , & ! an index
       imonth             , & ! 
       i                  , & ! an index
       j                  , & ! another index
       y                      ! one more

  !LOGICAL :: save_yearly = .true.   ! Set to false to save all data in one file         

  REAL, ALLOCATABLE ::       &
       year(:)             , & ! years
       month(:)            , & ! months
       fyear(:)            , & ! fractional year
       SO4(:,:)            , & ! sulfate data, loaded from file
       SO4_time(:)         , & ! time of sulfate data
       lat(:)              , & ! latitude of output
       reff  (:,:,:)       , & ! effective radius
       reff_vec(:)         , & ! dummy variable for call to EVA routine

       wl1_sw(:)           , & ! lower bound wavelength of shortwave spectral bands
       wl2_sw(:)           , & ! upper bound wavelength of shortwave spectral bands
       lambda_sw(:)        , & ! wavelengths of output
       ext550_sw(:,:,:)    , & ! extinction at 550 nm
       ext_sw(:,:,:,:)     , & ! extinction(time,z,lat,wl)
       ssa_sw(:,:,:,:)     , & ! ssa(time,z,lat,wl)
       asy_sw(:,:,:,:)     , & ! asy(time,z,lat,wl)
       aod_sw(:,:,:)       , & ! AOD, (time, lat, wl)
       ext550_sw_vec(:)    , & ! dummy variable for call to EVA routine
       ext_sw_vec(:,:)     , & ! dummy variable for call to EVA routine
       ssa_sw_vec(:,:)     , & ! dummy variable for call to EVA routine
       asy_sw_vec(:,:)     , & ! dummy variable for call to EVA routine
       aod_sw_dum(:)       , & ! dummy variable for call to EVA routine
       wl1_lw(:)           , & ! lower bound wavelength of longwave spectral bands
       wl2_lw(:)           , & ! upper bound wavelength of longwave spectral bands
       lambda_lw(:)        , & ! wavelengths of output
       ext550_lw(:,:,:)    , & ! extinction at 550 nm
       ext_lw(:,:,:,:)     , & ! extinction(time,z,lat,wl)
       ssa_lw(:,:,:,:)     , & ! ssa(time,z,lat,wl)
       asy_lw(:,:,:,:)     , & ! asy(time,z,lat,wl)
       aod_lw(:,:,:)       , & ! AOD, (time, lat, wl)
       ext550_lw_vec(:)    , & ! dummy variable for call to EVA routine
       ext_lw_vec(:,:)     , & ! dummy variable for call to EVA routine
       ssa_lw_vec(:,:)     , & ! dummy variable for call to EVA routine
       asy_lw_vec(:,:)     , & ! dummy variable for call to EVA routine
       aod_lw_dum(:)        ! dummy variable for call to EVA routine       
 
REAL ::                    &
       z(70)               , & ! height grid for EVA calculations
       so4_in(3)               ! sulfate triple to use in time loop

  character(8)  :: date
  character(10) :: time
  LOGICAL :: signed_year_output=.false.  ! signed years useful when extending into BCE
  real :: add_to_year_output=0        ! Sometimes used to avoid negative years in filename

  ! Input parameters list
  NAMELIST /FORCING_OUTPUT_FILES/ forcing_output_dir, forcing_file_savename, grid_filename, forcing_start_year, forcing_end_year
  OPEN (UNIT=30, FILE='eva_namelist', STATUS='OLD')
  READ (30, NML=FORCING_OUTPUT_FILES)
  CLOSE (30)

  ! Define time grid

  write(*,*) 'Building forcing'

  nyear=forcing_end_year-forcing_start_year+1
  ntime=nyear*12

  ALLOCATE(fyear(12))
  ALLOCATE(year(nyear))
  !ALLOCATE(month(ntime))


  i=1

  do iyear=1,nyear
    !do imonth=1,nmon
      year(i)  = forcing_start_year+iyear-1
    !  month(i) = mons(imonth)
      i=i+1
    !end do
  end do

  !fyear=year+(month-1)/12.

  ! Define EVA height grid
  z = (/ (i, i = 10, 79, 1) /) / 2.0
  nz = size(z)
 
  ! Read grid file for lat and wavelengths
  iret = nf90_open(TRIM(grid_filename), NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening grid file'
  iret = nf90_inq_dimid(ncid, "latitude", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlat)
  iret = nf90_inq_dimid(ncid, "solar_bands", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nsw)
  iret = nf90_inq_dimid(ncid, "terrestrial_bands", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlw)

  allocate(lat(nlat))
  allocate(reff(nz,nlat,12))
  allocate(reff_vec(nz))

  allocate(ext550_sw(nz,nlat,12))
  allocate(ext_sw(nz,nlat,nsw,12))
  allocate(ssa_sw(nz,nlat,nsw,12))
  allocate(asy_sw(nz,nlat,nsw,12))
  allocate(wl1_sw(nsw))
  allocate(wl2_sw(nsw))
  allocate(lambda_sw(nsw))
  allocate(ext550_sw_vec(nz))
  allocate(ext_sw_vec(nz,nsw))
  allocate(ssa_sw_vec(nz,nsw))
  allocate(asy_sw_vec(nz,nsw))
  allocate(aod_sw_dum(nsw))
  allocate(aod_sw(nsw,nlat,12))

  allocate(ext550_lw(nz,nlat,12))
  allocate(ext_lw(nz,nlat,nlw,12))
  allocate(ssa_lw(nz,nlat,nlw,12))
  allocate(asy_lw(nz,nlat,nlw,12))
  allocate(wl1_lw(nlw))
  allocate(wl2_lw(nlw))
  allocate(lambda_lw(nlw))
  allocate(ext550_lw_vec(nz))
  allocate(ext_lw_vec(nz,nlw))
  allocate(ssa_lw_vec(nz,nlw))
  allocate(asy_lw_vec(nz,nlw))
  allocate(aod_lw_dum(nlw))
  allocate(aod_lw(nlw,nlat,12))


  iret = nf90_inq_varid(ncid, "latitude", VarID)
  iret = nf90_get_var(ncid, VarID, lat(:)  , start=(/1/) ,count=(/nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile latitude'

  iret = nf90_inq_varid(ncid, "wl1_sun", VarID)
  iret = nf90_get_var(ncid, VarID, wl1_sw(:)  , start=(/1/) ,count=(/nsw/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile wavelengths'

  iret = nf90_inq_varid(ncid, "wl2_sun", VarID)
  iret = nf90_get_var(ncid, VarID, wl2_sw(:)  , start=(/1/) ,count=(/nsw/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile wavelengths'

  iret = nf90_inq_varid(ncid, "wl1_earth", VarID)
  iret = nf90_get_var(ncid, VarID, wl1_lw(:)  , start=(/1/) ,count=(/nlw/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile wavelengths'

  iret = nf90_inq_varid(ncid, "wl2_earth", VarID)
  iret = nf90_get_var(ncid, VarID, wl2_lw(:)  , start=(/1/) ,count=(/nlw/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile wavelengths'

  lambda_sw=( wl1_sw + wl2_sw ) / 2.0
  lambda_lw=( wl1_lw + wl2_lw ) / 2.0

  !write(*,*) lambda

  ! Input sulfate data

  iret = nf90_open(TRIM(sulfate_filename), NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening sulfate file'
  sulfate_ncid=ncid
  iret = nf90_inq_dimid(ncid, "time", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nSO4)
  iret = nf90_inq_dimid(ncid, "plume", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nplume)

  allocate(SO4(nplume,nSO4))
  allocate(SO4_time(nSO4))

  iret = nf90_inq_varid(ncid, "time", VarID)
  iret = nf90_get_var(ncid, VarID, SO4_time(:)  , start=(/1/) ,count=(/nSO4/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading sulfate file'

  iret = nf90_inq_varid(ncid, "SO4", VarID)
  iret = nf90_get_var(ncid, VarID, SO4(:,:)  , start=(/1,1/) ,count=(/nplume,nSO4/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading sulfate file'

  do y=1,nyear

     ! find index on SO4_time for t0

     SO4_i0=1
     do while ( SO4_time(SO4_i0) < year(y) )
       SO4_i0 = SO4_i0 + 1
     end do

     ! for each timestep
     do i=1,12
       so4_in=SO4(:,SO4_i0+i-1)
       do j=1,nlat
         call eva_ext_reff(lat,z,lat(j),so4_in,nlat,nz,ext550_sw_vec,reff_vec)
         call eva_aop_profile(lat,z,lat(j),so4_in,lambda_sw,nsw,nlat,nz,ext_sw_vec,ssa_sw_vec,asy_sw_vec)
         call eva_aod(lat,lat(j),so4_in,lambda_sw,nlat,nsw,nz,aod_sw_dum)

         ext550_sw(:,j,i)= ext550_sw_vec
         reff(:,j,i)  = reff_vec
         ext_sw(:,j,:,i) = ext_sw_vec
         ssa_sw(:,j,:,i) = ssa_sw_vec
         asy_sw(:,j,:,i) = asy_sw_vec
         aod_sw(:,j,i) = aod_sw_dum
         call eva_aop_profile(lat,z,lat(j),so4_in,lambda_lw,nlw,nlat,nz,ext_lw_vec,ssa_lw_vec,asy_lw_vec)
         call eva_aod(lat,lat(j),so4_in,lambda_lw,nlat,nlw,nz,aod_lw_dum)
         ext550_lw(:,j,i)= ext550_lw_vec
         ext_lw(:,j,:,i) = ext_lw_vec
         ssa_lw(:,j,:,i) = ssa_lw_vec
         asy_lw(:,j,:,i) = asy_lw_vec
         aod_lw(:,j,i) = aod_lw_dum
       end do
     end do

     !write(*,*) 'Writing forcing files'
 
     call date_and_time(DATE=date,TIME=time)
  
     ! save data in netcdf file
    !do i=1,nyear
    !ind=y*12+(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)

    !this_year=year(1+(i-1)*12)
    this_year=year(y)
    fyear=this_year+(mons-1)/12.0
    mons_since_1850=(mons-1)+(this_year-1850)*12
    this_year=this_year+add_to_year_output
    if (signed_year_output) then
       write ( yearstr , '(SP, i5.4)' ) this_year
    else
       write ( yearstr, '( i4.4)' ) this_year
    end if

    write(*,*) 'Writing forcing file: ', TRIM(forcing_output_dir)//'/'//TRIM(forcing_file_savename)//'_'//trim(yearstr)//'.nc'

    iret = NF90_NOERR
    iret = iret + nf90_create(TRIM(forcing_output_dir)//'/'//TRIM(forcing_file_savename)//'_'//trim(yearstr)//'.nc', &
        NF90_CLOBBER, ncid)
    iret = iret + nf90_def_dim(ncid, 'month' ,NF90_UNLIMITED    , timeID)
    iret = iret + nf90_def_dim(ncid, 'altitude'       ,nz    , zID)
    iret = iret + nf90_def_dim(ncid, 'latitude'       ,nlat  , latID)
    iret = iret + nf90_def_dim(ncid, 'solar_bands'    ,nsw   , wlID)
    iret = iret + nf90_def_dim(ncid, 'terrestrial_bands'    ,nlw   , lwID)
    IF (iret /= 6*NF90_NOERR) STOP 'Error in Creating File Dimensions'
    !
    iret = NF90_NOERR
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"title","EVA v1.2: stratospheric aerosol optical properties")
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"EVA_reference","Toohey, M., Stevens, B., Schmidt, H., and Timmreck,&
              C.: Easy Volcanic Aerosol (EVA v1.0): an idealized forcing generator for climate simulations, Geosci. &
              Model Dev., 9, 4049-4070, https://doi.org/10.5194/gmd-9-4049-2016, 2016.")
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"EVA_source_code","https://github.com/matthew2e/easy-volcanic-aerosol")
    iret = iret + nf90_copy_att(sulfate_ncid,NF90_GLOBAL,"input_vssi_file",ncid,NF90_GLOBAL)
    iret = iret + nf90_copy_att(sulfate_ncid,NF90_GLOBAL,"input_sulfate_parameter_file",ncid,NF90_GLOBAL)
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"input_forcing_parameter_file",TRIM(parameter_set_filename))
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"input_grid_file",TRIM(grid_filename))
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"input_Mie_file",TRIM(Lookuptable_filename))
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date(7:8)//'.'//date(5:6)//'.' &
                                       //date(1:4)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:6))
    IF (iret /= 7*NF90_NOERR) STOP 'Error in Creating File Attributes'

    iret = NF90_NOERR
    iret = iret + nf90_def_var(ncid, 'month'        , NF90_FLOAT, timeID,  var_t_ID)
    iret = iret + nf90_def_var(ncid, 'altitude'     , NF90_FLOAT, zID,     var_z_ID)
    iret = iret + nf90_def_var(ncid, 'latitude'     , NF90_FLOAT, latID,   var_lat_ID)
    iret = iret + nf90_def_var(ncid, 'dec_year'     , NF90_FLOAT, timeID,  var_dy_ID)
!    iret = iret + nf90_def_var(ncid, 'solar_bands'  , NF90_FLOAT, wlID,    var_sw_ID)
    iret = iret + nf90_def_var(ncid, 'wl1_sun'      , NF90_FLOAT, wlID,    var_wl1_sw_ID)
    iret = iret + nf90_def_var(ncid, 'wl2_sun'      , NF90_FLOAT, wlID,    var_wl2_sw_ID)
    iret = iret + nf90_def_var(ncid, 'wl1_earth'    , NF90_FLOAT, lwID, var_wl1_lw_ID)
    iret = iret + nf90_def_var(ncid, 'wl2_earth'    , NF90_FLOAT, lwID, var_wl2_lw_ID)
    iret = iret + nf90_def_var(ncid, 'ext_sun'      , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_ext_sw_ID)
    !print *, trim(NF90_STRERROR(iret)) 
    iret = iret + nf90_def_var(ncid, 'omega_sun'    , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_ssa_sw_ID)
    iret = iret + nf90_def_var(ncid, 'g_sun'        , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_asy_sw_ID)
    iret = iret + nf90_def_var(ncid, 'ext_earth'    , NF90_FLOAT, (/zID,latID,lwID,timeID/), var_ext_lw_ID)
    !print *, trim(NF90_STRERROR(iret)) 
    iret = iret + nf90_def_var(ncid, 'omega_earth'  , NF90_FLOAT, (/zID,latID,lwID,timeID/), var_ssa_lw_ID) 
    iret = iret + nf90_def_var(ncid, 'g_earth'      , NF90_FLOAT, (/zID,latID,lwID,timeID/), var_asy_lw_ID) 
    !iret = iret + nf90_def_var(ncid, 'reff'         , NF90_FLOAT, (/zID,latID,timeID/), var_reff_ID)
    iret = iret + nf90_def_var(ncid, 'aod_sun'      , NF90_FLOAT, (/wlID,latID,timeID/),     var_aod_sw_ID)
    iret = iret + nf90_def_var(ncid, 'aod_earth'    , NF90_FLOAT, (/lwID,latID,timeID/),     var_aod_lw_ID)
    IF (iret /= 17*NF90_NOERR) STOP 'Error in creating file variables'

    iret = NF90_NOERR
    iret = iret + nf90_put_att(ncid, var_t_ID     , "long_name", "month")
    iret = iret + nf90_put_att(ncid, var_t_ID     , "units"    , "months since 1-1850")
    iret = iret + nf90_put_att(ncid, var_t_ID     , "calendar" , "360_day")
    iret = iret + nf90_put_att(ncid, var_dy_ID    , "long_name", "decimal year")
    iret = iret + nf90_put_att(ncid, var_dy_ID    , "units"    , "years")
    iret = iret + nf90_put_att(ncid, var_z_ID     , "units"    , "km")
    iret = iret + nf90_put_att(ncid, var_lat_ID   , "long_name", "latitude")
    iret = iret + nf90_put_att(ncid, var_lat_ID   , "units"    , "degrees_north")
!    iret = iret + nf90_put_att(ncid, var_sw_ID    , "long_name", "wavelength")
!    iret = iret + nf90_put_att(ncid, var_sw_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_wl1_sw_ID    , "long_name", "wavelength")
    iret = iret + nf90_put_att(ncid, var_wl1_sw_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_wl2_sw_ID    , "long_name", "wavelength")
    iret = iret + nf90_put_att(ncid, var_wl2_sw_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_wl1_lw_ID    , "long_name", "wavelength")
    iret = iret + nf90_put_att(ncid, var_wl1_lw_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_wl2_lw_ID    , "long_name", "wavelength")
    iret = iret + nf90_put_att(ncid, var_wl2_lw_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_ext_sw_ID   , "long_name", "aerosol extinction")
    iret = iret + nf90_put_att(ncid, var_ext_sw_ID   , "units"    , "km**-1")
    iret = iret + nf90_put_att(ncid, var_ssa_sw_ID   , "long_name", "single scattering albedo")
    iret = iret + nf90_put_att(ncid, var_ssa_sw_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_asy_sw_ID   , "long_name", "aerosol scattering asymmtery factor")
    iret = iret + nf90_put_att(ncid, var_asy_sw_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_ext_lw_ID   , "long_name", "aerosol extinction")
    iret = iret + nf90_put_att(ncid, var_ext_lw_ID   , "units"    , "km**-1")
    iret = iret + nf90_put_att(ncid, var_ssa_lw_ID   , "long_name", "single scattering albedo")
    iret = iret + nf90_put_att(ncid, var_ssa_lw_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_asy_lw_ID   , "long_name", "aerosol scattering asymmtery factor")
    iret = iret + nf90_put_att(ncid, var_asy_lw_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_aod_sw_ID  , "long_name", "aerosol optical depth")
    iret = iret + nf90_put_att(ncid, var_aod_sw_ID  , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_aod_lw_ID  , "long_name", "aerosol optical depth")
    iret = iret + nf90_put_att(ncid, var_aod_lw_ID  , "units"    , "unitless")
    !iret = iret + nf90_put_att(ncid, var_reff_ID  , "long_name", "aerosol effective radius")
    !iret = iret + nf90_put_att(ncid, var_reff_ID  , "units"    , "mu m")
    iret = iret + nf90_enddef(ncid)
    IF (iret /= 32*NF90_NOERR) STOP 'Error in creating file variable attributes'
    !
    iret = NF90_NOERR  
    iret = iret + nf90_put_var(ncid, var_t_ID      , values=mons_since_1850)
    iret = iret + nf90_put_var(ncid, var_dy_ID     , values=fyear)
    iret = iret + nf90_put_var(ncid, var_z_ID      , values=z)
    iret = iret + nf90_put_var(ncid, var_lat_ID    , values=lat)
!    iret = iret + nf90_put_var(ncid, var_sw_ID     , values=lambda)
    iret = iret + nf90_put_var(ncid, var_wl1_sw_ID , values=wl1_sw)
    iret = iret + nf90_put_var(ncid, var_wl2_sw_ID , values=wl2_sw)
    iret = iret + nf90_put_var(ncid, var_wl1_lw_ID , values=wl1_lw)
    iret = iret + nf90_put_var(ncid, var_wl2_lw_ID , values=wl2_lw)
    iret = iret + nf90_put_var(ncid, var_ext_sw_ID    , values=ext_sw)
    iret = iret + nf90_put_var(ncid, var_ssa_sw_ID    , values=ssa_sw)
    iret = iret + nf90_put_var(ncid, var_asy_sw_ID    , values=asy_sw)
    iret = iret + nf90_put_var(ncid, var_ext_lw_ID    , values=ext_lw)
    iret = iret + nf90_put_var(ncid, var_ssa_lw_ID    , values=ssa_lw)
    iret = iret + nf90_put_var(ncid, var_asy_lw_ID    , values=asy_lw)
    iret = iret + nf90_put_var(ncid, var_aod_sw_ID  , values=aod_sw)
    iret = iret + nf90_put_var(ncid, var_aod_lw_ID  , values=aod_lw)
    !iret = iret + nf90_put_var(ncid, var_reff_ID   , values=reff)
    iret = iret + nf90_close(ncid)
    IF (iret /= 17*NF90_NOERR) STOP 'error writing data or in closing file'
  end do

END PROGRAM eva_build_forcing_file_cmip6 

