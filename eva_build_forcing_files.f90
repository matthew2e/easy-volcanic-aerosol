!! @brief Program to create volcanic aerosol forcing files with the EVA module 
!!
!! @author Matthew Toohey, MPI-M, Hamburg, GEOMAR, Kiel (2016-11-03)
!!
!! @par Copyright
!! 

PROGRAM eva_build_forcing_file 

  USE mo_eva
  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: &
       start_year = 1990, &
       end_year   = 1995

  CHARACTER(len=*), PARAMETER :: grid_filename         = "eva_gridfile_echam_T63_sw.nc"
  CHARACTER(len=*), PARAMETER :: sulfate_filename      = "eva_sulfate_timeseries.nc"
  CHARACTER(len=*), PARAMETER :: forcing_file_savename = "eva_forcing_echam_T63_sw"
  CHARACTER(len=4) :: yearstr

  INTEGER, PARAMETER :: nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)

  INTEGER ::          &
       nplume       , & ! number of plume
       nSO4         , & ! number of timesteps in SO4 file
       nyear        , & ! number of years to save forcing files
       ntime        , & ! number of time steps (nyear*12)
       nz           , & ! number of vertical levels
       nlat         , & ! number of latitudes
       nwl          , & ! number of wavelengths
       iret         , & ! netCDF reading return variable
       ncid         , & ! netCDF file ID
       VarID        , & ! pointer to generic dimension in netCDF file
       timeID       , & ! pointer to time dimension in netCDF file
       zID          , & ! pointer to z dimension in netCDF file
       latID        , & ! pointer to latitude dimension in netCDF file
       wlID         , & ! pointer to wavelength dimension in netCDF file
       var_t_ID     , & ! pointer to time variable in netCDF file
       var_z_ID     , & ! pointer to z variable in netCDF file
       var_lat_ID   , & ! pointer to lat variable in netCDF file
       var_wl_ID    , & ! pointer to wavelenth variable in netCDF file
       var_ext_ID   , & ! pointer to extinction variable in netCDF file
       var_ssa_ID   , & ! pointer to SSA variable in netCDF file
       var_asy_ID   , & ! pointer to ASY variable in netCDF file
       var_reff_ID  , & ! pointer to reff variable in netCDF file
       SO4_i0       , & ! index of time in SO4 file corresponding to first time of forcing file
       ind(12)      , & ! indices referencing years in full timeseries
       this_year    , & ! a single year, used in num2str 
       iyear        , & ! an index
       imonth       , & ! 
       i            , & ! an index
       j                ! another index

  !LOGICAL :: save_yearly = .true.   ! Set to false to save all data in one file         

  REAL, ALLOCATABLE ::       &
       year(:)             , & ! years
       month(:)            , & ! months
       fyear(:)            , & ! fractional year
       SO4(:,:)            , & ! sulfate data, loaded from file
       SO4_time(:)         , & ! time of sulfate data
       lat(:)              , & ! latitude of output
       lambda(:)           , & ! wavelengths of output
       ext550(:,:,:)       , & ! extinction at 550 nm
       reff  (:,:,:)       , & ! effective radius
       ext(:,:,:,:)        , & ! extinction(time,z,lat,wl)
       ssa(:,:,:,:)        , & ! ssa(time,z,lat,wl)
       asy(:,:,:,:)        , & ! asy(time,z,lat,wl)
       ext550_vec(:)       , & ! dummy variable for call to EVA routine
       reff_vec(:)         , & ! dummy variable for call to EVA routine
       ext_vec(:,:)        , & ! dummy variable for call to EVA routine
       ssa_vec(:,:)        , & ! dummy variable for call to EVA routine
       asy_vec(:,:)            ! dummy variable for call to EVA routine
  REAL ::                    &
       z(46)               , & ! height grid for EVA calculations
       so4_in(3)               ! sulfate triple to use in time loop

  character(8)  :: date
  character(10) :: time

  ! Define time grid

  write(*,*) 'Building forcing'

  nyear=end_year-start_year+1
  ntime=nyear*12

  ALLOCATE(fyear(ntime))
  ALLOCATE(year(ntime))
  ALLOCATE(month(ntime))

  i=1

  do iyear=1,nyear
    do imonth=1,nmon
      year(i)  = start_year+iyear-1
      month(i) = mons(imonth)
      i=i+1
    end do
  end do

  fyear=year+(month-1)/12.

  ! Define EVA height grid
  z = (/ (i, i = 5, 50, 1) /)
  nz = size(z)

  ! Read grid file for lat and wavelengths
  iret = nf90_open(grid_filename, NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening grid file'
  iret = nf90_inq_dimid(ncid, "lat", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlat)
  iret = nf90_inq_dimid(ncid, "nwl", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nwl)

  allocate(ext550(nz,nlat,ntime))
  allocate(ext(nz,nlat,nwl,ntime))
  allocate(ssa(nz,nlat,nwl,ntime))
  allocate(asy(nz,nlat,nwl,ntime))
  allocate(reff(nz,nlat,ntime))
  allocate(lambda(nwl))
  allocate(lat(nlat))
  allocate(ext550_vec(nz))
  allocate(reff_vec(nz))
  allocate(ext_vec(nz,nwl))
  allocate(ssa_vec(nz,nwl))
  allocate(asy_vec(nz,nwl))

  iret = nf90_inq_varid(ncid, "lat", VarID)
  iret = nf90_get_var(ncid, VarID, lat(:)  , start=(/1/) ,count=(/nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile latitude'

  iret = nf90_inq_varid(ncid, "wl_mid", VarID)
  iret = nf90_get_var(ncid, VarID, lambda(:)  , start=(/1/) ,count=(/nwl/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile wavelengths'

  ! Input sulfate data

  iret = nf90_open(sulfate_filename, NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening sulfate file'
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

  ! find index on SO4_time for t0

  SO4_i0=1 
  do while ( SO4_time(SO4_i0) < fyear(1) )
    SO4_i0 = SO4_i0 + 1
  end do

  ! for each timestep
  do i=1,ntime
    so4_in=SO4(:,SO4_i0+i-1)
    do j=1,nlat
      call eva_ext_reff(lat,z,lat(j),so4_in,nlat,nz,ext550_vec,reff_vec)
      call eva_aop_profile(lat,z,lat(j),so4_in,lambda,nwl,nlat,nz,ext_vec,ssa_vec,asy_vec)    
      ext550(:,j,i)= ext550_vec
      reff(:,j,i)  = reff_vec
      ext(:,j,:,i) = ext_vec
      ssa(:,j,:,i) = ssa_vec
      asy(:,j,:,i) = asy_vec
    end do
  end do

  write(*,*) 'Writing forcing files'
 
  call date_and_time(DATE=date,TIME=time)
  
  ! save data in netcdf file
  do i=1,nyear
    ind=(i-1)*12+(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
    this_year=year(1+(i-1)*12)
    write ( yearstr , '(i4)' ) this_year

    iret = NF90_NOERR
    iret = iret + nf90_create(forcing_file_savename//'_'//yearstr//'.nc', NF90_CLOBBER, ncid)
    iret = iret + nf90_def_dim(ncid, 'time' ,NF90_UNLIMITED    , timeID)
    iret = iret + nf90_def_dim(ncid, 'z'    ,nz    , zID)
    iret = iret + nf90_def_dim(ncid, 'lat'  ,nlat  , latID)
    iret = iret + nf90_def_dim(ncid, 'wl'   ,nwl   , wlID)
    IF (iret /= 6*NF90_NOERR) STOP 'Error in Creating File Dimensions'
    !
    iret = NF90_NOERR
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"title","EVA v1.0: stratospheric aerosol optical properties")
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date(7:8)//'.'//date(5:6)//'.'//date(1:4)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:6))
    IF (iret /= 2*NF90_NOERR) STOP 'Error in Creating File Attributes'

    iret = NF90_NOERR
    iret = iret + nf90_def_var(ncid, 'time'   , NF90_FLOAT, timeID,  var_t_ID)
    iret = iret + nf90_def_var(ncid, 'z'      , NF90_FLOAT, zID,     var_z_ID)
    iret = iret + nf90_def_var(ncid, 'lat'    , NF90_FLOAT, latID,   var_lat_ID)
    iret = iret + nf90_def_var(ncid, 'wl'     , NF90_FLOAT, wlID,    var_wl_ID)
    iret = iret + nf90_def_var(ncid, 'ext'    , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_ext_ID)
    !print *, trim(NF90_STRERROR(iret)) 
    iret = iret + nf90_def_var(ncid, 'ssa'    , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_ssa_ID)
    iret = iret + nf90_def_var(ncid, 'asy'    , NF90_FLOAT, (/zID,latID,wlID,timeID/), var_asy_ID)
    iret = iret + nf90_def_var(ncid, 'reff'   , NF90_FLOAT, (/zID,latID,timeID/), var_reff_ID)
    IF (iret /= 9*NF90_NOERR) STOP 'Error in creating file variables'

    iret = NF90_NOERR
    iret = iret + nf90_put_att(ncid, var_t_ID     , "long_name", "fractional year")
    iret = iret + nf90_put_att(ncid, var_t_ID     , "units"    , "year")
    iret = iret + nf90_put_att(ncid, var_z_ID     , "long_name", "altitude")
    iret = iret + nf90_put_att(ncid, var_z_ID     , "units"    , "km")
    iret = iret + nf90_put_att(ncid, var_lat_ID   , "long_name", "latitude")
    iret = iret + nf90_put_att(ncid, var_lat_ID   , "units"    , "degrees north")
    iret = iret + nf90_put_att(ncid, var_wl_ID    , "long_name", "wavelength")
    iret = iret + nf90_put_att(ncid, var_wl_ID    , "units"    , "mu m")
    iret = iret + nf90_put_att(ncid, var_ext_ID   , "long_name", "aerosol extinction")
    iret = iret + nf90_put_att(ncid, var_ext_ID   , "units"    , "km**-1")
    iret = iret + nf90_put_att(ncid, var_ssa_ID   , "long_name", "single scattering albedo")
    iret = iret + nf90_put_att(ncid, var_ssa_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_asy_ID   , "long_name", "aerosol scattering asymmtery factor")
    iret = iret + nf90_put_att(ncid, var_asy_ID   , "units"    , "unitless")
    iret = iret + nf90_put_att(ncid, var_reff_ID  , "long_name", "aerosol effective radius")
    iret = iret + nf90_put_att(ncid, var_reff_ID  , "units"    , "mu m")
    iret = iret + nf90_enddef(ncid)
    IF (iret /= 18*NF90_NOERR) STOP 'Error in creating file variable attributes'
    !
    iret = NF90_NOERR  
    iret = iret + nf90_put_var(ncid, var_t_ID    , values=fyear(ind))
    iret = iret + nf90_put_var(ncid, var_z_ID    , values=z)
    iret = iret + nf90_put_var(ncid, var_lat_ID  , values=lat)
    iret = iret + nf90_put_var(ncid, var_wl_ID   , values=lambda)
    iret = iret + nf90_put_var(ncid, var_ext_ID  , values=ext(:,:,:,ind))
    iret = iret + nf90_put_var(ncid, var_ssa_ID  , values=ssa(:,:,:,ind))
    iret = iret + nf90_put_var(ncid, var_asy_ID  , values=asy(:,:,:,ind))
    iret = iret + nf90_put_var(ncid, var_reff_ID , values=reff(:,:,ind))
    iret = iret + nf90_close(ncid)
    IF (iret /= 10*NF90_NOERR) STOP 'error writing data or in closing file'
  end do

END PROGRAM eva_build_forcing_file 

