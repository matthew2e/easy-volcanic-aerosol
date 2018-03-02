!>
!! @par Copyright
!!
!! @brief Program to create volcanic aerosol forcing files with the EVA module 
!!
!! @author Matthew Toohey, MPI-M, Hamburg (2016-04-08)
!!
!! $ID: n/a$
!!
!!
!! @par Copyright
!! 
!
PROGRAM eva_build_aod_file 

  USE mo_eva
  USE netcdf

  IMPLICIT NONE

  INTEGER :: &
       aod_start_year = 1990, &
       aod_end_year   = 1995

  CHARACTER(len=50) :: grid_filename         = "eva_gridfile_echam_T63_sw.nc"
  CHARACTER(len=50) :: sulfate_filename      = "eva_sulfate_timeseries.nc"
  CHARACTER(len=50) :: aod_output_dir        = "."
  CHARACTER(len=50) :: aod_file_savename     = "eva_aod_echam_T63_sw"
  CHARACTER(len=5)  :: start_yearstr, end_yearstr

  INTEGER, PARAMETER :: nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)

  INTEGER ::          &
       nplume       , & ! number of plume
       nSO4         , & ! number of timesteps in SO4 file
       nyear        , & ! number of years to save forcing files
       ntime        , & ! number of time steps (nyear*12)
       nlat         , & ! number of latitudes
       nz           , &
       iret         , & ! netCDF reading return variable
       ncid         , & ! netCDF file ID
       VarID        , & ! pointer to generic dimension in netCDF file
       timeID       , & ! pointer to time dimension in netCDF file
       latID        , & ! pointer to latitude dimension in netCDF file
       var_t_ID     , & ! pointer to time variable in netCDF file
       var_z_ID     , & ! pointer to z variable in netCDF file
       var_lat_ID   , & ! pointer to lat variable in netCDF file
       var_aod_ID   , & ! pointer to extinction variable in netCDF file
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
       aod550(:,:)         , &
       reff  (:,:)         , & ! effective radius
       ext550_vec(:)       , & ! dummy variable for call to EVA routine
       reff_vec(:)             ! dummy variable for call to EVA routine
  REAL ::                    &
       z(46)               , &
       so4_in(3)               ! sulfate triple to use in time loop

  character(8)  :: date
  character(10) :: time

  ! Input parameters from namelist
  NAMELIST /AOD_OUTPUT_FILES/ aod_output_dir, aod_file_savename, grid_filename, aod_start_year, aod_end_year
  OPEN (UNIT=30, FILE='eva_namelist_holo', STATUS='OLD')
  READ (30, NML=AOD_OUTPUT_FILES)
  CLOSE (30)

  ! Define time grid

  nyear=aod_end_year-aod_start_year+1
  ntime=nyear*12

  ALLOCATE(fyear(ntime))
  ALLOCATE(year(ntime))
  ALLOCATE(month(ntime))

  i=1

  do iyear=1,nyear
    do imonth=1,nmon
      year(i)  = aod_start_year+iyear-1
      month(i) = mons(imonth)
      i=i+1
    end do
  end do

  fyear=year+(month-1)/12.

  ! Define EVA height grid
  z = (/ (i, i = 5, 50, 1) /)
  nz = size(z)

  ! Read grid file for lat and wavelengths
  iret = nf90_open(TRIM(grid_filename), NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening grid file'
  iret = nf90_inq_dimid(ncid, "lat", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlat)

  allocate(aod550(nlat,ntime))
  allocate(reff(nlat,ntime))
  allocate(lat(nlat))
  allocate(ext550_vec(nz))
  allocate(reff_vec(nz))

  iret = nf90_inq_varid(ncid, "lat", VarID)
  iret = nf90_get_var(ncid, VarID, lat(:)  , start=(/1/) ,count=(/nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading gridfile latitude'

  ! Input sulfate data

  iret = nf90_open(TRIM(sulfate_filename), NF90_NOWRITE, ncid)
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
      aod550(j,i) = sum(ext550_vec)
      reff(j,i)  = reff_vec(1)
    end do
  end do

  ! save data in netcdf file

  write ( start_yearstr , '(i0)' ) aod_start_year
  write ( end_yearstr , '(i0)' ) aod_end_year

  write(*,*) end_yearstr, trim(end_yearstr)

  
  write(*,*) 'Writing forcing files'

  call date_and_time(DATE=date,TIME=time)

  iret = NF90_NOERR
  iret = iret + nf90_create(TRIM(aod_output_dir)//'/'//TRIM(aod_file_savename)//'_'//trim(start_yearstr)//'_'//trim(end_yearstr)//'.nc', NF90_CLOBBER, ncid)
  iret = iret + nf90_def_dim(ncid, 'time' ,NF90_UNLIMITED    , timeID)
  iret = iret + nf90_def_dim(ncid, 'lat'  ,nlat  , latID)
  IF (iret /= 4*NF90_NOERR) STOP 'Error in Creating File Dimensions'

   iret = NF90_NOERR
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,"title","EVA v1.0: stratospheric AOD")
    iret = iret + nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date(7:8)//'.'//date(5:6)//'.'//date(1:4)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:6))
    IF (iret /= 2*NF90_NOERR) STOP 'Error in Creating File Attributes' 
 
  iret = NF90_NOERR
  iret = iret + nf90_def_var(ncid, 'time'   , NF90_FLOAT, timeID,  var_t_ID)
  iret = iret + nf90_def_var(ncid, 'lat'    , NF90_FLOAT, latID,   var_lat_ID)
  iret = iret + nf90_def_var(ncid, 'aod550'    , NF90_FLOAT, (/latID,timeID/), var_aod_ID)
  !print *, trim(NF90_STRERROR(iret)) 
  iret = iret + nf90_def_var(ncid, 'reff'   , NF90_FLOAT, (/latID,timeID/), var_reff_ID)
  IF (iret /= 5*NF90_NOERR) STOP 'Error in creating file variables'

  iret = NF90_NOERR
  iret = iret + nf90_put_att(ncid, var_t_ID     , "long_name", "fractional year")
  iret = iret + nf90_put_att(ncid, var_t_ID     , "units"    , "year")
  iret = iret + nf90_put_att(ncid, var_lat_ID   , "long_name", "latitude")
  iret = iret + nf90_put_att(ncid, var_lat_ID   , "units"    , "degrees north")
  iret = iret + nf90_put_att(ncid, var_aod_ID   , "long_name", "aerosol optical depth")
  iret = iret + nf90_put_att(ncid, var_aod_ID   , "units"    , "")
  iret = iret + nf90_put_att(ncid, var_reff_ID  , "long_name", "aerosol effective radius")
  iret = iret + nf90_put_att(ncid, var_reff_ID  , "units"    , "um")
  iret = iret + nf90_enddef(ncid)

write(*,*) iret, NF90_NOERR
  IF (iret /= 10*NF90_NOERR) STOP 'Error in creating file variable attributes'
  !
  iret = NF90_NOERR  
  iret = iret + nf90_put_var(ncid, var_t_ID    , values=fyear)
  iret = iret + nf90_put_var(ncid, var_lat_ID  , values=lat)
  iret = iret + nf90_put_var(ncid, var_aod_ID  , values=aod550)
  iret = iret + nf90_put_var(ncid, var_reff_ID , values=reff)
  iret = iret + nf90_close(ncid)
  IF (iret /= 6*NF90_NOERR) STOP 'error writing data or in closing file'

END PROGRAM eva_build_aod_file 

