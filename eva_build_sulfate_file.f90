!! @brief EVA sulfate file builder 
!!
!! @author Matthew Toohey
!!
!! @par Copyright
!! 

PROGRAM eva_build_sulfate_file 

  USE mo_eva
  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: &
       start_year = 1960, &
       end_year   = 2015

  INTEGER, PARAMETER :: nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)

  INTEGER ::          &
       nplume       , &
       nyear        , &
       ntime        , & !< number of time steps
       iret         , & !< netCDF reading return variable
       ncid         , & !< netCDF file ID
       VarID        , & !< pointer to generic dimension in netCDF file
       plumeID      , &
       timeID       , &
       var_t_ID     , &
       var_so4_ID   , &
       iyear        , & !< index for year loop
       imonth       , &     !< index for month loop
       i             

  REAL, ALLOCATABLE ::       &
       year(:)             , &
       month(:)            , &
       fyear(:)            , &
       SO4(:,:)              !

  nplume = 3
  nyear=end_year-start_year+1
  ntime=nyear*12

  ALLOCATE(fyear(ntime))
  ALLOCATE(year(ntime))
  ALLOCATE(month(ntime))
  ALLOCATE(SO4(3,ntime))

  write(*,*) 'Building sulfate timeseries' 

  i=1

  do iyear=1,nyear
    do imonth=1,nmon
      year(i)  = start_year+iyear-1
      month(i) = mons(imonth)
      i=i+1
    end do
  end do

  CALL sulfate3_timeseries(year,month,ntime,SO4)

  write(*,*) 'Writing sulfate file'

  fyear=year+(month-1)/12.

  iret = NF90_NOERR
  iret = iret + nf90_create("./eva_sulfate_timeseries.nc", NF90_CLOBBER, ncid)
  iret = iret + nf90_def_dim(ncid, 'time'  ,NF90_UNLIMITED , timeID)  
  iret = iret + nf90_def_dim(ncid, 'plume'   ,nplume , plumeID)
  IF (iret /= 3*NF90_NOERR) STOP 'Error in Creating File Dimensions'
  !
  iret = NF90_NOERR
  iret = iret + nf90_def_var(ncid, 'time'         , NF90_FLOAT, timeID, var_t_ID)
  iret = iret + nf90_def_var(ncid, 'SO4'          , NF90_FLOAT, (/plumeID,timeID/), var_so4_ID)
  iret = iret + nf90_put_att(ncid, var_t_ID       , "long_name", "fractional year")
  iret = iret + nf90_put_att(ncid, var_t_ID       , "units"    , "year")
  iret = iret + nf90_put_att(ncid, var_so4_ID     , "long_name", "sulfate mass")
  iret = iret + nf90_put_att(ncid, var_so4_ID     , "units"    , "Tg S")
  iret = iret + nf90_enddef(ncid)
  IF (iret /= 7*NF90_NOERR) STOP 'Error in creating file variables'
  !
  iret = NF90_NOERR  
  iret = iret + nf90_put_var(ncid, var_t_ID     , values=fyear)
  iret = iret + nf90_put_var(ncid, var_so4_ID   , values=SO4)
  iret = iret + nf90_close(ncid)
  IF (iret /= 3*NF90_NOERR) STOP 'error writing data or in closing file'


END PROGRAM eva_build_sulfate_file 

