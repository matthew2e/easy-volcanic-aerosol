!>
!!
!! @brief Module mo_EVA: provides volcanic aerosol optical  
!properties as a function of lat,  height, time, and wavelength
!!
!! @remarks
!!
!! @author Matthew Toohey 
!!
!! $ID: n/a$
!!
!! @par Origin
!!
!! @par Copyright
!!

module mo_EVA

   USE netcdf

   implicit none     

   real, parameter :: pi = 3.1415927
   real, parameter :: deg2rad = pi/180.

   ! --------------------------------------------------------------------------------------
   ! Sulfate box model parameters

   !----------------------------------------------------------------------------------------
   ! Input files

   CHARACTER(len=100)    :: eruption_list_filename = 'Eruption_list_GMD_1960_2015.nc'
   CHARACTER(len=50)    :: parameter_set_filename = 'EVAv1_parameter_set_v1.1.nc'
   CHARACTER(len=50)    :: Lookuptable_filename   = 'eva_Mie_lookuptables.nc'

   ! --------------------------------------------------------------------------------------
   ! Sulfate box model parameters

   REAL, SAVE ::                            &
       background_SO2_in,                   &  ! background SO2 injection
       tau_loss,                            &  ! timescale of loss in extratropical plumes
       tau_prod,                            &  ! timescale of effective SO2->AOD production in ET
       tau_mix,                             &  ! timescale of mixing between EQ and ET plumes
       tau_res,                             &  ! timescale of residual mass circ from EQ-to-ET
       seasonal_amp,                        &  ! amplitude of seasonal cycle 
       tropical_edge                           ! edge of tropical pipe

   ! ----------------------------------------------------------------------------------------
   ! Spatial shape function parameters

   REAL, SAVE ::                            &
       v_width_ET,                          &  ! vertical width of the extratropical plumes
       v_width_EQ,                          &  ! vertical width of the equatorial plume
       v_offset_EQ,                         &  ! vertical offset of the equatorial plume
       h_center_ET,                         &  ! center latitude of the extratropical plumes
       h_width_ET,                          &  ! horizontal width of the extratropical plumes
       h_width_EQ                              ! horizontal width of the equatorial plum
   REAL, SAVE, allocatable ::               &
       cline(:),                            &  ! describes height of max ext wrt latitude, length nlat
       shape_EQ(:,:),                       &  ! 2D equatorial shape function
       shape_NH(:,:),                       &  ! 2D NH shape function
       shape_SH(:,:),                       &  ! 2D SH shape function
       v_shape_EQ(:,:),                     &  ! 2D equatorial vertical shape
       v_shape_ET(:,:),                     &  ! 2D extratropical vertical shape
       h_shape_EQ(:),                       &  ! equatorial horizontal shape
       h_shape_NH(:),                       &  ! NH horizontal shape
       h_shape_SH(:)                           ! SH horizontal shape
   ! ----------------------------------------------------------------------------------------
   ! Sulfate to AOD scaling parameters

   REAL, SAVE  ::                           &
       so4_to_aod_linear,                   &  ! Scaling factor for small eruptions
       so4_to_aod_nonlin,                   &  ! coefficient for 2/3 power realtionship for larger eruptions 
       so4_to_reff,                         &  ! Sacling factor for effective radius
       nonlin_threshold                        ! threshold for nonlinear aod to so4 scaling

   ! -----------------------------------------------------------------------------------------
   ! Mie look up tables (lut) parameters 

   REAL, ALLOCATABLE ::                     & 
       lut_extrat(:,:),                     &  ! extinction relative to ext at 550nm
       lut_ssa   (:,:),                     &  ! single scattering albedo
       lut_asy   (:,:),                     &  ! asymmetry factor
       lut_reff  (:),                       &  ! effective radius
       lut_lambda(:)                           ! wavelength

   INTEGER ::                               &
       nlut_reff,                           &  ! number of effective radii in lut
       nlut_wl                                 ! number of wavelengths in lut 

   ! NAMELISTS
   NAMELIST /EVA_INPUT/ eruption_list_filename, parameter_set_filename, Lookuptable_filename

   ! ------------------------------------------------------------------------------------------
   ! PROGRAM control parameters

   LOGICAL, SAVE     :: eva_initialized = .FALSE.

   ! -------------------------------------------------------------------------------------------

contains

   subroutine init_eva(lat, z, nlat, nz)
      real, intent(in) :: lat(nlat), z(nz)
      integer, intent(in) :: nlat, nz 
      write(*,*) 'Initializing EVA'

      call read_parameter_set
      call read_aero_volc_tables
      call def_vert_centerline(lat,nlat) 

      OPEN (UNIT=10, FILE='eva_namelist_holo', STATUS='OLD')
      READ (10, NML=EVA_INPUT)
      CLOSE (10)

      write(*,*) 'Center line defined'
      call eva_shape(lat,z)
      write(*,*) 'Shape defined'
      eva_initialized = .TRUE.

   end subroutine init_eva

   subroutine eva_aop_profile(lat,z,lat0,SO4,wl,nw,nlat,nz,ext,ssa,asy)

      real, intent(in), dimension(nlat) :: lat
      real, intent(in), dimension(nz) :: z
      real, intent(in) :: lat0
      real, intent(in) :: SO4(3)
      real, intent(in) :: wl(nw)
      integer, intent(in) :: nlat, nz, nw
      real, intent(out), dimension(nz,nw) :: ext, ssa, asy

      real, dimension(nz) :: ext550, reff  

      IF (.NOT.eva_initialized) CALL init_eva(lat, z, nlat, nz)

      call get_ext550(lat, z, lat0, SO4, nlat, nz, ext550)
      call get_reff(SO4,lat,lat0,nlat,nz,reff)
      call convert_aero_props_interp(ext550,reff,wl,nw,nz,ext,ssa,asy) 

   end subroutine eva_aop_profile    

   subroutine eva_ext_reff(lat,z,lat0,SO4,nlat,nz,ext550,reff)

      real, intent(in), dimension(nlat) :: lat
      real, intent(in), dimension(nz) :: z
      real, intent(in) :: lat0
      real, intent(in) :: SO4(3)
      integer, intent(in) :: nlat, nz
      real, intent(out), dimension(nz) :: ext550, reff

      IF (.NOT.eva_initialized) CALL init_eva(lat, z, nlat, nz)

      call get_ext550(lat, z, lat0, SO4, nlat, nz, ext550)
      call get_reff(SO4,lat,lat0,nlat,nz,reff)

   end subroutine eva_ext_reff

   subroutine eva_aod(lat,lat0,SO4,wl,nlat,nw,nz,aod)
      ! To do, get the AOD at lat0 given SO4 amount
      ! first get aod550, then convert aod using look up table
      real, intent(in), dimension(nlat) :: lat
      real, intent(in) :: lat0
      real, intent(in) :: SO4(3)
      real, intent(in) :: wl(nw)
      integer, intent(in) :: nlat, nw, nz
      real, intent(out), dimension(nw) :: aod
      ! local variables
      real :: aod550, reff(nz)

      call get_aod550(lat,lat0,SO4,nlat,aod550)
      call get_reff(SO4,lat,lat0,nlat,nz,reff)
      call convert_aod_interp(aod550,reff(1),wl,nw,aod)

   end subroutine eva_aod

   subroutine get_aod550(lat,lat0,SO4,nlat,aod550)

     ! get aod550 at given latitude lat0

     real, intent(in), dimension(nlat) :: lat
     real, intent(in) :: lat0
     real, intent(in) :: SO4(3)
     integer, intent(in) :: nlat
     real, intent(out) :: aod550

     ! local variables
     integer :: lat0_ind, i
     real :: sulfate_aerosol_mass(3)
     real :: SO4_lat

     nonlin_threshold = (so4_to_aod_nonlin/so4_to_aod_linear)**(3.)

     !write(*,*) z
      ! find index of lat0 in lat
      do i=1,nlat
         if ( lat(i) .eq. lat0 ) then
            lat0_ind=i
          end if
      end do

      SO4_lat = h_shape_EQ(lat0_ind)*SO4(2) &
                      + h_shape_SH(lat0_ind)*SO4(1) &
                      + h_shape_NH(lat0_ind)*SO4(3)

      if ( SO4_lat < nonlin_threshold ) then
         aod550=so4_to_aod_linear*SO4_lat
      else
         aod550=so4_to_aod_nonlin*(SO4_lat**(2.0/3.0))
      end if

   end subroutine get_aod550 


   subroutine get_ext550(lat,z,lat0,SO4,nlat,nz,ext550)

     ! get extinction profile at a given latitude lat0

     real, intent(in), dimension(nlat) :: lat
     real, intent(in), dimension(nz) :: z
     real, intent(in) :: lat0
     real, intent(in) :: SO4(3)
     integer, intent(in) :: nlat, nz
     real, intent(out), dimension(nz) ::ext550 

     ! local variables
     integer :: lat0_ind, i
     real :: sulfate_aerosol_mass(3)
     real :: aod, SO4_lat, v_integral
     real :: v_shape(nz), v_shape_norm(nz)

     nonlin_threshold = (so4_to_aod_nonlin/so4_to_aod_linear)**(3.)

     !write(*,*) z
      ! find index of lat0 in lat
      do i=1,nlat
         if ( lat(i) .eq. lat0 ) then
            lat0_ind=i
          end if
      end do

      SO4_lat = h_shape_EQ(lat0_ind)*SO4(2) &
                      + h_shape_SH(lat0_ind)*SO4(1) &
                      + h_shape_NH(lat0_ind)*SO4(3) 
     
      if ( SO4_lat < nonlin_threshold ) then
         aod=so4_to_aod_linear*SO4_lat
      else
         aod=so4_to_aod_nonlin*(SO4_lat**(2.0/3.0))
      end if

      v_shape=shape_EQ(lat0_ind,:)*SO4(2) &
                      + shape_SH(lat0_ind,:)*SO4(1) &
                      + shape_NH(lat0_ind,:)*SO4(3)
      ! normalize, needed for transition zones?
      v_integral=0
      do i=1,nz-1      
        v_integral=v_integral + 0.5*(v_shape(i)+v_shape(i+1))*(z(i+1)-z(i))
      end do
      if ( sum(v_shape) > 0 ) then 
         !v_shape_norm=v_shape/sum(v_shape*(z(2)-z(1)))
         v_shape_norm=v_shape/v_integral
      else
         v_shape_norm=0
         write(*,*) "normalization <0"
      end if

      ext550=aod*v_shape_norm

   end subroutine get_ext550 

   subroutine get_reff(SO4,lat,lat0,nlat,nz,reff)
      real, intent(in) :: SO4(3), lat(nlat), lat0
      integer, intent(in) :: nlat, nz
      real, intent(out) :: reff(nz)

      integer :: i, lat0_ind

      do i=1,nlat
         if ( lat(i) .eq. lat0 ) then
            lat0_ind=i
         end if
      end do

      reff = so4_to_reff*(shape_EQ(lat0_ind,:)*SO4(2) &
                      + shape_SH(lat0_ind,:)*SO4(1) &
                      + shape_NH(lat0_ind,:)*SO4(3))**(1./3.)
      reff=max(maxval(reff),0.2)

   end subroutine get_reff

   subroutine def_vert_centerline(lat,nlat)
      real, intent(in) :: lat(nlat)
      integer, intent(in) :: nlat
      integer :: i
      real :: cline_lat(37), cline_raw(37)
      INTEGER ::          &
         iret         , & !< netCDF reading return variable
         ncid         , & !< netCDF file ID
         VarID            !< pointer to generic dimension in netCDF file

      allocate(cline(nlat))

      iret = nf90_open(TRIM(parameter_set_filename), NF90_NOWRITE, ncid)
      IF (iret /= NF90_NOERR) STOP 'Error in opening parameter file'

      iret = nf90_inq_varid(ncid, "cline_lat", VarID)
      iret = nf90_get_var(ncid, VarID, cline_lat(:)  , start=(/1/) ,count=(/37/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading cline_lat'

      iret = nf90_inq_varid(ncid, "cline", VarID)
      iret = nf90_get_var(ncid, VarID, cline_raw(:)  , start=(/1/) ,count=(/37/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading cline'

      do i=1,nlat
         cline(i)=interp(cline_lat,cline_raw,lat(i),37)
      end do

   end subroutine def_vert_centerline

   subroutine EVA_shape(lat,z)

      ! input/output
      real, intent(in), dimension(:) :: lat, z

      ! local variables
      real, dimension(:), allocatable :: sinlat
      real :: h_sin_center_ET, h_sin_width_ET, h_sin_width_EQ, v_integral
      integer :: nlat, nz, i, j

      real, dimension(:,:), allocatable :: v_shape_EQ, v_shape_ET

      ! allocate arrays
      nlat = size(lat)
      nz=size(z)
      allocate(sinlat(nlat))
      allocate(h_shape_EQ(nlat))
      allocate(h_shape_NH(nlat))
      allocate(h_shape_SH(nlat))
      allocate(v_shape_EQ(nlat,nz))
      allocate(v_shape_ET(nlat,nz))
      allocate(shape_EQ(nlat,nz))
      allocate(shape_NH(nlat,nz))
      allocate(shape_SH(nlat,nz))

      ! Horizontal shape calculations

      sinlat=sin( lat * pi / 180. )
      h_sin_center_ET = sin( h_center_ET * pi / 180. )
      h_sin_width_ET = sin( h_width_ET * pi / 180. )
      h_sin_width_EQ = sin( h_width_EQ * pi / 180. )

      h_shape_EQ=exp( -(sinlat-0)**2. / (2. * h_sin_width_EQ )**2. )
      h_shape_NH=exp( -(sinlat-h_sin_center_ET)**2. / (2. * h_sin_width_ET)**2. )
      h_shape_SH=exp( -(sinlat+h_sin_center_ET)**2. / (2. * h_sin_width_ET)**2. )

      ! normalize (needed because extratropical shapes won't extend to
      ! full width. Also, normalizing by area puts shapes into units of
      ! km**-2)

      h_shape_EQ = h_shape_EQ / area_gmean(h_shape_EQ,lat)
      h_shape_NH = h_shape_NH / area_gmean(h_shape_NH,lat)
      h_shape_SH = h_shape_SH / area_gmean(h_shape_SH,lat)

      nlat=size(cline)

      ! build 3D shape

      do i=1,nlat
         v_shape_ET(i,:) = exp( -(z-cline(i))**2. / (2. * v_width_ET)**2.) &
                                  /(v_width_ET*sqrt(2.*pi))
         v_shape_EQ(i,:) = exp(-(z-cline(i)-v_offset_EQ)**2. / (2. * v_width_EQ)**2. ) &
                                  /(v_width_EQ*sqrt(2.*pi))

         ! normalize (assuming equally spaced z)
         v_integral=0
         do j=1,nz-1
           v_integral=v_integral + 0.5*(v_shape_ET(i,j)+v_shape_ET(i,j+1))*(z(j+1)-z(j))
         end do
         !v_shape_ET(i,:)=v_shape_ET(i,:)/sum(v_shape_ET(i,:)*(z(2)-z(1)))
         v_shape_ET(i,:)=v_shape_ET(i,:)/v_integral
         v_integral=0
         do j=1,nz-1
           v_integral=v_integral + 0.5*(v_shape_EQ(i,j)+v_shape_EQ(i,j+1))*(z(j+1)-z(j))
         end do
         !v_shape_EQ(i,:)=v_shape_EQ(i,:)/sum(v_shape_EQ(i,:)*(z(2)-z(1)))
         v_shape_EQ(i,:)=v_shape_EQ(i,:)/v_integral
         ! include horizontal shape
         shape_EQ(i,:)=v_shape_EQ(i,:)*h_shape_EQ(i)
         shape_NH(i,:)=v_shape_ET(i,:)*h_shape_NH(i)
         shape_SH(i,:)=v_shape_ET(i,:)*h_shape_SH(i)

      end do

      deallocate(sinlat)

   end subroutine EVA_shape

   subroutine sulfate3_timeseries(year,month,n,SO4)

   ! input/output
      real, intent(out), dimension(3,n) :: SO4
      real, intent(in), dimension(n) :: year, month
      integer, intent(in) :: n
      ! local variables
      integer :: ntime, i, j, k

      real, dimension(3,n) :: SO2_in, SO2, SO4_new
      real, dimension(:), allocatable :: erup_ssi, erup_lat, erup_hemi, erup_dur
      real :: C, tau_loss_EQ_use, tau_loss_ET_use, SO4_tot
      real, dimension(n) :: hemi_corr_nh, hemi_corr_sh
      integer, dimension(:), allocatable :: erup_region, erup_year, erup_month, erup_day
      real :: SHmix, NHmix, SHtrans, NHtrans, seas_asy
      real :: S_tot_t0, SO2_tot_t0, SO4_tot_t0, SO4_ET2EQ_ratio_t0

   INTEGER ::          &
       iret         , & !< netCDF reading return variable
       ncid         , & !< netCDF file ID
       VarID        , & !< pointer to generic dimension in netCDF file
       nerup        

      call read_parameter_set
      OPEN (UNIT=10, FILE='eva_namelist_holo', STATUS='OLD')
      READ (10, NML=EVA_INPUT)
      CLOSE (10)

      ntime=size(year)
      
      SO2_in=0
      SO2_in(1,:)=background_SO2_in/12.0/2.0
      SO2_in(3,:)=background_SO2_in/12.0/2.0
      SO2=0
      SO4=0
      SO4_new=0
      hemi_corr_nh=1.0
      hemi_corr_sh=1.0

      ! initialize SO4 based on background SO2_in (no spin-up needed)
      ! (this is very close to being right, SO4 still dips every so slighlty in
      ! first months)
      S_tot_t0=background_SO2_in/12.0*tau_loss
      SO4_tot_t0 = S_tot_t0/(1.0+tau_prod/tau_loss)
      SO2_tot_t0 = S_tot_t0/(1.0+tau_loss/tau_prod)
      SO4_ET2EQ_ratio_t0 = ((tau_mix*tau_res/2.0) + tau_mix*tau_loss + tau_loss*tau_res) / (tau_loss*tau_res)
      SO4(2,1)=SO4_tot_t0 / ( 2.0*SO4_ET2EQ_ratio_t0 + 1.0)
      SO4(1,1)=SO4_tot_t0 / ( 2.0 + 1.0/SO4_ET2EQ_ratio_t0)
      SO4(3,1)=SO4_tot_t0 / ( 2.0 + 1.0/SO4_ET2EQ_ratio_t0)
      SO2(2,1)=SO2_tot_t0 /  ( 2.0*SO4_ET2EQ_ratio_t0 + 1.0)
      SO2(1,1)=SO2_tot_t0 / ( 2.0 + 1.0/SO4_ET2EQ_ratio_t0) 
      SO2(3,1)=SO2_tot_t0 / ( 2.0 + 1.0/SO4_ET2EQ_ratio_t0)


      ! Read eruption list
      write(*,*) 'Reading file: ', TRIM(eruption_list_filename)
      iret = nf90_open(TRIM(eruption_list_filename), NF90_NOWRITE, ncid)
      IF (iret /= NF90_NOERR) STOP 'Error in opening eruption list file'
      iret = nf90_inq_dimid(ncid, "nerup", VarID)
      iret = nf90_inquire_dimension(ncid, VarID, len = nerup)

      allocate(erup_year(nerup))
      allocate(erup_month(nerup))
      allocate(erup_day(nerup))
      allocate(erup_lat(nerup))
      allocate(erup_ssi(nerup))
      allocate(erup_region(nerup))
      allocate(erup_hemi(nerup))
      allocate(erup_dur(nerup))

      iret = nf90_inq_varid(ncid, "year", VarID)
      iret = nf90_get_var(ncid, VarID, erup_year(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption year'

      iret = nf90_inq_varid(ncid, "month", VarID)
      iret = nf90_get_var(ncid, VarID, erup_month(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption month'

      iret = nf90_inq_varid(ncid, "day", VarID)
      iret = nf90_get_var(ncid, VarID, erup_day(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption day'

      iret = nf90_inq_varid(ncid, "lat", VarID)
      iret = nf90_get_var(ncid, VarID, erup_lat(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption lat'

      iret = nf90_inq_varid(ncid, "ssi", VarID)
      iret = nf90_get_var(ncid, VarID, erup_ssi(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption ssi'

      iret = nf90_inq_varid(ncid, "hemi", VarID)
      iret = nf90_get_var(ncid, VarID, erup_hemi(:)  , start=(/1/) ,count=(/nerup/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading eruption hemi'

      iret = nf90_inq_varid(ncid, "duration", VarID)
      IF (iret /= NF90_NOERR) THEN 
         erup_dur=0
      ELSE
         iret = nf90_get_var(ncid, VarID, erup_dur(:)  , start=(/1/) ,count=(/nerup/))
         IF (iret /= NF90_NOERR) STOP 'Error in reading eruption duration'
      END IF


      iret = nf90_close  (ncid)

      !write(*,*) erup_ssi
 
      do i=1,nerup
         IF ( erup_lat(i) < -tropical_edge ) THEN
            erup_region(i) = 1
         ELSE IF ( erup_lat(i) > tropical_edge ) THEN
            erup_region(i) = 3
         ELSE
            erup_region(i) = 2
         END IF
      end do

      ! Insert S injections into timeseries
      ! for each eruption in list
      do i=1,nerup
         ! if date is >= Year_start & date <= Year_end
         if (erup_year(i) .ge. year(1) .and. erup_year(i) .le. year(ntime) ) then
            ! loop through length of year/month to find matching date
            do j=1,ntime
               if (year(j) .eq. erup_year(i) .and. month(j) .eq. erup_month(i)) then
                  ! set SO2_in for this date = elist_SO2
                  if (erup_dur(i) .gt. 0) then
                     SO2_in(erup_region(i),j:j+erup_dur(i))=SO2_in(erup_region(i),j:j+erup_dur(i)) + erup_ssi(i)/erup_dur(i)
                  else
                     SO2_in(erup_region(i),j)=SO2_in(erup_region(i),j) + erup_ssi(i)
                  end if
                  seas_asy=0.4*seasonal_amp*cos(2.0*pi*(month(j)+4.0)/12.0)+1 ! a rough approximation of what the 
                                                                 ! seasonal paramterized transport produces
                  if ( erup_hemi(i) > 0.0 ) then
                     if ( abs(seas_asy - erup_hemi(i)) > 0.1 ) then
                        if ( erup_hemi(i) < 1.0 ) then
                           do k=j,min(j+17,ntime)
                              hemi_corr_nh(k)=seas_asy/erup_hemi(i) 
                           end do
                        else
                           do k=j,min(j+17,ntime)
                              hemi_corr_sh(k) = erup_hemi(i)/seas_asy
                           end do
                        end if
                     end if
                  end if
               end if
            end do
         end if
      end do
   
      ! Run box-model calculations to produce SO4 timeseries
     
      !SO2(:,1)=SO2_in(:,1)

      ! SO2 is converted into SO4
      do i=2,ntime
         SO2(1,i)=SO2(1,i-1)*exp(-1/tau_prod)*exp(-1/tau_loss) +SO2_in(1,i)
         SO2(2,i)=SO2(2,i-1)*exp(-1/tau_prod)*exp(-1/tau_loss) +SO2_in(2,i)
         SO2(3,i)=SO2(3,i-1)*exp(-1/tau_prod)*exp(-1/tau_loss) +SO2_in(3,i)
         SO4_new(1,i)=SO2(1,i-1)*(1-exp(-1/tau_prod))
         SO4_new(2,i)=SO2(2,i-1)*(1-exp(-1/tau_prod))
         SO4_new(3,i)=SO2(3,i-1)*(1-exp(-1/tau_prod))
      end do

      ! SO4 grows and mixes

      do i=2,ntime
         ! compute mixing and transport timescales with seasonal
         ! variability. Seasonality modeled as cosine, with maxima in
         ! January (NH) and July (SH).
         SHmix = (1.0/(tau_mix*hemi_corr_sh(i)))*(seasonal_amp*cos(2*pi*mod(i-1+6,12)/12)+1)
         NHmix = (1.0/(tau_mix*hemi_corr_nh(i)))*(seasonal_amp*cos(2*pi*mod(i-1,12)/12)+1)
         SHtrans = (1.0/(hemi_corr_sh(i)*tau_res))*(seasonal_amp*cos(2*pi*mod(i-1+6,12)/12)+1)
         NHtrans = (1.0/(hemi_corr_nh(i)*tau_res))*(seasonal_amp*cos(2*pi*mod(i-1,12)/12)+1)

         ! accelerate loss for large sulfur loadings (special case for
         ! 1257 eruption: relationship based roughly on ECHAM5-HAM
         ! simulations of large eruptions. But eventually this should be
         ! replace, the loss rate would be better related to Reff, this
         ! requires merging the sulfate and Reff growth routines
         SO4_tot=sum(SO4(:,i-1))
         if (SO4_tot > 10) then
            tau_loss_EQ_use=((tau_loss-6.0)/0.3679)*exp(-SO4_tot/10.0)+6.0
            tau_loss_ET_use=((tau_loss-6.0)/0.3679)*exp(-SO4_tot/10.0)+6.0
         else
            tau_loss_EQ_use=tau_loss
            tau_loss_ET_use=tau_loss
         end if


         SO4(2,i) = (SO4(2,i-1)+SO4_new(2,i))*exp(-1/tau_loss_EQ_use) &
                          ! add production and subtract STE loss
                  - (SO4(2,i-1) - SO4(3,i-1))*SHmix &
                           ! subtract loss by mixing to SH
                  - (SO4(2,i-1)- SO4(1,i-1))*NHmix &
                           ! subtract loss by mixing to NH
                   - SO4(2,i-1)*NHtrans - SO4(2,i-1)*SHtrans
                           ! subtract loss by residual circ to both hemispheres

         SO4(1,i) = (SO4(1,i-1)+SO4_new(1,i))*exp(-1/tau_loss_ET_use) &
                           ! add production and subtract STE loss
                  + (SO4(2,i-1) - SO4(1,i-1))*SHmix &
                           ! add mixing from EQ
                   + SO4(2,i-1)*SHtrans
                           ! add residual circ from EQ

         SO4(3,i) = (SO4(3,i-1)+SO4_new(3,i))*exp(-1/tau_loss_ET_use) &
                           ! add production and subtract STE loss
                  + (SO4(2,i-1) - SO4(3,i-1))*NHmix &
                           ! add mixing from EQ
                   + SO4(2,i-1)*NHtrans
                           ! add residual circ from EQ
 
      end do
      
      !write(*,*) SO4(:,3)

      deallocate(erup_ssi)
      deallocate(erup_region)
      deallocate(erup_year)
      deallocate(erup_month)

   end subroutine sulfate3_timeseries

   function area_gmean(x,lat)        
      real, dimension(:) :: x, lat          
      real, dimension(:), allocatable :: w, rlat
      real :: area_gmean, lat_bounds(2)
      integer :: i, nlat
      real, parameter :: pi = 3.1415927
      nlat=size(lat)
      allocate(w(nlat))
      allocate(rlat(nlat))
      rlat = (pi/180.)*lat
      w=cos(rlat)
      w = w / sum( w )
      area_gmean = sum ( w * x )

      deallocate(w)
      deallocate(rlat) 

   end function area_gmean

   subroutine read_parameter_set
      INTEGER ::          &
         iret         , & !< netCDF reading return variable
         ncid         , & !< netCDF file ID
         VarID            !< pointer to generic dimension in netCDF file

    write(*,*) 'Reading parameter set: ', TRIM(parameter_set_filename)     
 
      iret = nf90_open(TRIM(parameter_set_filename), NF90_NOWRITE, ncid)
      IF (iret /= NF90_NOERR) STOP 'Error in opening parameter file'
      iret = nf90_inq_varid(ncid, "version", VarID)

      iret = nf90_inq_varid(ncid, "background_SO2_in", VarID)
      iret = nf90_get_var(ncid, VarID, background_SO2_in)
      IF (iret /= NF90_NOERR) STOP 'Error in reading background_SO2_in'

      iret = nf90_inq_varid(ncid, "tau_loss", VarID)
      iret = nf90_get_var(ncid, VarID, tau_loss)
      IF (iret /= NF90_NOERR) STOP 'Error in reading tau_loss'

      iret = nf90_inq_varid(ncid, "tau_prod", VarID)
      iret = nf90_get_var(ncid, VarID, tau_prod)
      IF (iret /= NF90_NOERR) STOP 'Error in reading tau_prod'

      iret = nf90_inq_varid(ncid, "tau_mix", VarID)
      iret = nf90_get_var(ncid, VarID, tau_mix)
      IF (iret /= NF90_NOERR) STOP 'Error in reading tau_mix'

      iret = nf90_inq_varid(ncid, "tau_res", VarID)
      iret = nf90_get_var(ncid, VarID, tau_res)
      IF (iret /= NF90_NOERR) STOP 'Error in reading tau_res'

      iret = nf90_inq_varid(ncid, "seasonal_amp", VarID)
      iret = nf90_get_var(ncid, VarID, seasonal_amp)
      IF (iret /= NF90_NOERR) STOP 'Error in reading seasonal_amp'

      iret = nf90_inq_varid(ncid, "tropical_edge", VarID)
      iret = nf90_get_var(ncid, VarID, tropical_edge)
      IF (iret /= NF90_NOERR) STOP 'Error in reading tropical_edge'

      iret = nf90_inq_varid(ncid, "v_width_ET", VarID)
      iret = nf90_get_var(ncid, VarID, v_width_ET)
      IF (iret /= NF90_NOERR) STOP 'Error in reading v_width_ET'

      iret = nf90_inq_varid(ncid, "v_width_EQ", VarID)
      iret = nf90_get_var(ncid, VarID, v_width_EQ)
      IF (iret /= NF90_NOERR) STOP 'Error in reading v_width_EQ'

      iret = nf90_inq_varid(ncid, "v_offset_EQ", VarID)
      iret = nf90_get_var(ncid, VarID, v_offset_EQ)
      IF (iret /= NF90_NOERR) STOP 'Error in reading v_offset_EQ'

      iret = nf90_inq_varid(ncid, "h_center_ET", VarID)
      iret = nf90_get_var(ncid, VarID, h_center_ET )
      IF (iret /= NF90_NOERR) STOP 'Error in reading h_center_ET'

      iret = nf90_inq_varid(ncid, "h_width_ET", VarID)
      iret = nf90_get_var(ncid, VarID, h_width_ET)
      IF (iret /= NF90_NOERR) STOP 'Error in reading h_width_ET'

      iret = nf90_inq_varid(ncid, "h_width_EQ", VarID)
      iret = nf90_get_var(ncid, VarID, h_width_EQ)
      IF (iret /= NF90_NOERR) STOP 'Error in reading h_width_EQ'

      iret = nf90_inq_varid(ncid, "so4_to_aod_linear", VarID)
      iret = nf90_get_var(ncid, VarID, so4_to_aod_linear)
      IF (iret /= NF90_NOERR) STOP 'Error in reading so4_to_aod_linear'

      iret = nf90_inq_varid(ncid, "so4_to_aod_nonlin", VarID)
      iret = nf90_get_var(ncid, VarID, so4_to_aod_nonlin)
      IF (iret /= NF90_NOERR) STOP 'Error in reading so4_to_aod_nonlin'

      iret = nf90_inq_varid(ncid, "so4_to_reff", VarID)
      iret = nf90_get_var(ncid, VarID, so4_to_reff)
      IF (iret /= NF90_NOERR) STOP 'Error in reading so4_to_reff'

      !write(*,*) tau_loss, tau_prod, so4_to_reff

   end subroutine read_parameter_set


   subroutine read_aero_volc_tables
      ! read in lookup tables from netcdf file
      INTEGER ::          &
         iret         , & !< netCDF reading return variable
         ncid         , & !< netCDF file ID
         VarID        , & !< pointer to generic dimension in netCDF file
         nreff        , & !< length of effective radius vector in lookuptable
         nwl          
      write(*,*) 'Reading look up table:', TRIM(Lookuptable_filename)
      ! read in the lookup tables
      !
      iret = nf90_open(TRIM(Lookuptable_filename), NF90_NOWRITE, ncid)
      IF (iret /= NF90_NOERR) STOP 'Error in opening lookup tables file'
      iret = nf90_inq_dimid(ncid, "reff", VarID)
      iret = nf90_inquire_dimension(ncid, VarID, len = nreff)
      iret = nf90_inq_dimid(ncid, "wl", VarID)
      iret = nf90_inquire_dimension(ncid, VarID, len = nwl)

      allocate(lut_lambda(nwl))
      allocate(lut_reff(nreff))
      allocate(lut_extrat(nreff,nwl))
      allocate(lut_ssa(nreff,nwl))
      allocate(lut_asy(nreff,nwl))

      iret = nf90_inq_varid(ncid, "reff", VarID)
      iret = nf90_get_var(ncid, VarID, lut_reff(:)  , start=(/1/)  ,count=(/nreff/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading effective radii'
      iret = nf90_inq_varid(ncid, "wl", VarID)
      iret = nf90_get_var(ncid, VarID, lut_lambda(:)  , start=(/1/)  ,count=(/nwl/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading wavelengths'
      iret = nf90_inq_varid(ncid, "extrat", VarID)
      iret = nf90_get_var(ncid, VarID, lut_extrat(:,:), start=(/1,1/),count=(/nreff,nwl/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading extrat'
      iret = nf90_inq_varid(ncid, "ssa", VarID)
      iret = nf90_get_var(ncid, VarID, lut_ssa(:,:), start=(/1,1/),count=(/nreff,nwl/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading ssa'
      iret = nf90_inq_varid(ncid, "asy", VarID)
      iret = nf90_get_var(ncid, VarID, lut_asy(:,:), start=(/1,1/),count=(/nreff,nwl/))
      IF (iret /= NF90_NOERR) STOP 'Error in reading asy'

      iret = nf90_close  (ncid)

      nlut_reff=nreff
      nlut_wl=nwl
      !write(*,*) lut_lambda
   end subroutine read_aero_volc_tables

   subroutine convert_aero_props_interp(ext550, reff, wl, nw, nz, ext, ssa, asy)
      ! Conversion of ext550 and Reff to ext(wl), ssa(wl) and asy(wl)
      ! based on interpolation in wavelength to look up table values
      ! To do: we will replace the look up tables for this option with
      !  something more useful for other models, i.e., higher res

      real, intent(in) :: ext550(nz), reff(nz), wl(nw)
      integer, intent(in) :: nw, nz
      real, intent(out) :: ext(nz,nw), ssa(nz,nw), asy(nz,nw)
      integer :: i, j
      real :: reff_limit
      real :: reff_to_use(nz)

      reff_limit=1.3
      reff_to_use=min(reff,reff_limit)

      ! For each height and wavelength, find ext, ssa and asym by bilinear interpolation
      !write(*,*) lut_extrat
      do i=1,nz
         do j=1,nw
            ext(i,j) = ext550(i) * interp2(lut_reff, lut_lambda, lut_extrat, reff_to_use(i), wl(j), nlut_reff, nlut_wl) 
            ssa(i,j) = interp2(lut_reff, lut_lambda, lut_ssa, reff_to_use(i), wl(j), nlut_reff, nlut_wl) 
            asy(i,j) = interp2(lut_reff, lut_lambda, lut_asy, reff_to_use(i), wl(j), nlut_reff, nlut_wl) 
         end do
      end do

   end subroutine convert_aero_props_interp 

   subroutine convert_aod_interp(aod550, reff, wl, nw, aod)
      ! Conversion of ext550 and Reff to ext(wl), ssa(wl) and asy(wl)
      ! based on interpolation in wavelength to look up table values
      ! To do: we will replace the look up tables for this option with
      !  something more useful for other models, i.e., higher res

      real, intent(in) :: aod550, reff, wl(nw)
      integer, intent(in) :: nw 
      real, intent(out) :: aod(nw)
      integer :: j
      real :: reff_limit
      real :: reff_to_use

      reff_limit=1.3
      reff_to_use=min(reff,reff_limit)
      ! For each wavelength, find aod by bilinear
      ! interpolation
      !write(*,*) lut_extrat
      do j=1,nw
         aod(j) = aod550 * interp2(lut_reff, lut_lambda, lut_extrat, reff_to_use, wl(j), nlut_reff, nlut_wl)
      end do

   end subroutine convert_aod_interp

   function interp(xin, yin, xout, n)
      ! simple routine to linearly interpolate
      ! xin should be monotonically increasing
      ! xout should be a single value
      ! For now, for xout values outside the range of xin, 
      !  the function returns the "nearest neighbour" 

      real :: xin(n), yin(n), xin_rev(n), yin_rev(n)
      integer :: n
      real :: interp, xout, slope
      integer :: i

      !if (xin(2).gt.xin(1)) then
      !   do i=1,n
      !     xin_rev(n-i+1)=xin(i)
      !     yin_rev(n-i+1)=yin(i)
      !   end do
      !   xin=xin_rev
      !   yin=yin_rev
      !end if
      i=1
      if (xout < xin(1)) then
         interp=yin(1)
         return
      end if

      if (xout > xin(n)) then
         interp=yin(n)
         return
      end if
!
      do while (xin(i) < xout)
         i=i+1
      end do

      slope = (yin(i)-yin(i-1))/(xin(i)-xin(i-1))
      interp = yin(i-1) + slope * ( xout - xin(i-1) )

   end function interp


   function interp2(xin, yin, zin, xout, yout, nx, ny)
   ! simple routine to bilinearly interpolate
   ! xin, yin should be monotonically increasing
   ! xout, yout should be single values
   ! For now, for xout,yout values outside the range of xin, yin u
   !  the function returns 0 

      real :: xin(nx), yin(ny), zin(nx,ny), xout, yout
      real :: interp2, z_x_y1, z_x_y2
      integer :: i, j, nx, ny
      real :: xin_rev(nx), yin_rev(ny),zin_rev(nx,ny)

      !if (xin(2).gt.xin(1)) then
      !   do i=1,nx
      !      xin_rev(nx-i+1) = xin(i)
      !      zin_rev(nx-i+1,:) = zin(i,:)
      !   end do
      !   xin=xin_rev
      !   zin=zin_rev
      !end if
      !if (yin(2).gt.yin(1)) then
      !   do i=1,ny
      !      yin_rev(ny-i+1) = yin(i)
      !      zin_rev(:,ny-i+1) = zin(:,i)
      !   end do
      !   yin=yin_rev
      !   zin=zin_rev
      !end if

      if (xout < xin(1)) then
         interp2=0
         return
      end if

      if (xout > xin(nx)) then
         interp2=0
         return
      end if
      if (yout < yin(1)) then
         interp2=0
         return
      end if

      if (yout > yin(ny)) then
         interp2=0 
         return
      end if 
        
      j=1 
      do while (yin(j) < yout)
         j=j+1
      end do
    
      ! first interpolate in x at y1 and y2
      z_x_y1 = interp(xin, zin(:,j-1), xout, nx)
      z_x_y2 = interp(xin, zin(:,j), xout, nx)

      ! then interpolate in y    
      interp2 = interp(yin(j-1:j), (/ z_x_y1, z_x_y2 /), yout, 2)

   end function interp2

 
end module mo_EVA
