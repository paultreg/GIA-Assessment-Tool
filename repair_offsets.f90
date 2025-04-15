  program repair_offsets

! program to read in a co-ordinate time series, identify and repair offsets detected within the time series. The process involves:
!
! 1. using the MIDAS method (Blewitt et al., 2016) to detrend the time series for each co-ordinate
! 2. calculate a running mean, using a 12-month sliding window
! 3. calculate the correlations of this running mean with a 2-year function of expected change in the presence of an offset
! 4. identify the epochs of peaks of correlation (or anti-correlation)
! 5. adjust epochs by 1.5 years to identify the epochs of offsets
! 6. calculate the magnitude of the offsets at these epochs
! 7. output a time series of coordinates corrected for the identified offsets
!
! P. Tregoning
! 24 June 2020

  implicit none

  character*100 :: input_crd_file,output_crd_file     ! input/output coordinate file names
  real*8        :: corrl_limit                        ! threshold for identifying offset correlation limits

! variables related to midas calculations
  integer*4          :: n_obs          ! number of epochs of coordinate observations
  real*8,allocatable :: neu(:,:)       ! array for input coordinates
  real*8,allocatable :: neu_sig(:,:)   ! array for input coordinate sigmas
  real*8,allocatable :: midas_vel(:,:) ! array for midas velocity estimates
  real*8,allocatable :: epochs(:)      ! epochs of coordinate estimates
  integer*4          :: median_val(3)  ! median values of midas velocity computations for each coordinate

! counters
  integer*4 :: icomp,iobs

! variables related to smoothing computations
  real*8,allocatable :: midas_smooth(:,:)  ! smoothed estiamtes of midas velocities
  real*8             :: smooth_window      ! time interval (in decimal years) over which to smooth the midas velocity estimates

! variables related to computing the correlations and offsets
  real*8,allocatable :: corrl(:,:)         ! correlations of a synthetic triangular-shaped signal with the smoothed midas velocities
  real*8,allocatable :: offsets(:,:)       ! offsets for each component at each epoch

! other variables
  integer*4,parameter  :: luin=10,luout=11
  integer*4            :: ioerr
  character*100        :: line,arg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,input_crd_file)
! provide help if needed
  if(input_crd_file(1:1) == " ")then
    print*,'Runstring:  repair_offsets input_crd_file output_crd_file corrl_limit'
    stop
  else
    ! open the input file
    open(luin,file=input_crd_file,status='old',iostat=ioerr)
  endif
! output file
  call getarg(2,output_crd_file)
  open(luout,file=output_crd_file,status='unknown',iostat=ioerr)
! correlation threshold
  call getarg(3,arg)
  read(arg,*)corrl_limit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the coordinates from the input coordinate file.   ioerr = 0
  n_obs = 0
  do while (ioerr == 0)
    read(luin,*,iostat=ioerr,end=1000)line
    if(ioerr ==0)n_obs = n_obs+1
  enddo

1000 print*,'There are ',n_obs,' entries in the input coordinate file.'
  rewind(luin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate arrays
  allocate(neu(n_obs,3))
  allocate(neu_sig(n_obs,3))
  allocate(midas_vel(n_obs,4))
  allocate(epochs(n_obs))
  allocate(midas_smooth(n_obs,3))
  allocate(corrl(n_obs,3))
  allocate(offsets(n_obs,3))
  midas_vel = 0.0
  midas_smooth = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now read the coordinates. Data needs to be in (unformatted):
!
!  decimal_year   north   east   up   sigma_n   sigma_e   sigma_u
  do iobs = 1,n_obs
    read(luin,*)epochs(iobs),neu(iobs,1:3),neu_sig(iobs,1:3)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now calculate the midas velocity at each epoch, for each coordinate type (NEU)
    do icomp = 1,3
      call midas(n_obs,epochs,neu(:,icomp),midas_vel(:,icomp),midas_vel(:,4))   ! subroutine returns midas velocity calcs plus flag information

      ! calculate the median
      call find_median(n_obs,midas_vel(:,icomp),midas_vel(:,4),median_val(icomp))

      print*,'Midas estimate for component ',icomp,': ',median_val(icomp),midas_vel(median_val(icomp),icomp)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now smooth the MIDAS velocity estimates using a 365 day sliding window
  smooth_window = 1.0d0   ! smooth the data over a 1-year time interval
  print*,"Smoothing midas estimates"
  do icomp = 1,3
    call smooth_data(n_obs,smooth_window,epochs,midas_vel(:,icomp),midas_vel(:,4),midas_smooth(:,icomp))
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, compute the correlation between the smoothed midas time series and  a triangular-shaped
! synthetic time series. Near 100% correlation indicates the presence of an offset.
print*,'compute correlations'
  do icomp = 1,3
    call correlate_for_offsets(n_obs,epochs,midas_smooth(:,icomp),corrl(:,icomp))
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEBUG
!print*,'iobs,orig,smoothed,correlations'
!do iobs = 1,n_obs
!  print*,epochs(iobs),midas_vel(iobs,1),midas_smooth(iobs,1),corrl(iobs,1)," obs and smoothed"
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, identify the large correlation values as being the indicators of offsets. Pass the 
! original coordinate obs in so that the offset can be deduced from the original data.
  do icomp = 1,3
    call identify_offsets(n_obs,icomp,corrl_limit,midas_vel(median_val(icomp),icomp),epochs,corrl(:,icomp),midas_vel(:,icomp) &
                         ,midas_smooth(:,icomp),neu(:,icomp),offsets(:,icomp))
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finally, correct the original time series of coordinates for offsets and output the corrected values
  do iobs = 1,n_obs
    write(luout,100)epochs(iobs),neu(iobs,1)+sum(offsets(1:iobs,1)) &
                    ,neu(iobs,2)+sum(offsets(1:iobs,2)),neu(iobs,3)+sum(offsets(1:iobs,3)),neu_sig(iobs,:),neu(iobs,:) &
                    ,midas_vel(iobs,1:3),midas_smooth(iobs,1:3) &
                    ,corrl(iobs,:)
100 format(f15.5,3f15.3,3f10.3,3f15.3,3f10.3,3f10.3,3f15.7)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end



!************************!!!!!!!!!!!!!!!!!!!!!!**********************************
! Subroutines
!
! midas                 : controls the estimation of midas velocities, one component at a time
! calc_midas_vel        : compute the midas velocities from coordinate pairs separated by 365 days
! find_mediam           : finds the median value of a vector of numbers
! sort_vals             : perform a bubble sort
! smooth_data           : calculate running mean smoothing of a vector of data
! correlate_for_offsets : correlate smoothed midas time series with a triangular-shped function of an offset
! identify_offsets      : use the correlations to identify where the offsets occurred


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine identify_offsets(n_obs,icomp,corrl_limit,median_val,epochs,corrl,midas_vel,midas_smooth,neu,offsets)

! use the computed correlations to identify the epochs at which offsets occurred. These will be where the
! correlations between the obs and a synthetic triangular-shaped signal are really high .... (to be defined!)
!
! P. Tregoning
! 24 June 2020

  implicit none

! passed variables
  integer*4, intent(in)  :: n_obs                ! number of observations
  integer*4, intent(in)  :: icomp                ! coordinate component being assessed
  real*8,    intent(in)  :: corrl_limit          ! threshold over which the correlation must be greater for it to be an offset
  real*8,    intent(in)  :: median_val           ! the midas median value velocity estimates for one component
  real*8,    intent(in)  :: epochs(n_obs)        ! epochs (in decimal years) of observations
  real*8,    intent(in)  :: midas_vel(n_obs)     ! the midas velocities  for one component
  real*8,    intent(in)  :: midas_smooth(n_obs)  ! the smoothed midas velocities  for one component
  real*8,    intent(in)  :: neu(n_obs)           ! original coordinate observations for one component
  real*8,    intent(in)  :: corrl(n_obs)         ! the value of the correlations at each epoch for one component
  real*8,    intent(out) :: offsets(n_obs)       ! the value of the offset at each epoch for one component

! local variables
  integer*4 :: iepoch
  real*8    :: time_offset                         ! time offset between peak correlation and actual epoch of offset in time series
  real*8    :: offset                              ! temporary storage for calculated offset
  integer*4 :: offset_epoch                        ! epoch in coordinate data at which the offset occurred
  integer*4 :: n_corrl                             ! number of epochs where the correlation exceeds the threshold
  integer*4 :: corrl_epochs(1000)                  ! vector to store list of high correlation epochs
  integer*4 :: epoch1,epoch2,epoch3

!! set the value of the time offset between peak correlation and actual offset in time series.
!! value of 2 years (minus 1 day) if using the peak correlation epoch
!  time_offset = 2.d0 - 1.d0/365.d0
!! value of 1 year (minus 1 day) if using the peak midas_smoothed value
  time_offset = 1.d0 - 1.d0/365.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through epochs and identify where correlations indicate offsets
  offsets = 0.d0
  n_corrl = 1
  corrl_epochs(1) = 1
  do iepoch = 1,n_obs
    if(dabs(corrl(iepoch)) > corrl_limit)then
      n_corrl = n_corrl + 1
      corrl_epochs(n_corrl) = iepoch
    endif
  enddo

!DEBUG
!do iepoch = 1,n_corrl
!   print*,iepoch,corrl_epochs(iepoch),epochs(corrl_epochs(iepoch)),corrl(corrl_epochs(iepoch))
!enddo
 
! now, assess the identified epochs and use only the epoch of the peak value if there are neighbouring epochs
  do iepoch = 2,n_corrl-1
    epoch1 = corrl_epochs(iepoch-1)
    epoch2 = corrl_epochs(iepoch)
    epoch3 = corrl_epochs(iepoch+1)
!print*,epoch1,epoch2,epoch3,corrl(epoch1),corrl(epoch2),corrl(epoch3)
    if( (dabs(corrl(epoch3)) > dabs(corrl(epoch2)) .and. epochs(epoch3)-epochs(epoch2) < 0.003d0) &
       .or. (dabs(corrl(epoch1)) > dabs(corrl(epoch2)) .and. epochs(epoch2)-epochs(epoch1) < 0.003d0)) then
    else
      call calc_offset_size(n_obs,epoch2,corrl(epoch2),median_val,epochs,time_offset,midas_vel,midas_smooth,neu,offset &
                           ,offset_epoch)
      if(offset_epoch > 0)then
        offsets(offset_epoch) = offset      
        print*,'Identify_offsets: offset detected on component', icomp,epochs(offset_epoch),corrl(epoch2) &
              ,offsets(offset_epoch),offset_epoch
      endif
    endif

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

  return
  end subroutine identify_offsets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_offset_size(n_obs,iepoch,corrl,median_val,epochs,time_offset,midas_vel,midas_smooth,neu,offset,offset_epoch)

! subroutine to calculate the size of an offset at a particular, known, epoch given the origina
! coordinate information and the midas velocity estimates
!
! P. Tregoning
! 25 June 2020

  implicit none

! passed variables
  integer*4, intent(in)  :: n_obs                ! number of observations
  integer*4, intent(in)  :: iepoch               ! identified epoch of greatest correlation
  real*8,    intent(in)  :: corrl                ! the correlation value identified and passed in
  real*8,    intent(in)  :: median_val           ! midas median value velocity estimate for the component
  real*8,    intent(in)  :: epochs(n_obs)        ! epochs (in decimal years) of observations
  real*8,    intent(in)  :: time_offset          ! time difference between peak of correlation and the epoch of the offset
  real*8,    intent(in)  :: midas_vel(n_obs)     ! the midas velocities  for one component
  real*8,    intent(in)  :: midas_smooth(n_obs)  ! the midas velocities  for one component
  real*8,    intent(in)  :: neu(n_obs)           ! original coordinate observations for one component
  real*8,    intent(out) :: offset               ! the value of the offset at each epoch for one component
  integer*4, intent(out) :: offset_epoch         ! the epoch in the coordinate data at which the offset occurred

  real*8 :: mean
  real*8 :: rms_tmp,rms_best,offset_best
  real*8 :: offset_tmp,offset_tmp1,offset_tmp2
  integer*4 :: ndays,ndays_used, ndays_used1,ndays_used2,iepoch_tmp
  integer*4 :: isearch,dt,offset_epoch_best,iday
  real*8,allocatable :: crd_tmp(:)
  logical :: debug

! variables for determining more accurately the epoch of the offset
  integer*4 :: min_smooth_epoch, max_smooth_epoch,offset_epoch_median
  real*8    :: offset_median    ! difference between the smoothed midas velocity estimate and the median estimate at peak smoothed epoch


debug=.false.

! PT200630: can we identify the peak/trough of the smoothed midas values, somewhere around "iepoch", to identify
!           more accurately a best-guess for the epoch of the offset. It would occur 2 years before the offset., just
!           as the maximum correlation, but it might be more accurate in identifying the epoch
  call offset_epoch_from_smoothed(n_obs,iepoch,epochs,midas_smooth,min_smooth_epoch,max_smooth_epoch)
  if(corrl < 0.d0)then
    offset_epoch = min_smooth_epoch
    offset_epoch_median = min_smooth_epoch
  else
    offset_epoch = max_smooth_epoch
    offset_epoch_median = max_smooth_epoch
  endif
  iepoch_tmp = offset_epoch
  offset_median = midas_smooth(offset_epoch) - median_val

if(debug)print*,'iepoch,offset_epoch',iepoch,offset_epoch,epochs(iepoch),epochs(offset_epoch),corrl,offset_median

! the value of "iepoch" passed in is the epoch of greatest correlation, not the epoch of the offset (which
! occurs 2 years hence). We need to identify the index of the epoch when the offset occurred.
  do while (epochs(offset_epoch) < epochs(iepoch_tmp)+time_offset .and. epochs(offset_epoch) < epochs(n_obs) )
    if(dabs(epochs(offset_epoch+1) - epochs(iepoch_tmp)+time_offset) > 0.001d0) offset_epoch = offset_epoch + 1
  enddo
if(debug)print*,'best guess for actual offset epoch: ',offset_epoch,epochs(offset_epoch),offset_median
    offset_epoch_median = offset_epoch

! the original approach for identifying the offset epoch showed that sometimes it identified one day before/after the
! actual day of the offset. This then yields a "spike" in the corrected time series because the offset correction is not
! applied at the correct epoch.

! To attempt to counteract this, we will compute and apply offsets for "offset_epoch" as determined above, plus for one 
! day before and one day after. The epoch that yields the smallest RMS of the corected coordinate time series will be selected
! as having had the offset applied at the correct epoch. This offset estimate (and epoch) will be returned from this subroutine.

! estimate the size of the required offset by averaging 5 data points before and after, then calculating the difference

  ndays = 10 
  allocate(crd_tmp(2*ndays))
  rms_best = 1.d6
  if(offset_epoch + ndays+1 < n_obs)then  ! there are enough obs available to calculate an offset
    do dt = -3,3
      ! calculate the offset, being the difference in the mean value before/after for "ndays" of data
      ! PT200626: only include midas_vel values for which there was a computed midas velocity!
      offset_tmp1 = 0.d0
      offset_tmp2 = 0.d0
      ndays_used1 = 0
      ndays_used2 = 0
      do iday = 1,ndays    ! we want "ndays" worth of useable data before the chosen epoch 
        if(dabs(midas_vel(offset_epoch+dt-ndays+iday-1)) < 300.d0)then  ! it is a valid midas velocity estimate
          ndays_used1 = ndays_used1 + 1
          offset_tmp1 = offset_tmp1 + midas_vel(offset_epoch+dt-ndays+iday-1)
          crd_tmp(ndays_used1) = midas_vel(offset_epoch+dt-ndays+iday-1)
if(debug)then
  print*,iday,epochs(offset_epoch+dt),'pre-offset use ',epochs(offset_epoch+dt-ndays+iday-1) &
        ,midas_vel(offset_epoch+dt-ndays+iday-1),midas_smooth(offset_epoch+dt-ndays+iday-1)
endif
        endif
      enddo

      do iday = 1,ndays    ! we want "ndays" worth of useable data after and including the chosen epoch 
        if(dabs(midas_vel(offset_epoch+dt+iday-1)) < 300.d0)then  ! it is a valid midas velocity estimate
          ndays_used2 = ndays_used2 + 1
          offset_tmp2 = offset_tmp2 + midas_vel(offset_epoch+dt+iday-1)
          crd_tmp(ndays_used2+ndays_used1) = midas_vel(offset_epoch+dt+iday-1)
if(debug)print*,iday,epochs(offset_epoch+dt),'pre-offset use ',epochs(offset_epoch+dt+iday-1),midas_vel(offset_epoch+dt+iday-1)
        endif
      enddo

      ! the offset is the difference between the two calculations
      offset_tmp = offset_tmp2/dble(ndays_used2) - offset_tmp1/dble(ndays_used1)
if(debug)then
  print*,epochs(offset_epoch+dt),' offset calculation: ',offset_tmp,offset_tmp2/dble(ndays_used2),offset_tmp1/dble(ndays_used1)
endif
      ! add the offset from the offset epoch onwards
      do iday=1,ndays_used2
if(debug)print*,'crd_tmp pre-correction',iday,ndays_used1+iday,crd_tmp(ndays_used1+iday)
        crd_tmp(ndays_used1+iday) = crd_tmp(ndays_used1+iday) - offset_tmp  ! the midas offset is -ive for a +ive crd offset!
      enddo




      ! compute the mean value of this 10-day time series
      ndays_used = ndays_used1+ndays_used2
      mean = sum(crd_tmp(1:ndays_used))/dble(ndays_used)
      ! compute the rms
      rms_tmp = 0.d0
      do iday = 1,ndays_used
if(debug)print*,iday,'offset-corrected midas:',crd_tmp(iday)
        rms_tmp = rms_tmp + (crd_tmp(iday)-mean)**2
      enddo

      rms_tmp = dsqrt( rms_tmp/(dble(ndays_used)) )
if(debug)then
  print*,epochs(offset_epoch),epochs(offset_epoch+dt),dt,offset_tmp,' rms: ',rms_tmp
endif

      ! check whether it is the best of the three tested epochs
      if(rms_tmp < rms_best)then
        rms_best = rms_tmp
        offset_best = offset_tmp
        offset_epoch_best = offset_epoch+dt
      endif

    enddo
    if(debug)print*,'Offset of ',offset_best,' identified as being at epoch ',epochs(offset_epoch_best),offset_epoch_best

  else
    print*,'Offset epoch is too close to the end of the input data time span',offset_epoch,n_obs
    offset = 0.d0
    offset_epoch = 0
    return
  endif

! transfer the best estimate to the returned variable
  offset_epoch = offset_epoch_best
  offset = offset_best

! try using the offset at the peak of the smoothed midas
  offset = -1.d0*offset_median
  offset_epoch = offset_epoch_median

  deallocate(crd_tmp)

  return
  end subroutine calc_offset_size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine offset_epoch_from_smoothed(n_obs,iepoch,epochs,midas_smooth,min_epoch,max_epoch)

! subroutine to refine the best-guess of the epoch of the offset, making use of peak/trough of the smoothed
! midas velocities to identify the epoch exactly two years before the epoch of the offset. This should be more
! accurate than using the correlation max, which can be off by several epochs
!
! P. Tregoning
! 30 June 2020

  implicit none

! passed variables
  integer*4, intent(in)  :: n_obs                ! number of observations
  integer*4, intent(in)  :: iepoch               ! identified epoch of greatest correlation
  real*8,    intent(in)  :: epochs(n_obs)        ! epochs (in decimal years) of observations
  real*8,    intent(in)  :: midas_smooth(n_obs)  ! the smoothed midas velocities  for one component
  integer*4, intent(out) :: min_epoch            ! the minimum epoch in the coordinate data at which the smoothed midas peak/trough occurred
  integer*4, intent(out) :: max_epoch            ! the minimum epoch in the coordinate data at which the smoothed midas peak/trough occurred

  integer*4  :: epoch_span
  integer*4  :: i
  real*8     :: max_val,min_val

  logical :: debug

! we want to search through the smoothed midas velocities, centred on "iepoch" and find the maximum (or minimum) value, which 
! will be the epoch two years before the actual epoch of the offset.

! start 50 epochs before and go 50 epochs after "iepoch" and search for the highest (lowest) smoothed value
  epoch_span = 500
  min_val = 1.d6
  max_val = -1.d6
  do i = 1,epoch_span
    if(iepoch+i > 0 .and. iepoch+i < n_obs) then  
      if(midas_smooth(iepoch+i) > max_val .and. midas_smooth(iepoch+i) < 9998.d0)then
        max_epoch = iepoch+i
        max_val = midas_smooth(iepoch+i)
      endif
      if(midas_smooth(iepoch+i) < min_val .and. midas_smooth(iepoch+i) > -9998.d0)then
        min_epoch = iepoch+i
        min_val = midas_smooth(iepoch+i)
      endif
    endif
  enddo


!  print*,'min midas_smoothed epoch:',epochs(min_epoch),min_epoch,midas_smooth(min_epoch)
!  print*,'max midas_smoothed epoch:',epochs(max_epoch),max_epoch,midas_smooth(max_epoch)


  return
  end subroutine offset_epoch_from_smoothed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine correlate_for_offsets(n_obs,epochs,midas_obs,corrl)

! subroutine to generate a synthetic tirangular-shaped function and correlate it across the
! smoothed midas time series. WIth near 100% corrleation (or anti-correlation) we can identify
! the locations of offsets.
!
! P. Tregoning
! 24 June 2020

  implicit none

! passed variables
  integer*4, intent(in)  :: n_obs             ! number of observations
  real*8,    intent(in)  :: epochs(n_obs)     ! epochs (in decimal years) of observations
  real*8,    intent(in)  :: midas_obs(n_obs)  ! the smoothed midas velocities 
  real*8,    intent(out) :: corrl(n_obs)      ! the value of the correlations at each epoch

! local variables
  integer*4  :: ival1,ival2,i
  real*8     :: final_epoch,start_epoch,end_epoch
  real*8     :: midas_mean,synth_mean
  real*8,allocatable :: epochs2(:)
  real*8,allocatable :: synth_obs(:)
  real*8     :: dt,sum1,sum2,sum3

  allocate(epochs2(n_obs))
  allocate(synth_obs(n_obs))

  corrl = 0.d0
! store the epoch of the last measurement
  final_epoch = epochs(n_obs) -1.d0  ! we don't have midas estimates within 1 year of the end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through the epochs
  do ival1 = 1,n_obs-1
    start_epoch = epochs(ival1)
    end_epoch = start_epoch + 2.d0  ! use data of input time series up to 2 years hence
    ival2 = 1
    do while (epochs(ival1+ival2-1) <= end_epoch .and. epochs(ival1+ival2-1) < final_epoch)
      ! create a synthetic tiangular-shaped time series value for this epoch
      epochs2(ival2) = epochs(ival1+ival2-1)
      dt = epochs2(ival2)-start_epoch
      if( dt < 1.d0)then
         synth_obs(ival2) = dt * 5.d0
      else
         synth_obs(ival2) = 5.d0 - (dt-1.d0)*5.d0
      endif

      ival2 = ival2 + 1

    enddo
    ival2 = ival2-1

    ! calculate the means of the two time series
    midas_mean = sum(midas_obs(ival1:ival1+ival2-1))/(dble(ival2))
    synth_mean = sum(synth_obs(1:ival2))/(dble(ival2))
    sum1 = 0.d0
    sum2 = 0.d0
    sum3 = 0.d0
    do i=1,ival2   ! the number of observations in time series from ival1 to two years hence
      sum1 = sum1 + (midas_obs(ival1+i-1)-midas_mean)*(synth_obs(i)-synth_mean)
      sum2 = sum2 + (midas_obs(ival1+i-1)-midas_mean)**2
      sum3 = sum3 + (synth_obs(i)-synth_mean)**2
    enddo
    corrl(ival1) = sum1/(sqrt(sum2)*sqrt(sum3))
    
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine midas(n_obs,epochs,crds,midas_vel,midas_flag)

! subroutine to control the computation of the midas velocities for an input time series of values. Adapted
! from code in midas_v2.f90
!
! P. Tregoning
! 24 June 2020

  implicit none

  integer*4, intent(in)  :: n_obs              ! number of observations in time series
  real*8,    intent(in)  :: epochs(n_obs)      ! epochs (in decimal years) of observations
  real*8,    intent(in)  :: crds(n_obs)        ! coordinate time series
  
  real*8,    intent(out) :: midas_vel(n_obs)   ! estimated midas velocities 
  real*8,    intent(out) :: midas_flag(n_obs)  ! flags indicating status of midas velocity estimates

! local variables
  integer*4  :: iepoch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iepoch = 1,n_obs

! check that we are not into the last j*365 epochs
    if( epochs(iepoch) > epochs(n_obs) - 1.d0) then
      midas_vel(iepoch) = 9999.
      midas_flag(iepoch) = 9999.
    else
      call calc_midas_vel(iepoch,n_obs,epochs,crds,1,midas_vel(iepoch))
      ! update the flag if there was no midas velocity
      if(midas_vel(iepoch) < -9998.)then
        midas_flag(iepoch) = -9999.
      endif

    endif
  enddo  ! end of epoch loop

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!**********!!!!!!!!!!!***********!!!!!!!!!!!**********!!!!!!!!
  subroutine calc_midas_vel(iepoch,n_obs,epoch,neu,step,midas_vel)

! subroutine to step forward "step" year(s) and calculate the velocity from the difference
! in position. If there is no estimate exactly "step" years ahead then return a velocity of -9999.
!
! P. Tregoning
! 21 December 2015

  implicit none

  integer*4, intent(in) :: iepoch             ! particular epoch of the coordinate vector for which to compute the midas velocity
  integer*4, intent(in) :: step               ! time interval (in decimal years) over which to compute the coordinate difference
  integer*4, intent(in) :: n_obs              ! number of coordinate observations
  real*8    :: epoch(n_obs)                   ! array of epoch information (yr,doy,decimal_year)
  real*8    :: neu(n_obs)                     ! the input coordinate values
  real*8    :: midas_vel

! local variables
  integer*4 :: add_day
  real*8    :: dec_yr1,dt
  logical   :: found

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the decimal year epoch equivalent to "iepoch"
  dec_yr1 = epoch(iepoch)

! read through the values to see whether we can find one that is within one day of "step" years ahead
  found = .false.
  add_day = iepoch
  dt = 0.
  do while (.not. found .and. dt < step + 0.003 )  ! 0.0027 years is one day
    add_day = add_day + 1
    dt = (epoch(add_day) ) - dec_yr1
    if(dt >= step .and. dt < step + 0.003) then   ! we are "step" years ahead of our starting point
      found = .true.
      midas_vel = neu(add_day) - neu(iepoch) / step
    endif
  enddo

! did we find an epoch into the future
  if (.not. found) then
!    print*,'No value found at epoch ',epoch(iepoch+step,3),icomp
    midas_vel = -9999.
  endif

  return
  end
!!!!!!**********!!!!!!!!!!!***********!!!!!!!!!!!**********!!!!!!!!


!!!!!!**********!!!!!!!!!!!***********!!!!!!!!!!!**********!!!!!!!!
  subroutine find_median(nvals,vals,flags,median_val)

! subroutine to find the median of a set of numbers
!
! P. Tregoning
! 21 December 2015

  implicit none

  integer*4 , intent(in)  :: nvals               ! the number of values
  real*8    , intent(in)  :: vals(nvals)         ! the values
  real*8    , intent(in)  :: flags(nvals)        ! flags as to whether no value (-9999.), good (0.0), outside requested range (-888.) or too close to end (9999.)
  integer*4 , intent(out) :: median_val          ! pointer to the median value

! local variables
  real*8,allocatable    :: tmpvals(:),tmpflags(:)
  integer*4,allocatable :: pointers(:)
  integer*4             :: iepoch
  integer*4             :: n_midas(2)            ! number of epochs (1) without a midas velocity and (2) too close to the end

  allocate(tmpvals(nvals))
  allocate(tmpflags(nvals))
  allocate(pointers(nvals))

!  sort the values to arrange in ascending order
  tmpvals = vals
  tmpflags = flags
  do iepoch = 1,nvals
    pointers(iepoch) = iepoch
  enddo
  call sort_vals(nvals,tmpvals(1:nvals),pointers(1:nvals),tmpflags,n_midas)

!do iepoch=1,nvals
!  print*,'iepoch,vals(iepoch),tmpvals(iepoch),pointers(iepoch)',iepoch,vals(iepoch),tmpvals(iepoch),pointers(iepoch)
!enddo
!stop 'stop here'

! now, find the median value. We need to ignore values of -9999. (no midas_vel estimated) and 9999. (last year)
  median_val = pointers(int( (n_midas(1)+n_midas(2))/2.0))
!  print*,'pointer and median val',pointers(int( (n_midas(1)+n_midas(2))/2.0)),vals(median_val)

  deallocate(tmpvals)
  deallocate(tmpflags)
  deallocate(pointers)

  return
  end
!!!!!!**********!!!!!!!!!!!***********!!!!!!!!!!!**********!!!!!!!!



!!!!!!**********!!!!!!!!!!!***********!!!!!!!!!!!**********!!!!!!!!
  subroutine sort_vals(nvals,vals,pointers,flags,n_midas)

! write a bubble sort from scratch because I am in a plane and can't download a pre-written bubble sort ..... :-(
!
! P. Tregoning
! 21 December 2015

  implicit none

  integer*4  :: nvals
  integer*4  :: n_nomidas(2)    ! the number of epochs either a) without a midas velocity or b) within a year of the end
  real*8     :: vals(nvals)
  real*8     :: flags(nvals)    ! good (0), no midas vel (-9999.), outside range (-888.), too close to end (999.)
  real*8     :: tmpval,tmp_flag
  integer*4  :: pointers(nvals),tmp_pointer
  integer*4  :: iepoch,i,j
  integer*4  :: n_midas(2)       ! number of epochs with a midas velocity estimate


! try to shuffle things around ....
  do j=1,nvals
    do i=1,nvals-1
      if(vals(i) > vals(i+1))then
! swap the values
        tmpval    = vals(i+1)
        vals(i+1) = vals(i)
        vals(i)   = tmpval
! swap the pointers
        tmp_pointer  = pointers(i+1)
        pointers(i+1) = pointers(i)
        pointers(i)   = tmp_pointer
! swap the flags
        tmp_flag = flags(i+1)
        flags(i+1) = flags(i)
        flags(i) = tmp_flag
      endif
    enddo
  enddo

! now, strip out the -9999. and 9999. estimates.
  n_midas(1) = -9999
  n_midas(2) = 0

! PT160102: we now have -888. as an indicator of a velocity outside the epoch range in which to compute the midas velocity
  do iepoch = 1,nvals
    if(n_midas(1) == -9999 .and. flags(iepoch) > -888.)then
      n_midas(1) = iepoch
    endif
    if(n_midas(2) == 0 .and. flags(iepoch) > 9998.)then
      n_midas(2) = iepoch - 1
    endif
  enddo

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine smooth_data(n_obs,smooth_window,epochs,obs_orig,obs_flag,obs_smooth)

! subroutine to compute a sliding window smoothing of a time series of observations
!
! P. Tregoning
! 24 June 2020

  implicit none

! passed variables
  integer*4,  intent(in)  :: n_obs             ! number of observations
  real*8,     intent(in)  :: smooth_window     ! time interval over which to smooth the data points
  real*8,     intent(in)  :: epochs(n_obs)     ! epochs of observations (in decimal_years)
  real*8,     intent(in)  :: obs_orig(n_obs)   ! original, unsmoothed observations
  real*8,     intent(in)  :: obs_flag(n_obs)   ! flag of whether midas velocity estimate is trustworthy
  real*8,     intent(out) :: obs_smooth(n_obs) ! smoothed observations

! local variables
  integer*4 :: iobs,itmp,counter
  real*8    :: tmpsum
  real*8    :: min_smooth,max_smooth   ! track the max and min smoothed values
  integer*4 :: imax_smooth,imin_smooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  min_smooth = 1.d6
  max_smooth = -1.d6
  iobs = 1
  obs_smooth = 0.d0
  do while (epochs(iobs) < epochs(n_obs) - smooth_window)
!print*,'smooting observation window starting at ',iobs,epochs(iobs)
    itmp = iobs
    counter = 0
    tmpsum = 0.d0
    do while (epochs(itmp) < epochs(iobs)+smooth_window)
      if(epochs(itmp) < epochs(iobs)+smooth_window )then
        if(obs_orig(itmp) < 9999.d0 .and. obs_orig(itmp)> -9999.d0)then
          counter = counter + 1
          tmpsum = tmpsum + obs_orig(itmp)
        endif

        itmp = itmp + 1
          
      endif
    enddo
    obs_smooth(iobs) = tmpsum/dble(counter)

    ! check if it is a min/max value
    if(obs_smooth(iobs) > max_smooth)then
      max_smooth = obs_smooth(iobs)
      imax_smooth = iobs
    else if(obs_smooth(iobs) < min_smooth)then
      min_smooth = obs_smooth(iobs)
      imin_smooth = iobs
    endif

!print*,'smoothed: ',epochs(iobs),obs_smooth(iobs)
    ! increment to the next observation
    iobs = iobs + 1
  enddo

!  print*,'epochs of max and min smoothed values: ',imax_smooth,imin_smooth,epochs(imax_smooth),epochs(imin_smooth)
  return
  end subroutine smooth_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





