module set_para
!$$$   module documentation block
!                .      .    .                                       .
! module:  set_para
! prgmmr:  Xu Li          date: 2019-03-13
!
! abstract: This module contains variables to process the SST/Tf analysis
!
! program history log:
! 
! Subroutines Included:
!   sub init_grdmod   - initialize grided related variables to default values
!
! Variable Definitions:
!   def atime         - analysis time: yyyymmddhh
!   def lputsi        - logical, fill land grids or not
!   dsearch           - the search radius in KM 
!   sfcio             - I/O mode (nems or nc)
!   nx                - the x dimension of the rtg-like SST file
!   ny                - the y dimension of the rtg-like SST file
!$$$ end documentation block

  character(10) :: catime
  character(10) :: sfcio
  logical :: lputsi
  integer :: nx,ny
  real    :: dsearch
  
contains

  subroutine init_setup
!
! abstract:  set defaults for observation related variables
!
    implicit none

!   Initialize arrays used in namelist obs_input 
    catime='1981010100'
    sfcio='nemsio'
    lputsi=.false.
    dsearch = 400.0
    nx = 4320
    ny = 2160
  end subroutine init_setup
  
end module set_para


