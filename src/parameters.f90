module parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use kinds,only : dp
  use control
  use constants,only : maxlen
  implicit none
  
  integer,public    :: nband,nmode,ntota,ntotk,ntotq
  !体系能带数，声子谱支数，晶格总数，k点总数，q点总数
  integer,public    :: nbasis,nfreem,ibasis,ifreem
  !基矢总数，总的简正自由度
  !integer,public    :: naver,nsnap,nstep
  integer,public    :: iaver,isnap,istep
  real(kind=dp)     :: E_humo,E_lumo,E_max_V,E_min_C
  integer           :: index_humo,index_lumo,index_max_V,index_min_C
  integer           :: Nocc,Nunocc
  
  real(kind=dp),allocatable ::  ph_Q(:),ph_Q0(:),ph_P(:),ph_P0(:)
  !简正坐标，正则动量
  real(kind=dp),allocatable :: HH_e(:,:),HH_h(:,:)
  real(kind=dp),allocatable :: E_e(:),E_h(:),P_e(:,:),P_h(:,:)
  !电子空穴本征值与本征矢
  real(kind=dp),allocatable :: E0_e(:),E0_h(:),P0_e(:,:),P0_h(:,:)
                                   
  logical :: lelecsh,lholesh,lelecholesh,lexcitonsh
  
  
  character(len=maxlen) ::  ctmp
  integer               ::  itmp
  logical               ::  ltmp
  real(kind=dp)         ::  rtmp  
  !临时数据
  logical :: lallocate,lexist
  integer :: iostate,istate
  character(len=maxlen) :: msg
  character(len=9)  :: cdate,ctime
  real(kind=dp)     ::t0,t1,t2,t3,time0,time1,time2,time3  
  
	namelist / shinput / &
          shtype,dimention,na1,na2,na3,nk1,nk2,nk3,nq1,nq2,nq3,         &
          representation,nshiftstep,dtadq,pj_adjust,                    &
          temp,gamma,dt,nstep,nsnap,naver,                              &  
          initnmstat,Lrandomsita,L_hotphonon,hot_mode,hot_scal,         &
          Num_occupied,initehstat,initek,initeb,inithk,inithb,initeWF,  &
          inithWF,initeEN,inithen,initeES,inithES,                      &
          L_exciton,epsr,                                               &
          MSH,Ldecoherece,Tdecoherence,Lfeedback,                       &
          mkl_threads,Lrunsh

  contains
  
  subroutine treat_parameters()
    use constants
    implicit none	    
    !! Change to Hartree atomic units
    temp          = temp/Au2k
    gamma         = gamma*Au2ps
    kb            = k_B_SI/Au2J*Au2k !in Hartree atomic units kb = 1
    dt            = dt/Au2fs
    initeEN       = initeEN/Au2eV
    inithEN       = inithEN/Au2eV
    Nocc          = na1*na2*na3*Num_occupied
    Nunocc        = nbasis - Nocc
    
    if (shtype == "elec") then
      lelecsh=.TRUE.
      lholesh=.FALSE.
      lelecholesh=.FALSE.
      lexcitonsh = .FALSE.
    elseif(shtype == "hole") then
      lelecsh=.FALSE.
      lholesh=.TRUE.
      lelecholesh=.FALSE.
      lexcitonsh = .FALSE.
    elseif(shtype=="elechole") then
      lelecsh=.TRUE.
      lholesh=.TRUE.
      lelecholesh=.TRUE.
      lexcitonsh = .FALSE.
    elseif(shtype == "exciton") then
      lelecsh=.TRUE.
      lholesh=.TRUE.
      lelecholesh=.FALSE.
      lexcitonsh = .TRUE.
    endif
    
  end subroutine  treat_parameters  
    
  
end module parameters   
  