module surfacehopping
  use kinds,only : dp,dpc
  use constants
  use parameters
  use hamiltonian
  implicit none
  
  !integer :: iaver
  integer :: isurface,isurface_e,isurface_h
  integer :: isurface_a,isurface_b,isurface_j,isurface_k
  real(kind=dp),allocatable     :: S_ai(:),S_bi(:)
  !real(kind=dp),allocatable     :: E_e(:),E_h(:)
  real(kind=dp),allocatable     :: dE_dQ_e(:,:),dE_dQ_h(:,:)
  !电子和空穴的本征能量对简正坐标的导数
  complex(kind=dpc),allocatable :: C(:),W(:),W0(:)
  complex(kind=dpc),allocatable :: C_e(:),C_h(:),W_e(:),W_h(:),W0_e(:),W0_h(:)
  !电子和空穴的态在绝热和非绝热基矢下的展开系数
  real(kind=dp)    ,allocatable :: n_e(:),n_h(:)
  !电子和空穴在透热表象下的密度分布
  real(kind=dp)    ,allocatable :: d_e(:,:,:),d_h(:,:,:)
  real(kind=dp)    ,allocatable :: d0_e(:,:,:),d0_h(:,:,:)
  real(kind=dp)    ,allocatable :: g_e(:),g1_e(:),g_h(:),g1_h(:)
  !电子和空穴的非绝热耦合项
  
  !store information
  real(kind=dp),allocatable :: pes_elec(:,:,:),inf_elec(:,:,:),csit_elec(:,:),&
                                      wsit_elec(:,:),psit_elec(:,:),ipr_elec(:)                         
  real(kind=dp),allocatable :: pes_hole(:,:,:),inf_hole(:,:,:),csit_hole(:,:),&
                                      wsit_hole(:,:),psit_hole(:,:),ipr_hole(:)
  real(kind=dp),allocatable :: pes_exciton(:,:,:)  
  real(kind=dp),allocatable :: xsit(:,:),ksit(:,:)
  contains
  
  
  subroutine allocatesh(nfreem,nbasis)
    implicit none
    integer,intent(in) :: nfreem,nbasis
    lallocate = allocated(E_e)
    if(.not. lallocate) allocate(E_e(nbasis))
    lallocate = allocated(E0_e)
    if(.not. lallocate) allocate(E0_e(nbasis))
    lallocate = allocated(E_h)
    if(.not. lallocate) allocate(E_h(nbasis))
    lallocate = allocated(E0_h)
    if(.not. lallocate) allocate(E0_h(nbasis))
    lallocate = allocated(P0_e)
    if(.not. lallocate) allocate(P0_e(nbasis,nbasis))
    lallocate = allocated(P_e)
    if(.not. lallocate) allocate(P_e(nbasis,nbasis))    
    lallocate = allocated(P0_h)
    if(.not. lallocate) allocate(P0_h(nbasis,nbasis))    
    lallocate = allocated(P_h)
    if(.not. lallocate) allocate(P_h(nbasis,nbasis))      
    lallocate = allocated(dE_dQ_e)
    if(.not. lallocate) allocate(dE_dQ_e(nbasis,nfreem))
    lallocate = allocated(dE_dQ_h)
    if(.not. lallocate) allocate(dE_dQ_h(nbasis,nfreem))       
    lallocate = allocated(C)
    if(.not. lallocate) allocate(C(nbasis))
    lallocate = allocated(C_e)
    if(.not. lallocate) allocate(C_e(nbasis))
    lallocate = allocated(C_h)
    if(.not. lallocate) allocate(C_h(nbasis))    
    lallocate = allocated(W)
    if(.not. lallocate) allocate(W(nbasis)) 
    lallocate = allocated(W0)
    if(.not. lallocate) allocate(W0(nbasis)) 
    lallocate = allocated(W)
    if(.not. lallocate) allocate(W_e(nbasis)) 
    lallocate = allocated(W0_e)
    if(.not. lallocate) allocate(W0_e(nbasis))   
    lallocate = allocated(W_h)
    if(.not. lallocate) allocate(W_h(nbasis)) 
    lallocate = allocated(W0_h)
    if(.not. lallocate) allocate(W0_h(nbasis))
     lallocate = allocated(n_e)
    if(.not. lallocate) allocate(n_e(nbasis)) 
    lallocate = allocated(n_h)
    if(.not. lallocate) allocate(n_h(nbasis))
    lallocate = allocated(d_e)
    if(.not. lallocate) allocate(d_e(nbasis,nbasis,nfreem)) 
    lallocate = allocated(d_h)
    if(.not. lallocate) allocate(d_h(nbasis,nbasis,nfreem))
    lallocate = allocated(d0_e)
    if(.not. lallocate) allocate(d0_e(nbasis,nbasis,nfreem)) 
    lallocate = allocated(d0_h)
    if(.not. lallocate) allocate(d0_h(nbasis,nbasis,nfreem))  
    lallocate = allocated(g_e)
    if(.not. lallocate) allocate(g_e(nbasis))
    lallocate = allocated(g1_e)
    if(.not. lallocate) allocate(g1_e(nbasis)) 
    lallocate = allocated(g_h)
    if(.not. lallocate) allocate(g_h(nbasis))  
    lallocate = allocated(g1_h)
    if(.not. lallocate) allocate(g1_h(nbasis)) 
    lallocate = allocated(pes_elec)
    if(.not. lallocate) allocate(pes_elec(0:nbasis,1:nsnap,1:naver))
    lallocate = allocated(inf_elec)
    if(.not. lallocate) allocate(inf_elec(3,nsnap,naver))
    lallocate = allocated(csit_elec)
    if(.not. lallocate) allocate(csit_elec(nbasis,nsnap))
    lallocate = allocated(wsit_elec)
    if(.not. lallocate) allocate(wsit_elec(nbasis,nsnap))
    lallocate = allocated(psit_elec)
    if(.not. lallocate ) allocate(psit_elec(nbasis,nsnap))
    lallocate = allocated(pes_hole)
    if(.not. lallocate) allocate(pes_hole(0:nbasis,1:nsnap,1:naver))
    lallocate = allocated(inf_hole)
    if(.not. lallocate) allocate(inf_hole(3,nsnap,naver))
    lallocate = allocated(csit_hole)
    if(.not. lallocate) allocate(csit_hole(nbasis,nsnap))
    lallocate = allocated(wsit_hole)
    if(.not. lallocate) allocate(wsit_hole(nbasis,nsnap))
    lallocate = allocated(psit_hole)
    if(.not. lallocate ) allocate(psit_hole(nbasis,nsnap))
    lallocate = allocated(xsit)
    if(.not. lallocate) allocate(xsit(nfreem,nsnap))
    lallocate = allocated(ksit)
    if(.not. lallocate) allocate(ksit(nfreem,nsnap))
  end subroutine allocatesh
  !===================================!
  != calculate nonadiabatic coupling =!
  !===================================!
  != ref: notebook page 630          =!
  !===================================!

  subroutine calculate_nonadiabatic_coupling(nnbasis,nnfreem,ee,pp,Hep,dd,dE_dQ)
    implicit none
    integer,intent(in) :: nnfreem,nnbasis
    real(kind=dp),intent(in) :: ee(nnbasis),pp(nnbasis,nnbasis),&
                                Hep(nnbasis,nnbasis,nnfreem)
    real(kind=dp),intent(out):: dd(nnbasis,nnbasis,nnfreem),dE_dQ(nnbasis,nnfreem)
    integer:: ifreem
    integer:: ibasis,jbasis,ik1site,ik2site
    dd=0.0d0
    !dij_qv
    do ibasis=1,nnbasis
      do jbasis=1,ibasis
        if(abs(ee(ibasis)-ee(jbasis))<=1.0/au2ev) then
          do ifreem=1,nnfreem
            do ik1site = 1,nnbasis
              do ik2site = 1,ik1site
                dd(ibasis,jbasis,ifreem) = dd(ibasis,jbasis,ifreem)+pp(ik1site,ibasis)*pp(ik2site,jbasis)*Hep(ik1site,ik2site,ifreem)
              enddo
            enddo
          enddo
          if(jbasis /= ibasis) then
            dd(ibasis,jbasis,:)=dd(ibasis,jbasis,:)/(ee(jbasis)-ee(ibasis))
            dd(jbasis,ibasis,:)= - dd(ibasis,jbasis,:)
            !dij_dQq
          endif
        endif
      enddo
      dE_dQ(ibasis,:) = dd(ibasis,ibasis,:)
    enddo
    
  end subroutine calculate_nonadiabatic_coupling  
  
  !=================================!
  != calculate hopping probability =!
  !=================================!
  != ref: notebook page 631        =!
  !=================================!
  subroutine calculate_hopping_probability(nbasis,nfreem,ww,vv,dd,tt,gg,gg1,isurface)
    implicit none
    integer,intent(in)           :: nbasis,nfreem
    complex(kind=dpc),intent(in) :: ww(nbasis)
    real(kind=dp)    ,intent(in) :: vv(nfreem),dd(nbasis,nbasis,nfreem)
    real(kind=dp)    ,intent(out):: gg(nbasis),gg1(nbasis)
    real(kind=dp)    ,intent(in) :: tt
    integer          ,intent(in) :: isurface
    integer       :: ibasis,jbasis,ifreem
    real(kind=dp) :: sumvd

    gg=0.0d0
    gg1=0.0d0
    !求isurface-> ibasis 之间的跃迁几率
    do ibasis=1,nbasis
      if(ibasis /= isurface) then
        sumvd=0.0d0
        do ifreem=1,nfreem
          sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
        enddo
        !if(llhole) sumvd = -sumvd
        
        gg(ibasis)=2.0d0*tt*real(conjg(ww(isurface))*ww(ibasis))*&
                    sumvd/real(conjg(ww(isurface))*ww(isurface))
        gg1(ibasis)=gg(ibasis)
        if(gg(ibasis) < 0.0d0) gg(ibasis)=0.0d0
        !gg = max{ 0 , gg1}
      endif
    enddo
    
  endsubroutine calculate_hopping_probability  
  
  subroutine correct_hopping_probability(MMSH,nbasis,gg,gg1,E0,ww0,ww,pp0,pp,isurface)
    implicit none
    character(len=maxlen),intent(in) :: MMSH
    integer,intent(in)    :: nbasis,isurface
    real(kind=dp),intent(inout) :: gg(nbasis)
    real(kind=dp),intent(in)    :: gg1(nbasis)
    real(kind=dp),intent(in)    :: E0(nbasis)     
    complex(kind=dpc),intent(in):: ww0(nbasis),ww(nbasis)
    real(kind=dp),intent(in) :: pp0(nbasis,nbasis),pp(nbasis,nbasis)
    real(kind=dp) :: sumg0,sumg1
    real(kind=dp) :: min_dE
    integer       :: ibasis
    integer       :: itrival
    integer       :: max_Sai(1)
    real(kind=dp) :: S_aa
    sumg0=(abs(w0(isurface))**2-abs(w(isurface))**2)/abs(w0(isurface))**2      
    if(trim(adjustl(MMSH)) == 'SC-FSSH') then
      
      if(isurface == 1) then
        itrival = isurface+1
        min_dE  = e0(isurface+1)-e0(isurface)
      elseif(isurface == nbasis) then
        itrival = isurface-1
        min_dE  = e0(isurface)-e0(isurface-1)
      elseif((e0(isurface+1)-e0(isurface)) < (e0(isurface)-e0(isurface-1))) then
        itrival = isurface+1
        min_dE  =(e0(isurface+1)-e0(isurface))
      else
        itrival = isurface-1
        min_dE  = e0(isurface)-e0(isurface-1)
      endif              
      
      GG(itrival) = sumg0 - (SUM(GG1)-GG1(itrival))
      if(GG(itrival) < 0.0d0) GG(itrival) = 0.0d0
      if(SUM(GG) > 1.0d0) GG=GG/SUM(GG)
      
    elseif(trim(adjustl(MMSH)) == "CC-FSSH") then
      lallocate = allocated(S_ai)
      if(.not. lallocate) allocate(S_ai(nbasis))
      S_ai = cmplx_0
      isurface_a = isurface
      S_aa = SUM(pp0(:,isurface_a)*pp(:,isurface_a))
      !S_aa = SUM(CONJG(pp0(:,isurface_a))*pp(:,isurface_a))
      if(S_aa**2 >= 0.5) then
        isurface_j = isurface_a
      else
        do ibasis = 1,nbasis
          S_ai(ibasis) = SUM(pp0(:,isurface_a)*pp(:,ibasis))
          !S_ai(ibasis) = SUM(CONJG(pp0(:,isurface_a))*pp(:,ibasis))
        enddo
        S_ai = S_ai**2
        !S_ai = CONJG(S_ai)*S_ai
        max_Sai =  MAXLOC(S_ai)
        isurface_j = max_Sai(1)
      endif
      if(isurface_j == isurface) then
        gg = gg1
      else
        GG(itrival) = sumg0 - (SUM(GG1)-GG1(itrival))
        if(GG(itrival) < 0.0d0) GG(itrival) = 0.0d0
        if(SUM(GG) > 1.0d0) GG=GG/SUM(GG)
      endif
    
    endif    
     
  end subroutine correct_hopping_probability
  
  subroutine nonadiabatic_transition(MMSH,nbasis,nfreem,ee,pp0,pp,dd,isurface,gg,ww,vv,cc)
    use randoms
    implicit none
    character(len=maxlen),intent(in) :: MMSH
    integer,intent(in) :: nbasis,nfreem
    real(kind=dp),intent(in) :: ee(nbasis)
    real(kind=dp),intent(in) :: pp0(nbasis,nbasis),pp(nbasis,nbasis),dd(nbasis,nbasis,nfreem)
    integer ,intent(inout)   :: isurface
    real(kind=dp),intent(in) :: gg(nbasis)
    complex(kind=dpc),intent(in) :: ww(nbasis),cc(nbasis)
    real(kind=dp),intent(inout) :: vv(nfreem)
    integer:: ibasis,jbasis,ifreem
    integer :: max_Sbi(1)
    real(kind=dp)::sumvd,sumdd,sumgg,flagr,flagd

    call more_random()
    call random_number(flagr)
    sumgg=0.0d0
    do ibasis=1,nbasis
      if(ibasis /= isurface) then
        sumgg=sumgg+gg(ibasis)
        if(flagr < sumgg) then
          isurface_b = ibasis
          sumvd=0.0d0
          sumdd=0.0d0
          do ifreem=1,nfreem
            !sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
            !sumdd=sumdd+dd(isurface,ibasis,ifreem)**2
            sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
            sumdd=sumdd+dd(isurface,ibasis,ifreem)**2
          enddo
          
          !flagd=1.0d0+2.0d0*(ee(isurface)-ee(ibasis))*sumdd/mass/sumvd**2
          flagd=1.0d0+2.0d0*(ee(isurface)-ee(ibasis))*sumdd/sumvd**2
          
          if(trim(adjustl(MMSH))  == "SC-FSSH") then
            if(flagd >= 0.0d0) then
              flagd=(sumvd/sumdd)*(-1.0d0+dsqrt(flagd))
              do ifreem=1,nfreem
                vv(ifreem)=vv(ifreem)+flagd*dd(isurface,ibasis,ifreem)
              enddo
              isurface=ibasis
            endif
          elseif(trim(adjustl(MMSH)) == "CC-FSSH") then
            do jbasis = 1,nbasis
              S_bi(jbasis) = SUM(pp0(:,isurface_b)*pp(:,jbasis))
            enddo
            S_bi = S_bi**2
            max_Sbi = MAXLOC(S_bi)
            isurface_k = max_Sbi(1)
            if(isurface_j == isurface_a) then
              !for type 1 and 3
              if(flagd >= 0.0d0) then
                flagd=(sumvd/sumdd)*(-1.0d0+dsqrt(flagd))
                do ifreem=1,nfreem
                  vv(ifreem)=vv(ifreem)+flagd*dd(isurface,isurface_b,ifreem)
                enddo
                isurface=isurface_k
              endif
            elseif( isurface_b == isurface_k) then
              !for type 2
              if(isurface_b == isurface_j) then
                !type 2 and b = j
                isurface = isurface_j
              else
                !type 2 and b /= j
                if(flagd >= 0.0d0) then
                  flagd=(sumvd/sumdd)*(-1.0d0+dsqrt(flagd))
                  do ifreem=1,nfreem
                    vv(ifreem)=vv(ifreem)+flagd*dd(isurface,isurface_b,ifreem)
                  enddo
                  isurface=isurface_k
                else 
                  isurface = isurface_j
                endif
              endif
            else
              !fot type 4
              if(isurface_b == isurface_j) then
                !type 4 and b = j
                isurface = isurface_j
              else
                !type 4 and b /= j
                if(flagd >= 0.0d0) then
                  flagd=(sumvd/sumdd)*(-1.0d0+dsqrt(flagd))
                  do ifreem=1,nfreem
                    vv(ifreem)=vv(ifreem)+flagd*dd(isurface,isurface_b,ifreem)
                  enddo
                  isurface=isurface_k
                else 
                  isurface = isurface_j
                endif
              endif
              
            endif
                          
          endif

          exit
        endif
      endif
    enddo
    
  endsubroutine nonadiabatic_transition  
  
  
end module surfacehopping