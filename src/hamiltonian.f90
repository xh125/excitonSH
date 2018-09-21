module hamiltonian
  use kinds,only : dp,dpc
  use parameters
  use constants
  use readwannierfile
  use io,only : io_error,io_file_unit,open_file,close_file
  
  implicit none
  real(kind=dp)                 :: rdotk
  complex(kind=dpc)             :: fac  
  real(kind=dp),allocatable     :: H0(:,:),Hep(:,:,:),HH(:,:),H_tmp(:,:)
  !平衡位置哈密顿量，电声耦合项，总的H和临时H
  real(kind=dp),allocatable     :: E(:),P(:,:),E0(:),P0(:,:)
  !哈密顿量的本征值与本征矢
  
  real(kind=dp),allocatable     :: coulomb(:,:,:,:,:)
  
  real(kind=dp),allocatable     :: Hr_1d(:,:,:),Hr_2d(:,:,:,:),&
                                   Hr_3d(:,:,:,:,:)
  complex(kind=dpc),allocatable :: Hrc_1d(:,:,:),Hrc_2d(:,:,:,:),&
                                   Hrc_3d(:,:,:,:,:)
  real(kind=dp),allocatable     :: Hrep_1d(:,:,:,:),Hrep_2d(:,:,:,:,:),&
                                   Hrep_3d(:,:,:,:,:,:)
  complex(kind=dpc),allocatable :: Hrepc_1d(:,:,:,:),Hrepc_2d(:,:,:,:,:),&
                                   Hrepc_3d(:,:,:,:,:,:)                                   
  real(kind=dp),allocatable     :: Tij_1d_tmp(:,:,:),Tij_2d_tmp(:,:,:,:),&
                                   Tij_3d_tmp(:,:,:,:,:)
  real(kind=dp),allocatable     :: Tij_1d_0(:,:,:),Tij_1d_ep(:,:,:,:),&
                                   Tij_2d_0(:,:,:,:),Tij_2d_ep(:,:,:,:,:),&
                                   Tij_3d_0(:,:,:,:,:),Tij_3d_ep(:,:,:,:,:,:)
  logical,allocatable           :: adj_Tij_1d(:,:,:),adj_Tij_2d(:,:,:,:),&
                                   adj_Tij_3d(:,:,:,:,:)  
  
  
  contains
  
  !================================================!
  !根据dimention参数，将 Tij_3(21)d_tmp数组写入    !
  !Tij_name                                        !
  ! Tij_3(21)d_tmp ->Tij_name                      !
  !================================================!
  subroutine write_Tij(Tij_name,dimention)
    implicit none
    character(len=maxlen),intent(in) :: Tij_name
    integer,intent(in)               :: dimention
    integer :: Tij_unit
    integer :: m,n
    Tij_unit = io_file_unit()
    call open_file(Tij_name,Tij_unit)
    if(dimention==3) then
      write(Tij_unit,"(1X,5A5,A12)") &
            " ir1 "," ir2 "," ir3 ","  m  ","  n  "," Tij "
    elseif(dimention==2) then
      write(Tij_unit,"(1X,4A5,A12)") &
            " ir1 "," ir2 ","  m  ","  n  "," Tij " 
    elseif(dimention==1) then
      write(Tij_unit,"(1X,3A5,A12)") &
            " ir1 ","  m  ","  n  "," Tij "
    endif
    
    do ir3=-1,1
      do ir2=-1,1
        do ir1=-1,1
          do m=1,num_wann
            do n=1,num_wann
              if(dimention==3) then
                write(Tij_unit,"(1X,5I5,F12.6)") &
                    ir1,ir2,ir3,m,n,Tij_3d_tmp(n,m,ir1,ir2,ir3)
              elseif(dimention==2 .and. ir3==0) then
                write(Tij_unit,"(1X,4I5,F12.6)") &
                    ir1,ir2,m,n,Tij_2d_tmp(n,m,ir1,ir2)
              elseif(dimention==1 .and. ir3==0 .and. ir2==0) then
                write(Tij_unit,"(1X,3I5,F12.6)") ir1,m,n,Tij_1d_tmp(n,m,ir1)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
    
    call close_file(Tij_name,Tij_unit)
    
  end subroutine write_Tij
  
  !===================================================!
  ! 读取Tij文件，并将Tij参数写入Tij_3(21)d_tmp数组    !
  ! Tij_name -> Tij_3(21)d_tmp                        !
  !===================================================!
  subroutine read_Tij(Tij_name,dimention)
    implicit none
    character(len=maxlen),intent(in) :: Tij_name
    integer,intent(in) :: dimention
    integer :: Tij_unit
    integer :: m,n
    Tij_unit = io_file_unit()
    call open_file(Tij_name,Tij_unit)
    read(Tij_unit,*) 
    do ir3=-1,1
      do ir2=-1,1
        do ir1=-1,1
          do m=1,num_wann
            do n=1,num_wann
              if(dimention==3) then
                read(Tij_unit,"(1X,T27,F12.6)") Tij_3d_tmp(n,m,ir1,ir2,ir3)
              elseif(dimention==2 .and. ir3==0) then
                read(Tij_unit,"(1X,T22,F12.6)") Tij_2d_tmp(n,m,ir1,ir2)
              elseif(dimention==1 .and. ir3==0 .and. ir2==0) then
                read(Tij_unit,"(1X,T17,F12.6)") Tij_1d_tmp(n,m,ir1)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
    call close_file(Tij_name,Tij_unit)
    
  end subroutine read_Tij  
  
  !=======================================!
  !通过读取wannier90_hr文件得到Ham_r数组 =!
  !Hr_name ->Ham_r                        !
  !=======================================!
  subroutine getHam(Hr_name)
    use readwannierfile
    implicit none
    character(len=maxlen),intent(in) :: Hr_name
    
    call readwannhr(Hr_name)
    
  end subroutine getHam  

  !==============================================================================!
  !根据实空间的Ham_r(num_wann,num_wann,nrpts) 通过傅里叶变换得到倒空间某一K点处的!
  !哈密顿量 Ham_kprm(num_wann,num_wann)                                          !
  !Ham_r() + kpoint(3) -> Ham_kprm(num_wann,num_wann)                            !
  ! ham_kprm^k(i,j) = SUM_R( (e^i2pi*R*k)*ham_r(i,j,R)/ndegen(R) )               !
  !==============================================================================!
  subroutine Ham2Hamkprm(kkpoint)
    use readwannierfile
    implicit none
    real(kind=dp),intent(in) :: kkpoint(3)
    real(kind=dp)            :: SUM2U
    !real(kind=dp) :: rdotk
    !complex(kind=dpc) :: fac
    lallocate = Allocated(Ham_kprm)
    if(.not. lallocate) allocate(Ham_kprm(num_wann,num_wann))
    lallocate = Allocated(U_int)
    if(.not. lallocate) allocate(U_int(num_wann,num_wann))
    lallocate = Allocated(E_k)
    if(.not. lallocate ) allocate(E_k(num_wann))
    do irpt=1,nrpts
      rdotk=twopi*dot_product(kkpoint(:),irvec(:,irpt))
      fac=cmplx(cos(rdotk),sin(rdotk),dp)/real(ndegen(irpt),dp)
      ham_kprm=ham_kprm+fac*ham_r(:,:,irpt)
    end do
      U_int = ham_kprm
    call dia_heH(num_wann,ham_kprm,E_k,U_int)  
    SUM2U = SUM(CONJG(U_int(:,1))*U_int(:,1))
  end subroutine Ham2Hamkprm
  
  !==================================================!
  !将Ham数组根据dimention参数，设置Tij_3(21)d_tmp数组!
  ! Ham_r -> Tij_3(21)d_tmp                          !
  !==================================================!
  subroutine Ham2Tij(dimention)
    implicit none
    integer,intent(in) :: dimention
    integer :: m_wann,n_wann,m,n
    do irpt=1, nrpts       !R
      ir1=irvec(1,irpt)
      ir2=irvec(2,irpt)
      ir3=irvec(3,irpt)      
      do n=1, num_wann   !n
        do m=1, num_wann  !m
          ReH = real(Ham_r(m,n,irpt))
          ImH = Aimag(Ham_r(m,n,irpt))
          if(dimention==3) then
            if( abs(ir1)<=1 .and. abs(ir2)<=1 .and. abs(ir3)<=1 ) then
            Tij_3d_tmp(m,n,ir1,ir2,ir3) = ReH
            endif
            
          elseif(dimention==2) then
            if( abs(ir1)<=1 .and. abs(ir2)<=1 .and. ir3==0) then
              Tij_2d_tmp(m,n,ir1,ir2) = ReH
            endif
          
          elseif(dimention==1) then
            if( abs(ir1)<=1 .and. ir2==0 .and. ir3==0) then
              Tij_1d_tmp(m,n,ir1) = ReH
              !<m0|H|nR>
              !if(abs(ReH) >= pj_adjust) then
                !adj_Tij_1d(m,n,ir1) = .TRUE.
              !endif
            endif
          endif
        enddo
      enddo
    enddo                    
    
  end subroutine Ham2Tij 
    
  
  !========================================!
  !通过Ham_r数组，设置Hr_3(21)d数组        !
  !Ham_r ->  Hr_3(21)d                     !
  !========================================!
  subroutine Ham2Hr(dimention)
    implicit none
    integer,intent(in) :: dimention
    logical :: lexist
    integer :: n_wann,m_wann,n,m
    if(dimention == 3) then
      Hr_3d = 0.0d0
      Hrc_3d=cmplx_0
    elseif(dimention==2) then
      Hr_2d =0.0d0
      Hrc_2d=cmplx_0
    elseif(dimention==1) then
      Hr_1d = 0.0d0
      Hrc_1d=cmplx_0
    endif
    
    !read <m0|H|nR>
    !in wannier90_hr.dat R,m,n,HR,HI
    outer:do irpt=1, nrpts       !R
      ir1=irvec(1,irpt)
      ir2=irvec(2,irpt)
      ir3=irvec(3,irpt)      
      do n=1, num_wann   !n
        do m=1, num_wann  !m  
          if(dimention==3) then
            if(ir1>na1 .or. ir1<-na1 .or. ir2>na2 .or. ir2<-na2) cycle outer
            if(ir3>na3 .or. ir3<-na3) cycle outer
            Hr_3d(m,n,ir1,ir2,ir3) = real(Ham_r(m,n,irpt))
            Hrc_3d(m,n,ir1,ir2,ir3)= Ham_r(m,n,irpt)
            
          elseif(dimention==2) then
            if(ir1>na1 .or. ir1<-na1 .or. ir2>na2 .or. ir2<-na2) cycle outer
            Hr_2d(m,n,ir1,ir2) = real(Ham_r(m,n,irpt))
            Hrc_2d(m,n,ir1,ir2)= Ham_r(m,n,irpt)
            
          elseif(dimention==1) then
            if(ir1>na1 .or. ir1<-na1) cycle outer         
            Hr_1d(m,n,ir1) = real(Ham_r(m,n,irpt))
            Hrc_1d(m,n,ir1)= Ham_r(m,n,irpt)           
          endif
          
        enddo
      enddo
    enddo outer     
  end subroutine  
  
  !=========================================!
  !通过Hr_3(21)d数组,以及na1,na2,na3,设置   !
  !hamiltonian H_tmp 用于得到H0，Hep，HH    !
  ! Hr_3(21)d -> H_tmp(H0)
  !=========================================!
  subroutine Hr2Htmp(dimention)
    !! read Hr(Hr_1d,Hr_2d,Hr_3d) to H_tmp(H0 unit)
    use readwannierfile
    implicit none
    integer,intent(in)               :: dimention
    integer :: ia1_l,ia1_r,ia2_l,ia2_r,ia3_l,ia3_r
    integer :: n1,n2,n3,im,jm
      
      if(dimention == 3) then
        H_tmp    = 0.0
        do ia3_l=0,na3-1
          do ia3_r=0,na3-1
            n3=ia3_r-ia3_l        
            do ia2_l=0,na2-1
              do ia2_r=0,na2-1
                n2=ia2_r-ia2_l
                do ia1_l=0,na1-1
                  do ia1_r=0,na1-1
                    n1=ia1_r-ia1_l
              
                    im=(ia3_l*na2*na1+ia2_l*na1+ia1_l)*num_wann
                    jm=(ia3_r*na2*na1+ia2_r*na1+ia1_r)*num_wann
                    H_tmp(im+1:im+num_wann,jm+1:jm+num_wann)=Hr_3d(:,:,n1,n2,n3)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        
        
      elseif(dimention == 2 ) then       
        H_tmp    = 0.0
        !<m0|H|nR>
        !<m0|
        do ia2_l=0,na2-1
          do ia2_r=0,na2-1
            n2=ia2_r-ia2_l
            do ia1_l=0,na1-1
              do ia1_r=0,na1-1
                n1=ia1_r-ia1_l
              
                im=(ia2_l*na1+ia1_l)*num_wann
                jm=(ia2_r*na1+ia1_r)*num_wann
                H_tmp(im+1:im+num_wann,jm+1:jm+num_wann)=Hr_2d(:,:,n1,n2)
              enddo
            enddo
          enddo
        enddo
        
      elseif(dimention == 1) then
        H_tmp    = 0.0
        do ia1_l=0,na1-1
          do ia1_r=0,na1-1
            n1=ia1_r-ia1_l          
            im=(ia1_l)*num_wann
            jm=(ia1_r)*num_wann
            H_tmp(im+1:im+num_wann,jm+1:jm+num_wann)=Hr_1d(:,:,n1)
          enddo
        enddo   
      endif
        
  end subroutine Hr2Htmp
  
  
  !======================================!
  !用于得到平衡位置处的H0(nbasis,nbasis) !
  !in Hartrre unit                       !
  !======================================!
  subroutine getH0(nbasis,represtation,dimention)
    use readwannierfile
    implicit none
    integer,intent(out)              :: nbasis
    character(len=maxlen),intent(in) :: represtation
    integer,intent(in)               :: dimention
    character(len=maxlen)            :: H0parafile,wanncentfile
    integer :: ia1_l,ia1_r,ia2_l,ia2_r,ia3_l,ia3_r
    integer :: n1,n2,n3,im,jm
    if( trim(adjustl(represtation))== "wfstat" ) then
      !H0parafile -> H0
      H0parafile = './wannier/wannier90_hr.dat'
      wanncentfile = './wannier/wannier90_centres.xyz'
      inquire(directory = './Hr_wannier',exist=lexist)
      if(.not. lexist) call system('mkdir ./Hr_wannier')
      inquire(file="./Hr_wannier/wannier90_hr.dat",exist=lexist)
      if(.not. lexist) then
        call system('cp ./wannier/wannier90_hr.dat ./Hr_wannier/wannier90_hr.dat')
      endif
      H0parafile = "./Hr_wannier/wannier90_hr.dat"
      !H0parafile -> num_wann
      call set_num_wann(H0parafile,num_wann)
      !num_wann -> Hr Hrc
      if(dimention==3) then
        allocate(Hr_3d(num_wann,num_wann,-na1:na1,-na2:na2,-na3:na3))
        allocate(Hrc_3d(num_wann,num_wann,-na1:na1,-na2:na2,-na3:na3))    
        Hr_3d   = 0.0
        Hrc_3d  = cmplx_0        
      elseif(dimention==2) then
        allocate(Hr_2d(num_wann,num_wann,-na1:na1,-na2:na2))
        allocate(Hrc_2d(num_wann,num_wann,-na1:na1,-na2:na2))        
        Hr_2d  = 0.0
        Hrc_2d = cmplx_0          
      elseif(dimention==1) then
        allocate(Hr_1d(num_wann,num_wann,-na1:na1))    
        allocate(Hrc_1d(num_wann,num_wann,-na1:na1))
        Hr_1d  = 0.0 
        Hrc_1d = cmplx_0   
      endif                
      !wanncentfile -> Rwann
      call readWFcentre(wanncentfile)
      !读取wannier函数的中心
      nband  = num_wann
      nbasis = na1*na2*na3*num_wann
      allocate(H0(nbasis,nbasis),H_tmp(nbasis,nbasis),HH(nbasis,nbasis))
      allocate(E(nbasis),P(nbasis,nbasis))
      !H0parafile -> Ham_r
      call getHam(H0parafile)     ! get Ham_r
      allocate(Ham_r_0(num_wann,num_wann,nrpts))
      Ham_r_0=Ham_r
      !Ham_r ->  Hr_3(21)d 
      call Ham2Hr(dimention)      ! Ham_r -> Hr
      !通过读取wannier90_hr.dat 文件得到Ham_r(num_wann,num_wann,nrpts)
      !并根据dimention 和Ham_r 设置Hr和Hrc
      !Hr_3(21)d ->H_tmp
      call Hr2Htmp(dimention)     ! Hr -> H_tmp -> H0
      !根据Hr或者Hrc写出体系hamitonian人临时文件H_tmp
      H0 = H_tmp
      H0 = H0/AU2EV
      
    elseif( trim(adjustl(represtation))== "adiabatic" ) then
    
    elseif( trim(adjustl(represtation))== "blochstat" ) then
    
    elseif( trim(adjustl(represtation))== "atomob" ) then
    
    elseif( trim(adjustl(represtation))== "molerob" ) then
    
    else
    
    endif
  
  end subroutine getH0
  
  !=======================================!
  !用于得到Hep(nbasis,nbasis,nfreem)      !
  !in Hartrre unit                        !
  !=======================================!
  subroutine getHep(represtation,dimention,nfreem)
    use datafitting
    use readwannierfile
    implicit none
    character(len=maxlen),intent(in) :: represtation
    integer,intent(in)               :: dimention
    integer,intent(in)               :: nfreem
    integer                          :: ifreem
    character(len=maxlen)            :: ctmpfreem,ctmpldQ,Hr_name
    character(len=maxlen)            :: Hep_name,Hep_fitname
    integer                          :: Hep_unit
    integer:: istat
    integer:: ndQ,idQ
    complex(kind=dpc),allocatable    :: Ham_r_phQ(:,:,:,:)
    real(kind=dp),allocatable        :: Ham_r_slope(:,:,:),Ham_r_yint(:,:,:),&
                                        Ham_r_r(:,:,:)
    character(len=maxlen) :: err_msg
    integer :: ia1_l,ia1_r,ia2_l,ia2_r,ia3_l,ia3_r
    integer :: n1,n2,n3,im,jm
    real(kind=dp) :: slope,y_int,r
    real(kind=dp),allocatable :: ldQ(:),HdQ(:)
    
    allocate(Hep(nbasis,nbasis,nfreem),stat=istat,errmsg=err_msg)
    if(istat/=0) call io_error(err_msg)
    Hep = 0.0d0
    if( trim(adjustl(represtation))== "wfstat" ) then
      do ifreem=1,nfreem
        write(ctmpfreem,*) ifreem
        Hep_name = "./Hr_wannier/wannier90_hr_"//trim(adjustl(ctmpfreem))//".dat"
        inquire(file=Hep_name,exist=lexist)
        if(lexist) then
          call getHam(Hep_name)
          call Ham2Hr(dimention)
          call Hr2Htmp(dimention)
          Hep(:,:,ifreem) = H_tmp
        else
          !建议修改为通过不同的dQ的值，得到的Ham_r进行线性拟合得到的电声耦合常数。
          !进行线性拟合
          ndQ=2*nshiftstep+1
          lallocate = allocated(Ham_r_phQ)
          if(.not. lallocate) allocate(Ham_r_phQ(num_wann,num_wann,nrpts,ndQ))
          lallocate = allocated(Ham_r_slope)
          if(.not. lallocate) allocate(Ham_r_slope(num_wann,num_wann,nrpts))
          lallocate = allocated(Ham_r_yint)
          if(.not. lallocate) allocate(Ham_r_yint(num_wann,num_wann,nrpts))  
          lallocate = allocated(Ham_r_r)
          if(.not. lallocate) allocate(Ham_r_r(num_wann,num_wann,nrpts))          
          lallocate = allocated(ldQ)
          if(.not. lallocate) allocate(ldQ(ndQ))
          lallocate = allocated(HdQ)
          if(.not. lallocate) allocate(HdQ(ndQ))
          Ham_r_phQ = cmplx_0
          do idQ=1,ndQ
            ldQ(idQ)=(idQ-nshiftstep-1)*dtadQ
            write(ctmpldQ,"(F8.4)") ldQ(idQ)
            Hr_name= "./nomashift/noma_"//trim(adjustl(ctmpfreem))//&
                  "/shift_"//trim(adjustl(ctmpldQ))//"/wannier90_hr.dat"
            call getHam(Hr_name)
            Ham_r_phQ(:,:,:,idQ) = Ham_r
          enddo
          
          do irpt=1,nrpts
            do n2=1,num_wann
              do n1=1,num_wann
                HdQ = real(Ham_r_phQ(n1,n2,irpt,:))
                call linefitting(ndQ,ldQ,HdQ,slope,y_int,r)
                Ham_r_slope(n1,n2,irpt) = slope
                Ham_r_yint(n1,n2,irpt)  = y_int
                Ham_r_r(n1,n2,irpt)     = r
                Ham_r(n1,n2,irpt)       = CMPLX(slope)
              enddo
            enddo
          enddo
          
          
          call writewannhr(Hep_name,num_wann,nrpts,ndegen,irvec,Ham_r)
          call Ham2Hr(dimention)
          call Hr2Htmp(dimention)
          Hep(:,:,ifreem) = H_tmp
          
          Hep_fitname = "./Hr_wannier/wannier90_hr_linefit"//trim(adjustl(ctmpfreem))//".dat"   
          call writeHepfit(Hep_fitname,num_wann,nrpts,ndegen,irvec,Ham_r_slope,Ham_r_yint,Ham_r_r)
        endif
      enddo
      
      Hep = Hep/AU2EV*Au2ang*dsqrt(au2amu)
      lallocate = allocated(Ham_r_phQ)
      if(lallocate) deallocate(Ham_r_phQ)
      lallocate = allocated(Ham_r_slope)
      if(lallocate) deallocate(Ham_r_slope)
      lallocate = allocated(Ham_r_yint)
      if(lallocate) deallocate(Ham_r_yint)  
      lallocate = allocated(Ham_r_r)
      if(lallocate) deallocate(Ham_r_r)          
      lallocate = allocated(ldQ)
      if(lallocate) deallocate(ldQ)
      lallocate = allocated(HdQ)
      if(lallocate) deallocate(HdQ)
      
    elseif( trim(adjustl(represtation))== "adiabatic" ) then
    
    elseif( trim(adjustl(represtation))== "blochstat" ) then
    
    elseif( trim(adjustl(represtation))== "atomob" ) then
    
    elseif( trim(adjustl(represtation))== "molerob" ) then
    
    else
    
    endif

  end subroutine getHep
  
  !========================================!
  !用于得到不同site处的电子与空穴之间的库仑!
  !相互作用                                !
  !========================================!
  subroutine get_Coulomb(nbasis,epsr)
    use readposcar
    implicit none
    integer,intent(in)::nbasis
    real(kind=dp),intent(in) :: epsr
    !integer :: ia1_l,ia1_r,ia2_l,ia2_r,ia3_l,ia3_r
    integer :: i_WF,j_WF
    integer :: n1,n2,n3!,im,jm,imwf,jmwf
    real(kind=dp) :: R1(3),R2(3),R3(3)
    real(kind=dp) :: a(3,3)
    real(kind=dp) :: Relec_hole(3),Reh
    a = real_lattice
    lallocate = allocated(coulomb)
    if(.not. lallocate) allocate(coulomb(num_wann,num_wann,1-na1:na1-1,1-na2:na2-1,1-na3:na3-1))
    coulomb = 0.0d0
    do n3=1-na3,na3-1
      R3(:) = n3*real_lattice(3,:)
      do n2=1-na2,na2-1
        R2(:) = n2*real_lattice(2,:)
        do n1=1-na1,na1-1
          R1(:) = n1*real_lattice(1,:)             
          do i_WF=1,num_wann
            do j_WF=1,num_wann
              Relec_hole(:) =R3(:)+R2(:)+R1(:)+ Rwann(:,i_WF)-Rwann(:,j_WF)
              Reh = sqrt(Sum(Relec_hole**2))
              if(Reh < a_lattice ) then
                coulomb(j_WF,i_WF,n1,n2,n3) = 1.0/epsr*a_lattice
              else
                coulomb(j_WF,i_WF,n1,n2,n3) = 1.0/epsr*Reh
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  
  end subroutine get_Coulomb
  
  !===========================================!
  ! HH = H0+ SUM(phQ(ifreem)*Hep(:,:,ifreem))=!
  !===========================================!
  subroutine set_HH(nfreem,nbasis,phQ,H0,Hep,HH)
    implicit none
    integer,intent(in)       :: nfreem,nbasis
    real(kind=dp),intent(in) :: phQ(nfreem),H0(nbasis,nbasis),Hep(nbasis,nbasis,nfreem)
    real(kind=dp),intent(out):: HH(nbasis,nbasis)
    integer :: ifreem
    HH = 0.0d0
    HH = H0
    do ifreem=1,nfreem
      HH = HH + phQ(ifreem)*Hep(:,:,ifreem)
    enddo
    
  end subroutine set_HH
   
  subroutine set_H_with_Coulomb(nbasis,nn_elec,nn_hole,HH_e,HH_h)
    implicit none
    integer,intent(in)          :: nbasis
    real(kind=dp),intent(in)    :: nn_elec(nbasis),nn_hole(nbasis)
    real(kind=dp),intent(inout) :: HH_e(nbasis,nbasis),HH_h(nbasis,nbasis)
    integer :: ia1_l,ia1_r,ia2_l,ia2_r,ia3_l,ia3_r
    integer :: ibasis,jbasis
    integer :: n1,n2,n3
    integer :: i_WF,j_WF
    !lallocate = allocated(coulomb)
    !if(.not. lallocate) allocate(coulomb(num_wann,num_wann,na1,na2,na3))    

    do ia3_l=0,na3-1
      do ia3_r=0,na3-1
        n3=ia3_r-ia3_l 
        do ia2_l=0,na2-1
          do ia2_r=0,na2-1
            n2=ia2_r-ia2_l
            do ia1_l=0,na1-1
              do ia1_r=0,na1-1
                n1=ia1_r-ia1_l
                ibasis = (ia3_l*na2*na1+ia2_l*na1+ia1_l)*num_wann
                jbasis = (ia3_r*na2*na1+ia2_r*na1+ia1_r)*num_wann
                do i_WF=1,num_wann
                  do j_WF=1,num_wann   
                    HH_e(ibasis+i_WF,ibasis+i_WF) = HH_e(ibasis+i_WF,ibasis+i_WF)-&
                    0.5*(nn_hole(jbasis+j_WF)-nn_elec(jbasis+j_WF))*coulomb(j_WF,i_WF,n1,n2,n3)
                    HH_h(ibasis+i_WF,ibasis+i_WF) = HH_h(ibasis+i_WF,ibasis+i_WF)+&
                    0.5*(nn_hole(jbasis+j_WF)-nn_elec(jbasis+j_WF))*coulomb(j_WF,i_WF,n1,n2,n3)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    
  end subroutine set_H_with_Coulomb
  
  !==================================================!
  !得到Tij_3d_0和adj_Tij_3d                          !
  !==================================================!
  subroutine get_Tij_0(dimention)
    use io
    implicit none
    integer,intent(in)  :: dimention
    character(len=maxlen) :: Hr_name
    character(len=maxlen) :: Tij_name
    if(dimention==3) then
      allocate(Tij_3d_tmp(num_wann,num_wann,-1:1,-1:1,-1:1))      
      allocate(Tij_3d_0(num_wann,num_wann,-1:1,-1:1,-1:1))
      allocate(adj_Tij_3d(num_wann,num_wann,-1:1,-1:1,-1:1))  
      Tij_3d_0  = 0.0d0
      adj_Tij_3d     = .False.              
    elseif(dimention==2) then
      allocate(Tij_2d_tmp(num_wann,num_wann,-1:1,-1:1))    
      allocate(Tij_2d_0(num_wann,num_wann,-1:1,-1:1))
      allocate(adj_Tij_2d(num_wann,num_wann,-1:1,-1:1))
      Tij_2d_0  = 0.0d0
      adj_Tij_2d     = .False.
    elseif(dimention==1) then
      allocate(Tij_1d_tmp(num_wann,num_wann,-1:1))    
      allocate(Tij_1d_0(num_wann,num_wann,-1:1))
      allocate(adj_Tij_1d(num_wann,num_wann,-1:1))
      Tij_1d_0  = 0.0d0
      adj_Tij_1d     = .False. 
    endif    
    
    Tij_name = "./Tij_parameter/Tij_0"
    inquire(directory = './Tij_parameter',exist=lexist)
    if(.not. lexist) call system('mkdir ./Tij_parameter')    
    !inquire(file="./Tij_parameter/Tij_0",exist=lexist)
    inquire(file=Tij_name,exist=lexist)
    
    if (lexist) then
      !Tij_name = "./Tij_parameter/Tij_0"
      call read_Tij(Tij_name,dimention)
    else
      !!读取wannier文件，并将Tij参数写入Tij_3(21)d_tmp数组
      Hr_name = "./wannier/wannier90_hr.dat"
      call readwannhr(Hr_name)
      !根据dimention参数，由Ham_r设置Tij_3(21)d_tmp数组
      call Ham2Tij(dimention)
      call write_Tij(Tij_name,dimention)
    endif

    if(dimention==3) then
      Tij_3d_0 = Tij_3d_tmp/AU2EV
      adj_Tij_3d = (Tij_3d_0 > pj_adjust/AU2EV)      
    elseif(dimention==2) then
      Tij_2d_0 = Tij_2d_tmp/AU2EV
      
      adj_Tij_2d = (Tij_2d_0 > pj_adjust/AU2EV)
    elseif(dimention ==1 ) then
      Tij_1d_0 = Tij_1d_tmp/AU2EV
      adj_Tij_1d = (Tij_1d_0 > pj_adjust/AU2EV)
    endif    
      
  end subroutine get_Tij_0
  
  !==================================================!
  !得到 Tij_3(21)d_ep                                !
  !==================================================!
  subroutine get_Tij_dQ(nfreem,dimention)
    use io
    implicit none
    integer,intent(in)  :: dimention
    integer,intent(in)  :: nfreem
    integer             :: ifreem
    integer             :: ndQ,idQ
    real(kind=dp)       :: ldQ
    character(len=maxlen) :: Hr_name,ctmpfreem
    character(len=maxlen) :: Tij_name,ctmpldQ
    if(dimention==3) then  
      allocate(Tij_3d_ep(num_wann,num_wann,-1:1,-1:1,-1:1,1:nfreem)) 
      Tij_3d_ep  = 0.0d0            
    elseif(dimention==2) then
      allocate(Tij_2d_ep(num_wann,num_wann,-1:1,-1:1,1:nfreem))
      Tij_2d_ep  = 0.0d0
    elseif(dimention==1) then
      allocate(Tij_1d_ep(num_wann,num_wann,-1:1,1:nfreem))    
      Tij_1d_ep  = 0.0d0
    endif    
    
    inquire(directory = './Tij_parameter',exist=lexist)
    if(.not. lexist) call system('mkdir ./Tij_parameter')
    ndQ   = 2*nshiftstep+1
    allocate(Ham_r_dQ(num_wann,num_wann,nrpts,ndQ))
    Ham_r_dQ = cmplx_0
    do ifreem=1,nfreem
      write(ctmpfreem,*) ifreem
      Tij_name= "./Tij_parameter/Tij_"//trim(adjustl(ctmpfreem))
      inquire(file=trim(adjustl(Tij_name)),exist=lexist)
      if(lexist) then
        call read_Tij(Tij_name,dimention)
      else
        do idQ=1,ndQ
          ldQ=(idQ-nshiftstep-1)*dtadQ
          write(ctmpldQ,"(F8.4)") ldQ
        
          Hr_name= "./nomashift/noma_"//trim(adjustl(ctmpfreem))//&
                  "/shift_"//trim(adjustl(ctmpldQ))//"/wannier90_hr.dat"
          Ham_r = cmplx_0
          call getHam(Hr_name)
          Ham_r_dQ(:,:,:,idQ) = Ham_r
        enddo
        Ham_r = Ham_r_dQ(:,:,:,ndQ) - Ham_r_dQ(:,:,:,1)
        Ham_r = Ham_r/(2.0*real(nshiftstep)*dtadQ)
        call Ham2Tij(dimention)
        call write_Tij(Tij_name,dimention)
      endif
      
      if(dimention==3) then
        Tij_3d_ep(:,:,:,:,:,ifreem) = Tij_3d_tmp/AU2EV*Au2ang*dsqrt(au2amu)
      elseif(dimention==2) then
        Tij_2d_ep(:,:,:,:,ifreem) = Tij_2d_tmp/AU2EV*Au2ang*dsqrt(au2amu)
      elseif(dimention ==1 ) then
        Tij_1d_ep(:,:,:,ifreem) = Tij_1d_tmp/AU2EV*Au2ang*dsqrt(au2amu)
      endif 
      
    end do      
    
    
  
  end subroutine
  
  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!
  subroutine dia_syH(nbasis,HH,EE,PP)
    use f95_precision
    use lapack95
    implicit none
    !!use LAPACK with Fortran f95 interface
    integer ,intent(in)       :: nbasis
    real(kind=dp),intent(in)  :: hh(nbasis,nbasis)
    real(kind=dp),intent(out) :: ee(nbasis),pp(nbasis,nbasis)
    !pp(:,ibasis) 本征态
  
    pp = hh
    !!     USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    call syev(pp,ee,'V','U') 
    !!On exit, hh array is overwritten
  end subroutine dia_syH  

  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!
  subroutine dia_heH(nbasis,HH,EE,PP)
    use f95_precision
    use lapack95
    implicit none
    !!use LAPACK with Fortran f95 interface
    integer ,intent(in)       :: nbasis
    complex(kind=dpc),intent(in)   ::  hh(nbasis,nbasis)
    complex(kind=dpc),intent(out)  ::  pp(nbasis,nbasis)
    real(kind=dp),intent(out) :: ee(nbasis)
    !pp(:,ibasis) 本征态
  
    pp=hh
    !!     USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    call heev(pp,ee,'V','U')    !P1143 MKL
    !!On exit, hh array is overwritten
  end subroutine dia_heH   
  

  
  !=========================================================!
  != convert wavefunction from diabatic to adiabatic basis =!
  !=========================================================!
  subroutine convert_diabatic_adiabatic(nnbasis,pp,cc,ww)
  
  implicit none
    integer,intent(in)            :: nnbasis
    integer                       :: ibasis,jbasis
    real(kind=dp),intent(in)      :: pp(1:nnbasis,1:nnbasis)
    complex(kind=dpc),intent(in)  :: cc(1:nnbasis)
    complex(kind=dpc),intent(out) :: ww(1:nnbasis)
    real(kind=dp)                 :: sumw2,sumC2
    sumC2 = SUM(CONJG(cc)*cc)
    ww=cmplx_0
    do ibasis=1,nnbasis
      do jbasis=1,nnbasis
        ww(ibasis)=ww(ibasis)+pp(jbasis,ibasis)*cc(jbasis)
      enddo
    enddo
    sumw2 = SUM(CONJG(ww)*ww)
    ww = ww /sqrt(sumw2)
    
  endsubroutine convert_diabatic_adiabatic
  
  
end module hamiltonian