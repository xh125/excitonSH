module result_io
  implicit none

  integer ,parameter :: dp=kind(1.0d0)
  integer, parameter, public      :: maxlen = 120
  real(kind=dp),parameter         :: Au2fs  = 2.418884326505d-2
  integer :: nnode,inode
  integer :: ncore,icore
  !integer :: iaver,naver
  integer, public, save           :: stdout
  character(len=maxlen)           :: stdfile
  !! Unit on which stdout is written
  !! Max column width of input file
  real(kind=dp) :: dt
  integer       :: nstep,istep
  integer       :: iaver,naver
  integer       :: isnap,nsnap
  integer       :: ibasis,nbasis
  integer       :: ifreem,nfreem
  
  integer                         :: ierr
  logical                         :: lexist
  character(len=maxlen)           :: msg
  character(len=maxlen)           :: home_dir
  
  
  real(kind=dp),allocatable :: pes_exciton(:,:,:),inf_elec(:,:,:),inf_hole(:,:,:)
  real(kind=dp),allocatable :: csit_elec(:,:),csit_elec_tmp(:,:)
  real(kind=dp),allocatable :: csit_hole(:,:),csit_hole_tmp(:,:)
  real(kind=dp),allocatable :: wsit_elec(:,:),wsit_elec_tmp(:,:)
  real(kind=dp),allocatable :: wsit_hole(:,:),wsit_hole_tmp(:,:)
  real(kind=dp),allocatable :: psit_elec(:,:),psit_elec_tmp(:,:)
  real(kind=dp),allocatable :: psit_hole(:,:),psit_hole_tmp(:,:)  
  real(kind=dp),allocatable :: xsit(:,:),xsit_tmp(:,:)
  real(kind=dp),allocatable :: ksit(:,:),ksit_tmp(:,:)
  real(kind=dp),allocatable :: ipr_elec(:),ipr_elec_tmp(:)
  real(kind=dp),allocatable :: ipr_hole(:),ipr_hole_tmp(:)
  
  character(len=maxlen) ::  pes_name,csit_name,wsit_name,xsit_name,psit_name,ksit_name
  integer               ::  pes_unit,csit_unit,wsit_unit,xsit_unit,psit_unit,ksit_unit
  character(len=maxlen) ::  inf_name,msd_name,msds_name,ipr_name
  integer               ::  inf_unit,msd_unit,msds_unit,ipr_unit  
  
  namelist / saveinput / nnode,ncore,dt,nstep,nsnap,nbasis,nfreem
  
  contains
  
  !========================================
  subroutine io_error ( error_msg )
  !========================================
  !! Abort the code giving an error message 
  !========================================

    implicit none
    character(len=*), intent(in) :: error_msg

    write(stdout,*)  'Exiting.......' 
    write(stdout, '(1x,a)') trim(error_msg)    
    close(stdout)    
    write(*, '(1x,a)') trim(error_msg)
    write(*,'(A)') "Error: examine the output/error file for details" 
    STOP
         
  end subroutine io_error  
  
  function io_file_unit() !得到一个当前未使用的unit，用于打开文件
  !===========================================                                     
  !! Returns an unused unit number
  !! so we can later open a file on that unit.                                       
  !===========================================
  implicit none

    integer :: io_file_unit,unit_index
    logical :: file_open

    unit_index = 9
    file_open  = .true.
    do while ( file_open )
      unit_index = unit_index + 1
      inquire( unit=unit_index, OPENED = file_open ) !用于检查文件状态，参考P.536
    end do
    
    io_file_unit = unit_index

    return
  
  end function io_file_unit
  
  subroutine open_file(file_name,file_unit)
    implicit none
    
    character(len=*),intent(in) :: file_name
    integer,intent(in)          :: file_unit
    open(unit=file_unit, file=file_name,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem opening "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif
  end subroutine open_file
  
  subroutine close_file(file_name,file_unit)
    implicit none
    
    integer,intent(in)          :: file_unit
    character(len=*),intent(in) :: file_name
    close(file_unit,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem close "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif   
  end subroutine close_file
  
  subroutine saveresult()
  implicit none
    
    inquire(directory = './result',exist=lexist)
    if (.not. lexist) call system('mkdir ./result')
    
    pes_unit = io_file_unit()
    pes_name = './result/pes.out'
    call open_file(pes_name,pes_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(pes_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(pes_exciton(ibasis,isnap,iaver),ibasis=-1,nbasis)
      enddo
    enddo
    call close_file(pes_name,pes_unit)
   
    inf_unit = io_file_unit()
    inf_name = './result/inf_elec.out'
    call open_file(inf_name,inf_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(inf_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(inf_elec(ibasis,isnap,iaver),ibasis=1,3)
      enddo
    enddo
    call close_file(inf_name,inf_unit)
    
    inf_unit = io_file_unit()
    inf_name = './result/inf_hole.out'
    call open_file(inf_name,inf_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(inf_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(inf_hole(ibasis,isnap,iaver),ibasis=1,3)
      enddo
    enddo
    call close_file(inf_name,inf_unit)
    
    csit_unit = io_file_unit()
    csit_name = './result/csit_elec.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(csit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    csit_unit = io_file_unit()
    csit_name = './result/csit_hole.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(csit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    wsit_unit = io_file_unit()
    wsit_name = './result/wsit_elec.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(wsit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)

    wsit_unit = io_file_unit()
    wsit_name = './result/wsit_hole.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(wsit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)
    
    psit_unit = io_file_unit()
    psit_name = './result/psit_elec.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(psit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    psit_unit = io_file_unit()
    psit_name = './result/psit_hole.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(psit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    xsit_unit = io_file_unit()
    xsit_name = './result/xsit.out'
    call open_file(xsit_name,xsit_unit)
    do isnap=1,nsnap
      write(xsit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(xsit(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(xsit_name,xsit_unit)
    
    ksit_unit = io_file_unit()
    ksit_name = './result/ksit.out'
    call open_file(ksit_name,ksit_unit)
    do isnap=1,nsnap
      write(ksit_unit,'(999999e12.5)') dt*nstep*isnap*Au2fs,(ksit(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call  close_file(ksit_name,ksit_unit)
    
    !msd_unit = io_file_unit()
    !msd_name = 'msd.out'
    !call open_file(msd_name,msd_unit)
    !do isnap=1,nsnap
    !  write(msd_unit,'(999999e12.5)') dt*nstep*isnap*au2fs,msd(isnap)
    !enddo
    !call close_file(msd_name,msd_unit)
   
    ipr_unit = io_file_unit()
    ipr_name = './result/ipr_elec.out'
    call open_file(ipr_name,ipr_unit)
    do isnap=1,nsnap
      write(ipr_unit,'(999999e12.5)') dt*nstep*isnap,ipr_elec(isnap)
    enddo
    call close_file(ipr_name,ipr_unit)
    
    ipr_unit = io_file_unit()
    ipr_name = './result/ipr_hole.out'
    call open_file(ipr_name,ipr_unit)
    do isnap=1,nsnap
      write(ipr_unit,'(999999e12.5)') dt*nstep*isnap,ipr_hole(isnap)
    enddo
    call close_file(ipr_name,ipr_unit)
  
  end subroutine saveresult  
  
end module

program main
  use result_io
  implicit none
  
  integer :: itmp,itmp1,itmp2,itmp3
  character(len=maxlen) :: ctmp,ctmp1,ctmp2,ctmp3
  logical:: ltmp
  real(kind=dp) :: rtmp
  integer :: incar_unit
  character(len=maxlen) :: incar_name
  
  stdout = io_file_unit()
  stdfile= "saveresult.out"
  call open_file(stdfile,stdout)
  
  !! initial input parameter
  ncore = 20
  nnode = 40
  dt    = 0.1
  nstep = 10
  nsnap = 2000
  nbasis= 3*3*(3*3*2)
  nfreem= 3*3*2-3
  
  incar_unit = io_file_unit()
  incar_name = "saveincar"
  call open_file(incar_name,incar_unit)
  read(unit=incar_unit,nml=saveinput,iostat=ierr,iomsg=msg)
  if(ierr /= 0) then
    call io_error('Error: Problem reading namelist file SHIN')
    call io_error(msg)
  endif
  
  naver = ncore*nnode
  allocate(pes_exciton(-1:nbasis,1:nsnap,1:naver))
  allocate(inf_elec(1:3,1:nsnap,1:naver),inf_hole(1:3,1:nsnap,1:naver))
  allocate(csit_elec(nbasis,nsnap),csit_elec_tmp(nbasis,nsnap))
  allocate(csit_hole(nbasis,nsnap),csit_hole_tmp(nbasis,nsnap))
  allocate(wsit_elec(nbasis,nsnap),wsit_elec_tmp(nbasis,nsnap))
  allocate(wsit_hole(nbasis,nsnap),wsit_hole_tmp(nbasis,nsnap))  
  allocate(psit_elec(nbasis,nsnap),psit_elec_tmp(nbasis,nsnap))
  allocate(psit_hole(nbasis,nsnap),psit_hole_tmp(nbasis,nsnap))
  allocate(xsit(nfreem,nsnap),ksit(nfreem,nsnap))
  allocate(ipr_elec(nsnap),ipr_elec_tmp(nsnap))
  allocate(ipr_hole(nsnap),ipr_hole_tmp(nsnap))
  
  !!initial arrays
  pes_exciton   = 0.0
  inf_elec      = 0.0
  inf_hole      = 0.0
  csit_elec     = 0.0
  csit_elec_tmp = 0.0
  csit_hole     = 0.0
  csit_hole_tmp = 0.0
  wsit_elec     = 0.0
  wsit_elec_tmp = 0.0
  wsit_hole     = 0.0
  wsit_hole_tmp = 0.0
  psit_elec     = 0.0
  psit_elec_tmp = 0.0
  psit_hole     = 0.0
  psit_hole_tmp = 0.0
  xsit          = 0.0
  ksit          = 0.0
  ipr_elec      = 0.0
  ipr_elec_tmp  = 0.0
  ipr_hole      = 0.0
  ipr_hole_tmp  = 0.0
  
  do inode=1,nnode
    write(ctmp1,*) inode
    do icore=1,ncore
      write(ctmp2,*) icore
      iaver = (inode-1)*ncore+icore
      
      pes_unit=io_file_unit()
      pes_name="node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/pes.out"
      call open_file(pes_name,pes_unit)
      write(ctmp,*) nbasis+3
      ctmp = "("//trim(adjustl(ctmp))//"e12.5"//")"
      !ctmp = "("//trim(adjustl(ctmp))//"(e12.5,1X)"//")"
      do isnap=1,nsnap
        read(pes_unit,ctmp) rtmp,(pes_exciton(ibasis,isnap,iaver),ibasis=-1,nbasis)
        write(ctmp2,*) isnap
        write(*,*)"read isnap="//trim(adjustl(ctmp2))//" is OK!"
      enddo
      call close_file(pes_name,pes_unit)
      
      inf_unit = io_file_unit()
      inf_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/inf_elec.out"
      call open_file(inf_name,inf_unit)
      do isnap=1,nsnap
        read(inf_unit,'(4e12.5)') rtmp,(inf_elec(ibasis,isnap,iaver),ibasis=1,3)
      enddo
      call close_file(inf_name,inf_unit)
      
      inf_unit = io_file_unit()
      inf_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/inf_hole.out"
      call open_file(inf_name,inf_unit)
      do isnap=1,nsnap
        read(inf_unit,'(4e12.5)') rtmp,(inf_hole(ibasis,isnap,iaver),ibasis=1,3)
      enddo
      call close_file(inf_name,inf_unit)      
      
      csit_unit = io_file_unit()
      csit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/csit_elec.out"
      call open_file(csit_name,csit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(csit_unit,*) rtmp,(csit_elec_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(csit_name,csit_unit)
      csit_elec = csit_elec + csit_elec_tmp
      
      csit_unit = io_file_unit()
      csit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/csit_hole.out"
      call open_file(csit_name,csit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(csit_unit,*) rtmp,(csit_hole_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(csit_name,csit_unit)     
      csit_hole = csit_hole + csit_elec_tmp
      
      wsit_unit = io_file_unit()
      wsit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/wsit_elec.out"
      call open_file(wsit_name,wsit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(wsit_unit,*) rtmp,(wsit_elec_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(wsit_name,wsit_unit)
      wsit_elec = wsit_elec + wsit_elec_tmp
     
      wsit_unit = io_file_unit()
      wsit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/wsit_hole.out"
      call open_file(wsit_name,wsit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(wsit_unit,*) rtmp,(wsit_hole_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(csit_name,csit_unit)           
      wsit_hole = wsit_hole + wsit_hole_tmp
      
      psit_unit = io_file_unit()
      psit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/psit_elec.out"
      call open_file(psit_name,psit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(wsit_unit,*) rtmp,(psit_elec_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(psit_name,psit_unit)
      psit_elec = psit_elec + psit_elec_tmp
     
      psit_unit = io_file_unit()
      psit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/psit_hole.out"
      call open_file(psit_name,psit_unit)
      write(ctmp,*) nbasis+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(wsit_unit,*) rtmp,(psit_hole_tmp(ibasis,isnap),ibasis=1,nbasis)
      enddo
      call close_file(csit_name,csit_unit)           
      psit_hole = psit_hole + psit_hole_tmp      
      
      xsit_unit = io_file_unit()
      xsit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/xsit.out"
      call open_file(xsit_name,xsit_unit)
      write(ctmp,*) nfreem+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(xsit_unit,*) rtmp,(xsit_tmp(ifreem,isnap),ifreem=1,nfreem)
      enddo
      call close_file(xsit_name,xsit_unit)
      xsit = xsit + xsit_tmp
      
      ksit_unit = io_file_unit()
      ksit_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/ksit.out"
      call open_file(ksit_name,ksit_unit)
      write(ctmp,*) nfreem+1
      ctmp = "'("//trim(adjustl(ctmp))//"e12.5"//")'"
      do isnap=1,nsnap
        read(ksit_unit,*) rtmp,(ksit(ibasis,isnap),ifreem=1,nfreem)
      enddo
      call close_file(ksit_name,ksit_unit)          
      ksit = ksit + ksit_tmp
      
      ipr_unit = io_file_unit()
      ipr_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/ipr_elec.out"
      call open_file(ipr_name,ipr_unit)
      do isnap=1,nsnap
        read(ipr_unit,'(2e12.5)') rtmp,ipr_elec_tmp(isnap)
      enddo
      call close_file(ipr_name,ipr_unit)
      ipr_elec = ipr_elec + ipr_elec_tmp
      
      ipr_unit = io_file_unit()
      ipr_name = "node"//trim(adjustl(ctmp1))//"/core"//trim(adjustl(ctmp2))//"/result/ipr_hole.out"
      call open_file(ipr_name,ipr_unit)
      do isnap=1,nsnap
        read(ipr_unit,'(2e12.5)') rtmp,ipr_hole_tmp(isnap)
      enddo
      call close_file(ipr_name,ipr_unit)
      ipr_hole = ipr_hole + ipr_elec_tmp      
      
    enddo
  enddo
  
  csit_elec=csit_elec/naver
  csit_hole=csit_hole/naver
  wsit_elec=wsit_elec/naver
  wsit_hole=wsit_hole/naver
  psit_elec=psit_elec/naver
  psit_hole=psit_hole/naver
  xsit     =xsit/naver
  ksit     =ksit/naver
  ipr_elec =ipr_elec/naver
  ipr_hole =ipr_hole/naver  
  
  call saveresult()
  
end program