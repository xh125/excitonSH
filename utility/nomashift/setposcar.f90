!#########################setposcar ##############################
!####version 0.01                                                #
!####Last update time 2018.7.23                                  #
!####Write by xiehua                                             #
!####Email : xh125@mail.ustc.edu.cn                              #
!####Function: This program is to shift the POSCAR atoms position#
!####          as (noma mode Vecter /atom_mass**0.5)             #
!####          and then ,calculate the wannier function and TB   #
!####          parameter as a functon of noma shift inorder to   #
!####          get the electro-phonon couping comparaters.       #
!####cankao: J. Chem. Theory Comput. 2018, 14, 3752−3762         #
!#################################################################
module poscarparameter
  implicit none
  
  integer,parameter     ::  dp=kind(1.0D0)
  integer,parameter     ::  maxlen=80
  integer               ::  poscar_unit
  character(len=maxlen) ::  poscar_name
  character(len=maxlen) ::  poscarcomment          !POSCAR  line=1
  real(kind=dp)         ::  scalel                 !POSCAR line=2
  real(kind=dp)         ::  real_lattice(3,3)      !POSCAR line=3,4,5
                                                   !格矢，a11,a12,a13...
  integer               ::  num_atoms,num_species  !POSCAR line=6,7
                                                   !原子总数，元素总数
  character(len=2),allocatable::  atoms_symbol(:)  !line=6,元素符号
  integer,allocatable         ::  atoms_species_num(:) !line=7,各元素原子个数
  character(len=maxlen)       :: Charac            !line=8 原子位置坐标类型
  real(kind=dp),allocatable   ::  amass(:)         !各元素的原子质量
  
  character(len=maxlen) ::  ctmp     !临时字符串
  integer :: ierr,itmp               !临时整数
  integer :: nmode,nstep,i,j,k,m     !ndQ
  
  integer :: iatomstype,iatom,imode,istep
  real(kind=dp)  :: ldQ
  real(kind=dp)  :: dQ
  character(len=maxlen) ::  ctmpstep,ctmpldQ,ctmpmode
  character(len=maxlen) ::  comment,vecterfile
  logical               ::  lexist,lbsub
  
  real(kind=dp),allocatable::positiondQ0(:,:)
  real(kind=dp),allocatable::positiondQ(:,:)
  real(kind=dp),allocatable::phonovecter(:,:,:)
  
  namelist /input/ dQ,nstep,lbsub,amass

  contains
  
  function io_file_unit() !得到一个当前未使用的unit，用于打开文件
    !===========================================
    !                                          
    !! Returns an unused unit number
    !! so we can later open a file on that unit.
    !                                           
    !===========================================
    implicit none

    integer :: io_file_unit,unit_index
    logical :: file_open

    unit_index = 9
    file_open  = .true.
    do while ( file_open )
      unit_index = unit_index + 1
      inquire( unit_index, OPENED = file_open ) !用于检查文件状态，参考P.536
    end do

    io_file_unit = unit_index

    return
  
  end function io_file_unit
  
end module

program mkPOSCAR
  use poscarparameter
  implicit none
  call readPOSCAR()
  call getUserInp()
  call readVecter()
  
  inquire(directory = './nomashift',exist=lexist)
  if(lexist) call system('rm -rf ./nomashift ')
  call system('wait')
  call system('mkdir ./nomashift')
  do imode=1,nmode
    ! mkdir noma_imode
    write(ctmpmode,*) imode
    ctmp = "mkdir ./nomashift/noma_"// trim(adjustl(ctmpmode))
    call system(ctmp)
    write(*,*) ctmp
    do istep=1,2*nstep+1
      ldQ=(istep-1)*dQ-nstep*dQ
      write(ctmpldQ,"(F8.4)") ldQ
      ctmp = "mkdir ./nomashift/noma_"// trim(adjustl(ctmpmode))//&
                    "/shift_"//trim(adjustl(ctmpldQ))
      call system(ctmp)
      ctmp = "cp ./wannier/POTCAR ./nomashift/noma_"// trim(adjustl(ctmpmode))//&
                 "/shift_"//trim(adjustl(ctmpldQ))
      call system(ctmp)
      ctmp = "cp ./wannier/KPOINTS ./nomashift/noma_"// trim(adjustl(ctmpmode))//&
                 "/shift_"//trim(adjustl(ctmpldQ))
      call system(ctmp)      
      ctmp = "cp ./wannierinput/* ./nomashift/noma_"// trim(adjustl(ctmpmode))//&
                 "/shift_"//trim(adjustl(ctmpldQ))
      call system(ctmp)
      ctmp = 'sed  -i "8s/wannier/'//trim(adjustl(ctmpmode))//'_'//trim(adjustl(ctmpldQ))//'/g" '// &
      './nomashift/noma_'// trim(adjustl(ctmpmode))//'/shift_'//trim(adjustl(ctmpldQ))//'/vasp.bsub'
      write(*,*) ctmp
      call system(ctmp)
      call Write_nomal_shift_POSCAR(imode,istep)      
      ctmp = "cd ./nomashift/noma_"// trim(adjustl(ctmpmode))//&
            "/shift_"//trim(adjustl(ctmpldQ))//";bsub < vasp.bsub;cd -"
      if(lbsub) then
        if(istep /= nstep+1) then
          call system(ctmp)
        else
          ctmp = "cp ./wannier/wannier90_hr.dat "//"./nomashift/noma_"// &
          trim(adjustl(ctmpmode))//"/shift_"//trim(adjustl(ctmpldQ))
          call system(ctmp)
        endif
      else
        write(*,*) "NOT calculate nomashift wannier !"
      endif
      
    enddo
  enddo
  
  call freememmery()
end program
  
  subroutine readPOSCAR()
    use poscarparameter
    implicit none         
    
    poscar_name  = './phonon-gamma/POSCAR'
    poscar_unit  = io_file_unit()
    call open_file(poscar_name,poscar_unit)
    read(poscar_unit,*) poscarcomment         !line 1
    read(poscar_unit,"(F16.8)") scalel        !line 2
    !line 3-5
    read(poscar_unit,'(3F20.12)') ((real_lattice(I,J),J=1,3),I=1,3)
    real_lattice = real_lattice*scalel
    
    !line 6-7
    num_atoms    = 0
    num_species  = 1
    ierr = 0
    read(poscar_unit,*)
    do while(ierr <= 0)
      allocate(atoms_species_num(num_species))
      read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
      num_species = num_species + 1
      deallocate (atoms_species_num)
      backspace(poscar_unit)
    end do
    num_species = num_species - 2
    
    allocate(atoms_symbol(num_species))
    allocate(atoms_species_num(num_species))
    backspace(poscar_unit)
    backspace(poscar_unit)
    !line 6-7
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_symbol(i),i=1,num_species)
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
    num_atoms = sum(atoms_species_num)
    nmode = 3 * (num_atoms-1)
    allocate(amass(num_species))
    allocate(positiondQ0(3,num_atoms))
    allocate(positiondQ(3,num_atoms))
    !line 8
    read(poscar_unit,*) Charac
    !line 9-~
    do iatom=1,num_atoms
      read(poscar_unit,*) positiondQ0(:,iatom)
    enddo
    call close_file(poscar_name,poscar_unit)
  
  end subroutine readPOSCAR
  
  subroutine Write_nomal_shift_POSCAR(iimode,iistep)
    use poscarparameter
    implicit none
    integer,intent(in)::iimode,iistep
    poscar_name = "./nomashift/noma_"// trim(adjustl(ctmpmode))//&
                  "/shift_"//trim(adjustl(ctmpldQ))//"/POSCAR"
    poscar_unit = io_file_unit()     
    call open_file(poscar_name,poscar_unit)
    
    !line 1
    write(poscar_unit,'(A30,1X,F8.4)') trim(adjustl(poscarcomment)),ldQ
    !line 2
    write(poscar_unit,"(F8.4)") scalel
    !line 3~5
    write(poscar_unit,'(3F20.12)') ((real_lattice(I,J),J=1,3),I=1,3)
    !do i=1,3
      !write(15,"(3F20.12)") (real_lattice(i,j),j=1,3)
    !end do
    !line 6
    write(ctmp,'(I5)') num_species
    ctmp = "("//trim(adjustl(ctmp))//"A5 )"
    write(poscar_unit,ctmp) (atoms_symbol(i),i=1,num_species)
    !line 7
    write(ctmp,'(I5)') num_species
    ctmp = "("//trim(adjustl(ctmp))//"I5)"
    write(poscar_unit,ctmp) (atoms_species_num(i),i=1,num_species)
    !line 8
    if(Charac(1:1)=='C' .or. Charac(1:1)=='c') then
      write(poscar_unit,*) trim(adjustl(Charac))
      positiondQ(:,:)=positiondQ0(:,:)+ldQ*phonovecter(:,:,iimode)
      do iatom=1,num_atoms
        write(poscar_unit,"(3F16.9)") (positiondQ(j,iatom),j=1,3)
      enddo
      call close_file(poscar_name,poscar_unit)
    else 
    write(poscar_unit,*) "Direct"
      call close_file(poscar_name,poscar_unit)
      write(*,*) "POSCAR error!! need Caracter POSCAR"
      stop
    endif       
  
  end subroutine Write_nomal_shift_POSCAR
  
  subroutine getUserInp()
    use poscarparameter
    implicit none
    integer               ::inp_unit
    character(len=maxlen) ::inp_name
    ! set default values for namelist parameters
    nstep = 5
    dQ    = 0.0020
    lbsub = .False.
    amass(:) = 1.0
    !
    inp_unit  = io_file_unit()
    inp_name  = 'shiftinp'
    !call open_file(inp_name,inp_unit)
    open(unit=inp_unit,file=inp_name,status='unknown',action='read',iostat=ierr)
    if(ierr/=0) then
      write(*,*) "I/O error with input file : 'shiftinp'"
      stop
    endif
    read(inp_unit,nml=input)
    call close_file(inp_name,inp_unit)      
    allocate(phonovecter(3,num_atoms,nmode))
    
  end subroutine getUserInp
  
  subroutine readVecter()
    use poscarparameter
    implicit none
    character(len=maxlen) :: wvecter_name
    integer               :: wvecter_unit
    inquire(file='./phonon-gamma/wvecter.txt',exist=lexist)
    if(.NOT. lexist) then
      itmp = (num_atoms+3)*num_atoms*3+50
      write(ctmp,*) itmp
      ctmp ='grep -A'//trim(adjustl(ctmp))//&
            ' "dynamical matrix" ./phonon-gamma/OUTCAR > ./phonon-gamma/wvecter.txt'
      write(*,*) ctmp
      call system(ctmp)
      !call system('cd ./phono-gamma')
      !call system('grep -A200 "dynamical matrix" ./phonon-gamma/OUTCAR>./phono-gamma/wvecter.txt')
      !call system('wait')
      !call system('cd -')
    endif
    
    wvecter_unit = io_file_unit()
    wvecter_name = "./phonon-gamma/wvecter.txt"
    call open_file(wvecter_name,wvecter_unit)
    !read the first 3 lines.
    read(wvecter_unit,"(1X,//,A50)") comment

    do imode=1,nmode
      read(wvecter_unit,"(1X,//,A50)") comment
      !do iatom=1,num_atoms
        !read(wvecter_unit,*) positiondQ0(:,iatom),phonovecter(:,iatom,imode)
      !end do
      iatom = 0 
      do i=1,num_species
        do j=1,atoms_species_num(i)
          iatom = iatom + 1
          read(wvecter_unit,*) positiondQ0(:,iatom),phonovecter(:,iatom,imode)
          phonovecter(:,iatom,imode) = phonovecter(:,iatom,imode)/dsqrt(amass(i))
        enddo
      enddo
    end do
        
    call close_file(wvecter_name,wvecter_unit)
  
  end subroutine readVecter      
    
  subroutine freememmery()
    use poscarparameter
    implicit none
    
    deallocate(phonovecter)
    deallocate(atoms_symbol,atoms_species_num)
    deallocate(positiondQ0)
    deallocate(positiondQ)
    
  end subroutine freememmery    
  
  subroutine open_file(file_name,file_unit)
    implicit none
    integer,parameter::maxlen=80
    character(len=*),intent(in) :: file_name
    integer,intent(in)          :: file_unit
    integer::ierr
    character(len=maxlen)::msg
    
    open(unit=file_unit, file=file_name,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem opening "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif
  end subroutine open_file
  
  subroutine close_file(file_name,file_unit)    
    implicit none
    integer,parameter::maxlen=80
    integer::ierr
    character(len=maxlen)::msg
    integer,intent(in)          :: file_unit
    character(len=*),intent(in) :: file_name
    close(file_unit,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem close "'//trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif
    
  end subroutine close_file
  
  subroutine io_error ( error_msg )
  !========================================
  !! Abort the code giving an error message 
  !========================================

    implicit none
    character(len=*), intent(in) :: error_msg

    !write(stdout,*)  'Exiting.......' 
    !write(stdout, '(1x,a)') trim(error_msg)    
    !close(stdout)    
    write(*, '(1x,a)') trim(error_msg)
    write(*,'(A)') "Error: examine the output/error file for details" 
    
    STOP
         
  end subroutine io_error

  
  
