module readinput
  use kinds,only : dp
  use constants,only : maxlen
  use parameters
  implicit none
    integer               :: in_unit,tot_num_lines,ierr,loop,in1,in2
    integer               :: num_lines,line_counter
    character(len=maxlen),allocatable:: in_data(:) 
    character(len=maxlen) :: dummy
    integer               :: ipos
    character, parameter  :: TABCHAR = char(9) !char(9)为制表符TAB
  
  contains  
  
  subroutine get_inputfile()
    implicit none
    call treat_inputfile()
    call read_namelist()
  end subroutine

  !=======================================!
  subroutine treat_inputfile()  
  !用于将输入文件中的每一条写入字符串文件，并且改为小写，去除注释          
  !=======================================!
  !! Load the shin file into a character  
  !! array in_file, ignoring comments and  
  !! blank lines and converting everything 
  !! to lowercase characters               
  !=======================================!

    use io,        only : io_file_unit,io_error
    use utility,   only : utility_lowercase
    implicit none
    
    !ia = ichar('a')
    !iz = ichar('z')
    
    in_unit=io_file_unit( )
    open (unit=in_unit, file='SHIN',form='formatted',status='old',iostat=ierr)
    if(ierr /= 0) then
      call io_error('Error: Problem opening input file SHIN')
    endif
    
    num_lines=0;tot_num_lines=0;ierr=0
    do while( ierr == 0 )
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy   !参考p.177 P.540
      if(ierr > 0 ) then
        call io_error('Error: Problem reading input file SHIN')
        call io_error(msg)
      elseif(ierr == 0 )then
        
        ! convert all tabulation characters to spaces
        ipos = index(dummy,TABCHAR) !查询字符串在字符串中出现的位置,并将制表符改为空格
        do while (ipos /= 0)
          dummy(ipos:ipos) = ' '
          ipos = index(dummy,TABCHAR)
        end do
        ! 
        dummy=adjustl(dummy)
        
        tot_num_lines=tot_num_lines+1
        if( dummy(1:1)/='!'  .and. dummy(1:1)/='#' ) then
          if(len_trim(adjustl(dummy)) > 0 ) num_lines=num_lines+1
        endif
      endif
    end do
    !得到SHIN文件中总的行数tot_num_lines以及非注释和空行 num_lines

    rewind(in_unit)

    allocate(in_data(num_lines),stat=ierr)  !字符串数组，内部文件 line=449
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

    line_counter=0
    do loop=1,tot_num_lines
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy
        if(ierr /= 0) then
          call io_error('Error: Problem opening input file SHIN')
          call io_error(msg)
        endif
      !I convert all tabulation characters to spaces
      ipos = index(dummy,TABCHAR)
      do while (ipos /= 0)
        dummy(ipos:ipos) = ' '
        ipos = index(dummy,TABCHAR)
      end do
      !
      dummy=utility_lowercase(dummy) !将字符串中大写字母全部改为小写
      dummy=trim(adjustl(dummy))     !
      if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
      if(len(trim(dummy)) == 0 ) cycle
      if(index(dummy,'=') <=1 )  cycle  !当该行中没有‘=’ 或‘=’前没有内容则跳过该行
      line_counter=line_counter+1
      
      !去除有效行信息中的注释部分，注释可以采用 ！或者 #
      in1=index(dummy,'!')
      in2=index(dummy,'#')
      if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
      !不存在'!'与'#'
      if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
      if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
      if(in2> 0 .and. in1>0 )  in_data(line_counter)=dummy(:min(in1,in2)-1)
      
      !如果输入参数为字符串，则给字符串加上"*"
      !itmp = index(dummy,'=')
      !ctmp = dummy(:itmp-1)
      !ctmp1= dummy(itmp+1:)
      !ctmp1= trim(adjustl(ctmp1))
      !itmp = ichar(ctmp1(1:1))
      !if(itmp>=ia .and. itmp<=iz) then
      !  dummy = ctmp//"='"//ctmp1//"'"
      !endif
      
    end do
    !得到包含有效信息的行数line_counter,和相应的数据in_data(line_counter)

    close(in_unit)

  end subroutine treat_inputfile 
  
  subroutine read_namelist()
    !use parameters
    use io
    implicit none
    integer::incar_unit,i
    
    !!write input file to namelist input file
    incar_unit = io_file_unit()
    open(unit=incar_unit,status='SCRATCH',iostat=ierr,iomsg=msg)
    if(ierr > 0 ) then
      call io_error('Error: Problem reading SCRATCH input namelist file')
      call io_error(msg)
    elseif(ierr == 0 )then    
      write(incar_unit,*)"&shinput" 
      do i=1,line_counter
        write(incar_unit,*) trim(adjustl(in_data(i)))
      enddo
      !write(incar_unit,"(A1)") "/"
      write(incar_unit,*) "/"
    endif
    rewind(incar_unit)
    
    !   set default values for variables in namelist
    Lrunsh        = .True.
    shtype        = "elec"
    dimention     = 1
    na1           = 1
    na2           = 1
    na3           = 1
    nk1           = 1
    nk2           = 1
    nk3           = 1
    nq1           = 1
    nq2           = 1
    nq3           = 1
    representation= "wfstat"
    nshiftstep    = 3
    dtadq         = 0.1
    pj_adjust     = 0.04
    temp          = 300
    gamma         = 20
    dt            = 0.1
    nstep         = 100
    nsnap         = 100
    naver         = 100
    initnmstat    = "quantum"
    Lrandomsita   = .TRUE.
    L_hotphonon   = .TRUE.
    hot_mode      = 1
    hot_scal      = 1.5
    Num_occupied  = 2
    initehstat    = "EN"
    initek        = 37
    initeb        = 9
    inithk        = 37
    inithb        = 10
    initeWF       = 1
    inithWF       = 1
    initeEN       = 0.75
    inithEN       = 0.75
    initeES       = 15
    inithES       = 16
    L_exciton     = .FALSE.
    epsr          = 2.4
    MSH           = "SC-FSSH"
    Ldecoherece   = .FALSE.
    Tdecoherence  = "enbased"
    Lfeedback     = .TRUE.
    mkl_threads   = 1 
    !!!end
    
    write(stdout,"(/,1X,A67)")   repeat("=",67)
    write(stdout,"(1X,10X,A)") "The namelist file as follows"
    write(stdout,"(1X,A67)")   repeat("=",67)
    do i=1,line_counter+2
      read(incar_unit,"(A80)") ctmp
      write(stdout,"(A80)") ctmp
    enddo
    rewind(incar_unit)
    read(UNIT=incar_unit,nml=shinput,iostat=ierr,iomsg=msg)
    if(ierr /= 0) then
      call io_error('Error: Problem reading namelist file SHIN')
      call io_error(msg)
    endif  
    close(incar_unit)
    write(stdout,"(/,1X,A)") "Read parameter Successful!" 
    write(stdout,*)   repeat("=",67)
    
  end subroutine read_namelist
  
  
end module readinput