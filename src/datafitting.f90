module datafitting
  use kinds
  implicit none
  
  contains
  subroutine linefitting(nstep,x,y,slope,y_int,r)
    implicit none
    integer,intent(in)       :: nstep
    real(kind=dp),intent(in) :: x(nstep),y(nstep)
    real(kind=dp),intent(out):: slope,y_int,r
    integer                  :: istep
    real(kind=dp)            :: sum_x,sum_y,sum_x2,sum_y2,sum_xy
    real(kind=dp)            :: x_bar,y_bar
    
    !initial
    slope = 0.0
    y_int = 0.0
    r     = 0.0
    sum_x = 0.0
    sum_x2= 0.0
    sum_y = 0.0
    sum_y2= 0.0
    sum_xy= 0.0
    x_bar = 0.0
    y_bar = 0.0
    
    sum_x = SUM(x)
    sum_x2= SUM(x**2)
    sum_y = SUM(y)
    sum_y2= SUM(y**2)
    sum_xy= SUM(x*y)
    x_bar = sum_x/real(nstep)
    y_bar = sum_y/real(nstep)
    
    slope = (sum_xy - sum_x * y_bar)/(sum_x2 - sum_x * x_bar)
    y_int = y_bar-slope * x_bar
    r= ( nstep*sum_xy-sum_x*sum_y )/sqrt( (nstep*sum_x2-sum_x**2)*(nstep*sum_y2-sum_y**2) )
    
  end subroutine linefitting

  subroutine writeHepfit(fitname,num_wann,nrpts,ndegen,irvec,Ham_r_slope,Ham_r_yint,Ham_r_r)
    use io
    use constants,only : maxlen  
    implicit none
    character(len=maxlen),intent(in) :: fitname
    integer,intent(in)               :: num_wann
    integer,intent(in)               :: nrpts
    integer,intent(in)               :: ndegen(nrpts)
    integer,intent(in)               :: irvec(3,nrpts)
    real(kind=dp),intent(in)         :: Ham_r_slope(num_wann,num_wann,nrpts),&
                                        Ham_r_yint(num_wann,num_wann,nrpts),&
                                        Ham_r_r(num_wann,num_wann,nrpts)
    integer :: fitunit
    integer :: irpt
    character(len=9)  :: cdate,ctime 
    character(len=maxlen) :: header
    integer :: i,j
    logical :: lexist
    fitunit = io_file_unit()
    call open_file(fitname,fitunit)
    call io_date(cdate,ctime)
    header = 'written on '//cdate//' at '//ctime
    write(fitunit,*) header
    write(fitunit,*) num_wann
    write(fitunit,*) nrpts
    write(fitunit,'(15I5)') (ndegen(irpt),irpt=1,nrpts)
    write(fitunit,"(1X,5A5,3A12)") " R1 ", " R2 "," R3 "," n "," m "," slope "," fitH0 "," fitcorr "
    do irpt=1,nrpts
      do i=1,num_wann
        do j=1,num_wann
          write(fitunit,'(5I5,3F12.6)') irvec(:,irpt),j,i,Ham_r_slope(j,i,irpt),&
                                        Ham_r_yint(j,i,irpt),Ham_r_r(j,i,irpt)
        enddo
      enddo
    enddo
    call close_file(fitname,fitunit)
    
  end subroutine writeHepfit
  
end module datafitting
    
    