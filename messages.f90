module messages
use types, only : rp => rprec
implicit none
$if ($MPI)
  ! use mpi
  include "mpif.h"
$endif
! implicit none

save
private

public :: enter_sub, exit_sub, error, warn,mesg,output_intm,output_slice,output_point,output_3point,output_intmz
public :: msg, n_l
public :: output_slice_reduce,output_x

integer, parameter :: n_blanks = 32
integer, parameter :: n_msg = 1024

character (n_blanks), parameter :: blanks = repeat (' ', n_blanks)
character (2), parameter :: n_l = achar (10) // ' '
                            ! carriage return and a space

character (n_msg) :: msg  ! experiment
character (64) :: fmt

integer, parameter :: lun = 6  ! system dependent

integer :: call_level = 0
$if ($MPI)
  integer :: ierr
$endif

interface error
  module procedure error_a, error_ai, error_al, error_aia, error_ai_array,  &
                   error_aiar, error_ar
end interface

interface mesg
  module procedure message_a, message_ai, message_aiai, message_aiar,  &
                   message_al, message_ar, message_aii, message_air,   &
                   message_ai_array, message_aiai_array,               &
                   message_ar_array, message_aiar_array
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine enter_sub (name)

character (*), intent (in) :: name

integer :: n

!---------------------------------------------------------------------

call_level = call_level + 1

n = min (n_blanks, call_level-1)

write (lun, '(1x,a)') blanks(1:n) // name // ': started'


end subroutine enter_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exit_sub (name)

character (*), intent (in) :: name

integer :: n

!---------------------------------------------------------------------

n = min (n_blanks, call_level-1)

write (lun, '(1x,a)') blanks(1:n) // name // ': done'

call_level = call_level - 1

end subroutine exit_sub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_a (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') name // ': ' // trim (msg)

end subroutine message_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ai (name, msg, i)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i

write (lun, '(1x,a,1x,i0)') name // ': ' // trim (msg), i

end subroutine message_ai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiai (name, msg1, i1, msg2, i2)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i1, i2

write (lun, '(2(1x,a,1x,i0))') name // ': ' // trim (msg1), i1,  &
                                               trim (msg2), i2

end subroutine message_aiai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiar (name, msg1, i, msg2, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r

fmt = '(1x,a,1x,i0,1x,a,1x,es11.4)'

write (lun, fmt) name // ': ' // trim (msg1), i, trim (msg2), r

end subroutine message_aiar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_al (name, msg, l)

character (*), intent (in) :: name
character (*), intent (in) :: msg
logical, intent (in) :: l

write (lun, '(1x,a,1x,l1)') name // ': ' // trim (msg), l

end subroutine message_al

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aii (name, msg, i, j)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i, j

fmt = '(1x,a,1x,i0,",",1x,i0)'
                    !--comma to separate output
write (lun, fmt) name // ': ' // trim (msg), i, j

end subroutine message_aii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_air (name, msg, i, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i
real (rp), intent (in) :: r

fmt = '(1x,a,1x,i0,",",1x,es11.4)'
                    !--comma to separate output 
write (lun, fmt) name // ': ' // trim (msg), i, r

end subroutine message_air

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ai_array (name, msg, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,', n, '(1x,i0))'
write (lun, fmt) name // ': ' // trim (msg), i_arr

end subroutine message_ai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiai_array (name, msg1, i, msg2, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,i0,a,', n, '(1x,i0))'
write (lun, fmt) name // ': ' // trim (msg1), i, trim(msg2), i_arr

end subroutine message_aiai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ar (name, msg, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r

write (lun, '(1x,a,1x,es11.4)') name // ': ' // trim (msg), r

end subroutine message_ar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_ar_array (name, msg, r_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (r_arr)
write (fmt, *) '(1x,a,', n, '(1x,es11.4))'
write (lun, fmt) name // ': ' // trim (msg), r_arr

end subroutine message_ar_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine message_aiar_array (name, msg1, i, msg2, r_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (r_arr)
write (fmt, *) '(1x,a,i0,a,', n, '(1x,es11.4))'
write (lun, fmt) name // ': ' // trim (msg1), i, trim (msg2), r_arr

end subroutine message_aiar_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine warn (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') '*****WARNING*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a)') trim (msg)
write (lun, '(1x,a)') '*****************'

end subroutine warn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_a (name, msg)

character (*), intent (in) :: name
character (*), intent (in) :: msg

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a)') trim (msg)
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ai (name, msg, i)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0)') trim (msg), i
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_ai

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ai_array (name, msg, i_arr)

character (*), intent (in) :: name
character (*), intent (in) :: msg
integer, intent (in) :: i_arr(:)

integer :: n

!---------------------------------------------------------------------

n = size (i_arr)
write (fmt, *) '(1x,a,', n, '(1x,i0))'

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, fmt) trim (msg), i_arr
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_ai_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_aia (name, msg1, i, msg2)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0,1x,a)') trim (msg1), i, trim (msg2)
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_aia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_aiar (name, msg1, i, msg2, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg1, msg2
integer, intent (in) :: i
real (rp), intent (in) :: r

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,i0,1x,a,1x,es11.4)') trim (msg1), i, trim (msg2), r
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_aiar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_al (name, msg, l)

character (*), intent (in) :: name
character (*), intent (in) :: msg
logical, intent (in) :: l

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,l1)') trim (msg), l
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_al

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error_ar (name, msg, r)

character (*), intent (in) :: name
character (*), intent (in) :: msg
real (rp), intent (in) :: r

write (lun, '(1x,a)') '*****ERROR*****'
write (lun, '(1x,a)') 'In ' // name // ':'
write (lun, '(1x,a,1x,es11.4)') trim (msg), r
write (lun, '(1x,a)') '***************'
write (lun, '(1x,a)') 'Program aborted'

$if ($MPI)
  !call mpi_finalize (ierr)
$endif

stop

end subroutine error_ar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_intm (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,$lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,i1.1,a)') filename,coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,i2.1,a)') filename,coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
 do k = 1, nz-1
  do j = 1, ny
 write(filenum,5168) (x(i,j,k),i=1,nx)
  enddo
 enddo
 endif
 close(filenum)

end if
5168     format(1400(E14.5))
end subroutine output_intm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_intmz (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension($lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,i1.1,a)') filename,coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,i2.1,a)') filename,coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
 do k = 1, nz-1
 write(filenum,5168) (x(k))
 enddo
 endif
 close(filenum)

end if
5168     format(1400(E14.5))
end subroutine output_intmz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_slice_reduce (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,$lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,i1.1,a)') filename,coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,i2.1,a)') filename,coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
  do j = 1, ny/3
 write(filenum,5168) (x((i-1)*3+1,(j-1)*3+1,nz-1),i=1,nx/3)
  enddo
 endif
 close(filenum)

end if
5168     format(1400(E14.5))
end subroutine output_slice_reduce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_x (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,$lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,i1.1,a)') filename,coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,i2.1,a)') filename,coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
  do j = ny/2, ny/2
 write(filenum,5168) (x(i,j,nz-1),i=1,nx)
  enddo
 endif
 close(filenum)

end if
5168     format(1400(E14.5))
end subroutine output_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_slice (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,1)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,i1.1,a)') filename,coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,i2.1,a)') filename,coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
  do j = 1, ny
 write(filenum,5168) (x(i,j,nz-1),i=1,nx)
  enddo
 endif
 close(filenum)

end if
5168     format(1400(E14.5))
end subroutine output_slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_3point (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1,ffname2,ffname3,ffname4,ffname5,ffname6

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,$lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,a,i1.1,a)') filename,'c',coord,'.dat'
  write(ffname2,'(a,a,i1.1,a)') filename,'d',coord,'.dat'
  write(ffname3,'(a,a,i1.1,a)') filename,'e',coord,'.dat'
  write(ffname4,'(a,a,i1.1,a)') filename,'f',coord,'.dat'
  write(ffname5,'(a,a,i1.1,a)') filename,'g',coord,'.dat'
  write(ffname6,'(a,a,i1.1,a)') filename,'h',coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,a,i2.1,a)') filename,'c',coord,'.dat'
  write(ffname2,'(a,a,i2.1,a)') filename,'d',coord,'.dat'
  write(ffname3,'(a,a,i2.1,a)') filename,'e',coord,'.dat'
  write(ffname4,'(a,a,i2.1,a)') filename,'f',coord,'.dat'
  write(ffname5,'(a,a,i2.1,a)') filename,'g',coord,'.dat'
  write(ffname6,'(a,a,i2.1,a)') filename,'h',coord,'.dat'
  endif
 !endif
!c
 open(filenum, file=path//ffname1, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(12+(i-1)*32,4+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(10+(i-1)*32,4+(j-1)*16,k),i=1,6)!P4
 !write(filenum,5168) (x(13+(i-1)*32,4+(j-1)*16,k),i=1,6)!P2
 !write(filenum,5168) (x(18+(i-1)*40,3+(j-1)*20,k),i=1,5)!F1
 write(filenum,5168) (x(3+(i-1)*40,18+(j-1)*20,k),i=1,5)!F6
     enddo
  enddo
 close(filenum)
!d
 open(filenum, file=path//ffname2, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(28+(i-1)*32,12+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(26+(i-1)*32,12+(j-1)*16,k),i=1,6) !P4
 !write(filenum,5168) (x(29+(i-1)*32,12+(j-1)*16,k),i=1,6) !P2
 !write(filenum,5168) (x(38+(i-1)*40,13+(j-1)*20,k),i=1,5) !F1
 write(filenum,5168) (x(13+(i-1)*40,38+(j-1)*20,k),i=1,5)  !F6
     enddo
  enddo
 close(filenum)
!e
 open(filenum, file=path//ffname3, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(12+(i-1)*32,12+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(10+(i-1)*32,12+(j-1)*16,k),i=1,6) !P4
 !write(filenum,5168) (x(13+(i-1)*32,12+(j-1)*16,k),i=1,6) !P2
 !write(filenum,5168) (x(18+(i-1)*40,13+(j-1)*20,k),i=1,5) !F1
 write(filenum,5168) (x(13+(i-1)*40,18+(j-1)*20,k),i=1,5)  !F6
     enddo
  enddo
 close(filenum)
!f
 open(filenum, file=path//ffname4, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(28+(i-1)*32,4+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(28+(i-1)*32,4+(j-1)*16,k),i=1,6) !P4
 !write(filenum,5168) (x(29+(i-1)*32,4+(j-1)*16,k),i=1,6) !P2
 !write(filenum,5168) (x(38+(i-1)*40,3+(j-1)*20,k),i=1,5) !F1
 write(filenum,5168) (x(3+(i-1)*40,38+(j-1)*20,k),i=1,5) !F6
     enddo
  enddo
 close(filenum)
!g
 open(filenum, file=path//ffname5, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(4+(i-1)*32,12+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(2+(i-1)*32,12+(j-1)*16,k),i=1,6) !P4
 !write(filenum,5168) (x(5+(i-1)*32,12+(j-1)*16,k),i=1,6) !P2
 !write(filenum,5168) (x(8+(i-1)*40,13+(j-1)*20,k),i=1,5) !F1
 write(filenum,5168) (x(13+(i-1)*40,8+(j-1)*20,k),i=1,5)  !F6
     enddo
  enddo
 close(filenum)
!h
 open(filenum, file=path//ffname6, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(20+(i-1)*32,4+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(18+(i-1)*32,4+(j-1)*16,k),i=1,6) !P4
 !write(filenum,5168) (x(21+(i-1)*32,4+(j-1)*16,k),i=1,6) !P2
 !write(filenum,5168) (x(28+(i-1)*40,3+(j-1)*20,k),i=1,5) !F1
 write(filenum,5168) (x(3+(i-1)*40,28+(j-1)*20,k),i=1,5) !F6
     enddo
  enddo
 close(filenum)
 endif
end if
5168     format(1400(E14.5))
end subroutine output_3point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_point (filename,filenum,coordnum,outtime,x)
use param
character (*), intent (in) :: filename
integer, intent (in) :: filenum
integer, intent (in) :: coordnum,outtime
character (len=1024) :: ffname1,ffname2

integer :: i,j,k
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rp),dimension(ld,ny,$lbz:nz)::x


if ( mod(jt,outtime) ==0) then
 if (coord == coordnum ) then
  if (coord .lt. 10) then
  write(ffname1,'(a,a,i1.1,a)') filename,'a',coord,'.dat'
  write(ffname2,'(a,a,i1.1,a)') filename,'b',coord,'.dat'
  elseif  ((coord .lt. 100) .and. (coord .ge. 10)) then
  write(ffname1,'(a,a,i2.1,a)') filename,'a',coord,'.dat'
  write(ffname2,'(a,a,i2.1,a)') filename,'b',coord,'.dat'
  endif
 !endif
 open(filenum, file=path//ffname1, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(4+(i-1)*32,4+(j-1)*16,k),i=1,6) !RSP5
 !write(filenum,5168) (x(2+(i-1)*32,4+(j-1)*16,k),i=1,6) !RSP4
 !write(filenum,5168) (x(5+(i-1)*32,4+(j-1)*16,k),i=1,6) !RSP2
 !write(filenum,5168) (x(8+(i-1)*40,3+(j-1)*20,k),i=1,5) !RSF1
 write(filenum,5168) (x(3+(i-1)*40,8+(j-1)*20,k),i=1,5) !RSF6
     enddo
  enddo
 close(filenum)

 open(filenum, file=path//ffname2, status='unknown',position='append')
  do k = 1, nz-1
     do j = 1,5
 !write(filenum,5168) (x(20+(i-1)*32,12+(j-1)*16,k),i=1,6) !P5
 !write(filenum,5168) (x(18+(i-1)*32,12+(j-1)*16,k),i=1,6)!P4
 !write(filenum,5168) (x(21+(i-1)*32,12+(j-1)*16,k),i=1,6)!P2
 !write(filenum,5168) (x(28+(i-1)*40,13+(j-1)*20,k),i=1,5)!F1
 write(filenum,5168) (x(13+(i-1)*40,28+(j-1)*20,k),i=1,5)!F6
     enddo
  enddo
 close(filenum)
endif
end if
5168     format(1400(E14.5))
end subroutine output_point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module messages
