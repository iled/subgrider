program teste_arrays

implicit none

character(len=256)::ficheiro
logical::header
integer::dg(3),id,i
real::nd,timer,start,finish


type::grid
    integer::dx,dy,dz
    real::nd
    real,allocatable::val(:)
end type grid

type(grid)::grida

ficheiro='shah.prn'
header=.TRUE.
dg=(/ 334, 134, 93 /)
nd=-999.0
id=10

print *,"a ler ",trim(ficheiro),"..."
call abre(ficheiro,header,dg,nd,grida,id,timer)
print *,tempo(timer)
print *,"a escrever sequencial"
call cpu_time(start)
open (id+1,file='teste_array_seq.txt',action='write')
do i=1,product(dg)
write(id+1,*) grida%val(i)
end do
call cpu_time(finish)
print *,tempo(finish-start)
print *,"a escrever directo"
call cpu_time(start)
open (id+2,file='teste_array_dir.txt',action='write')
write (id+2,*) grida%val(:)
call cpu_time(finish)
print *,tempo(finish-start)
stop

contains

subroutine abre(ficheiro,header,dg,nd,grida,id,timer)
character(len=256), intent(in) :: ficheiro
integer, intent(in) :: dg(3),id
logical,intent(in)::header
type(grid),intent(out)::grida
real,intent(in)::nd
real,intent(out)::timer
real::start,finish
integer::i
call cpu_time(start)
grida%dx=dg(1)
grida%dy=dg(2)
grida%dz=dg(3)
grida%nd=nd
allocate(grida%val(product(dg)))
open (id, file=ficheiro, action='read')
if (header) call header_skip(id)
do i=1,size(grida%val)
    read(id,*) grida%val(i)
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine abre

! devolve o tempo, em string, gasto pela ultima subrotina executada
function tempo(timer) result(time)
real,intent(in)::timer
real::te
!character(10),allocatable,dimension(:)::time
character(32)::time
character(len=6)::format,uni
if (timer<1) then
    te=timer*1000
    format="(f7.3)"
    uni=" ms"
elseif (timer>=1 .and. timer<60) then
    te=timer
    format="(f6.3)"
    uni=" s"
elseif (timer>=60 .and. timer<3600) then
    te=timer/60
    format="(f6.3)"
    uni=" min"
else
    te=timer/3600
    format="(f7.3)"
    uni=" h"
end if
write (time,format) te
time=time(1:len_trim(time))//trim(uni)
end function tempo

! salta o cabecario de um ficheiro
subroutine header_skip(id)
integer,intent(in)::id
integer::nvar,i
read (id,*)
read (id,*) nvar
do i=1,nvar
    read (id,*)
end do
end subroutine header_skip

end program teste_arrays
