! subgrider :: Julio Caineta, 2010
! sg_gridutils:: modulo de dependencias para manipulacao de grids (sg_griders)

module sg_gridutils

use qsort_c
implicit none

public::checkfile,del_sorted_nd,header_skip,itochar,media,percentil,percentil_sims,percentil_updt,rtochar,variancia,gridcoord,front_trim,get_3dcoord

private

! tipo de dados grid
type, public::grid
    integer::dx,dy,dz
    real::nd
    real,allocatable::val(:)
end type grid

contains

! verifica se um ficheiro existe
subroutine checkfile(file)
character(len=*),intent(in)::file
logical::file_exists
inquire(file=file,exist=file_exists)
if (.not.file_exists) then
    print *,"ficheiro ",trim(file)," nao encontrado."
    stop
end if
end subroutine checkfile

! funcao para converter inteiro para caracter(256)
function itochar(i) result(x)
integer,intent(in)::i
character(len=256)::x
write (x,*) i
end function itochar

! funcao para converter real para caracter(256)
function rtochar(r) result(x)
real,intent(in)::r
character(len=256)::x
write (x,*) r
end function rtochar

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

! calcula a media de um vector
function media(data) result(med)
real,intent(in),dimension(:)::data
real::med
med=sum(data)/size(data)
end function

! calcula a variancia de um vector
function variancia(data) result(var)
real,intent(in),dimension(:)::data
real::var
var=sum((data-media(data))**2)/size(data)
!var=media(data**2)-media(data)**2
end function variancia

! calcula o percentil p de uma grid
subroutine percentil(res,p,q,timer)
type(grid),intent(in)::res
real,intent(inout)::p
real,intent(out)::q
real,intent(out)::timer
type(grid)::aux,aux2
integer::i
real::start,finish,t
call cpu_time(start)
aux=res
if (p>1) p=p/100
call QsortC(aux%val)
call del_sorted_nd(aux%val,aux2%val,aux%nd,t)
deallocate(aux%val)
i=size(aux2%val)*p
q=aux2%val(i)
call cpu_time(finish)
timer=finish-start
end subroutine percentil

! percentil para matriz de simulacoes
subroutine percentil_sims(sims,p,q,id,header,output,timer)
real,dimension(:,:),intent(in)::sims
real,intent(in), dimension(2)::p, q
integer,intent(in)::id
logical,intent(in)::header
character(len=256)::output
real,intent(out)::timer
type(grid)::aux
integer::i
real,allocatable,dimension(:)::ps1,ps2,qs
real::start,finish,t,p1,p2,psim,q1,q2,pp(2)
call cpu_time(start)
open (id,file=output)
if (header) call  header_skip(id)
allocate(ps1(size(sims,1)))
allocate(ps2(size(sims,1)))
allocate(qs(size(sims,1)))
allocate(aux%val(size(sims,2)))
pp=p
do i=1,size(sims,1)
    aux%val=sims(i,:)
    call percentil(aux,pp(1),q1,t)
    call percentil(aux,pp(2),q2,t)
    call percentil_updt(aux,q(1),p1,t)
    call percentil_updt(aux,q(2),p2,t)
    psim=p2-p1
    ps1(i)=q1
    ps2(i)=q2
    qs(i)=psim
    write (id,*) q1,q2,psim
end do
write (id,*) "media ",media(ps1),media(ps2),media(qs)
write (id,*) "variancia ",variancia(ps1),variancia(ps2),variancia(qs)
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine percentil_sims

! calcula o percentil correspondente a um determinado valor
subroutine percentil_updt(res,q,p,timer)
type(grid),intent(in)::res
real,intent(in)::q
real,intent(out)::p
real,intent(out)::timer
type(grid)::aux,aux2
real::i,start,finish,t
call cpu_time(start)
aux=res
call QsortC(aux%val)
call del_sorted_nd(aux%val,aux2%val,res%nd,t)
deallocate(aux%val)
do i=1,size(aux2%val)
    if (aux2%val(i)>=q) exit
end do
p=i/size(aux2%val)
call cpu_time(finish)
timer=finish-start
end subroutine percentil_updt

! elimina os no-data's num vector ordenado; auxiliar para percentil() e percentil_updt()
subroutine del_sorted_nd(data,aux,nd,timer)
real,dimension(:),intent(in)::data
real,allocatable,dimension(:),intent(out)::aux
real,intent(in)::nd
real,intent(out)::timer
integer::i,a
real::start,finish
call cpu_time(start)
do i=1,size(data)
    if (data(i) .ne. nd) exit
end do
a=size(data)
allocate (aux(size(data)-i+1))
aux=data(i:a)
call cpu_time(finish)
timer=finish-start
end subroutine del_sorted_nd

! recebe o numero de uma linha de um grid e devolve as coordenadas associadas
! deprecated version
function gridcoord(xi,yi,zi,bl_dim,dx,dy,dz,l) result(coord)
real, intent(in):: xi,yi,zi,bl_dim(3)
integer, intent(in)::l,dx,dy,dz
real::coord(3)
real(kind=8)::A,B,C
C=bl_dim(1)/real((dx*dy*dz))
A=floor(l/real((dx*dy))-C)
B=l-dx*dy*A
coord(1)=xi+bl_dim(1)*(l-dx*floor(l/real(dx)-C))
coord(2)=yi+bl_dim(2)*(floor(B/real(dx)-C)+1)
coord(3)=zi+bl_dim(3)*(A+1)
end function gridcoord

function get_3dcoord(xi,yi,zi,bl_dim,dx,dy,dz,l) result(coord)
integer,intent(in)::dx,dy,dz,l !ind,xlen,xylen
real,intent(in)::xi,yi,zi,bl_dim(3)
real::x,y,z,coord(3)
integer::t
z=(l-1)/(dx*dy)+1
t=(mod(l-1,(dx*dy)))
y=t/dx+1
x=mod(l-1,dx)+1
coord(1)=xi+bl_dim(1)*x
coord(2)=yi+bl_dim(2)*y
coord(3)=zi+bl_dim(3)*z
end function get_3dcoord

! devolve o numero da primeira entrada de uma string diferente de espaco
function front_trim(str) result(ftlen)
character(len=*), intent(in)::str
integer::ftlen,i
do i=1,len_trim(str)
    if (str(i:i).ne.' ') then
        ftlen=i
        return
    end if
end do
end function front_trim

end module sg_gridutils
