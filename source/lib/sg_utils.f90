! subgrider :: Julio Caineta, 2010
! sg_utils :: modulo de dependencias do interface (sg_main)

module sg_utils

use sg_griders
use sg_gridutils
implicit none

public::abre_pocos,abre_sims,header_ask,novo,rand_int,rand_well,tempo

private

contains

! pergunta se o ficheiro a ler/escrever tem cabecario
subroutine header_ask(tipo,header,nvar)
character(len=1),intent(in)::tipo
logical,intent(out)::header
integer,intent(out)::nvar
character(len=1)::resp
do
    if (tipo=="r") print *,"ficheiro de input tem cabecario? (s/n)"
    if (tipo=="w") print *,"ficheiro de output com cabecario? (s/n)"
    read *,resp
    if (resp=="s" .or. resp=="S") then
        header=.TRUE.
        print *,"numero de variaveis"
        read *,nvar
        exit
    elseif (resp=="n" .or. resp=="N") then
        header=.FALSE.
        exit
    else
        print *,"erro - resposta invalida"
    end if
end do
end subroutine header_ask

! carrega o ficheiro com as coordenadas dos pocos
subroutine abre_pocos(pocos,px,py,timer)
character(len=256),intent(in)::pocos
integer,allocatable,intent(out)::px(:),py(:)
integer::i,n
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (19,file=pocos,action='read')
read (19,*) n
allocate(px(n),py(n))
do i=1,n
    read (19,*) px(i),py(i)
end do
close(19)
call cpu_time(finish)
timer=finish-start
end subroutine abre_pocos

! cria um novo ficheiro
subroutine novo(header,nvar,id,newfile,batch,timer)
logical,intent(in)::header
integer,intent(in)::nvar,batch
integer,intent(inout)::id
character(len=256),intent(out)::newfile
character(len=50)::nv,novonome
real,intent(out)::timer
real::start,finish
integer::i
call cpu_time(start)
id=id+1
!if (batch>=0) then
!    newfile='wells.prn'
    !newfile=ficheiro(1:len_trim(ficheiro)-4)//'_wells.prn'
!else
    print *,"nome do ficheiro"
    read *,newfile
!end if
open (id,file=newfile,action='write')
if (header) then
    if (batch>=0) then
        novonome=newfile(1:len_trim(newfile)-4)
        write(id,*) trim(novonome)
        write (id,*) trim(itochar(nvar))
        nv='var'
        do i=1,nvar
            write (id,*) trim(nv)
        end do
    else
        print *,"nome do conjunto de dados"
        read *,novonome
        write(id,*) trim(novonome)
        write (id,*) trim(itochar(nvar))
        print *,"nomes das variaveis (espacados com enter)"
        do i=1,nvar
            read *,nv
            write (id,*) trim(nv)
        end do
    end if
end if
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine novo

! carrega ficheiros de simulacoes para uma matriz
subroutine abre_sims(nr_sims,sims,dg,base,ext,header,timer)
integer,intent(in)::nr_sims(2),dg(3)
logical,intent(in)::header
real,allocatable,intent(out)::sims(:,:)
character(len=256),intent(in)::base
character(len=4),intent(in)::ext
character(len=256)::sim
character(len=4)::simN,format
type(grid)::gridsim
integer::i,id
real,intent(out)::timer
real::start,finish,t
call cpu_time(start)
allocate(sims(nr_sims(2)-nr_sims(1)+1,product(dg)))
id=20
do i=nr_sims(1),nr_sims(2)
    if (i<10) format="(I1)"
    if (i>=10 .and. i<100) format="(I2)"
    if (i>=100 .and. i<1000) format="(I3)"
    write (simN,format) i
    sim=trim(base)//trim(simN)//trim(ext)
    call checkfile(sim)
    call abre(sim,header,dg,-999.0,gridsim,id,t)
    sims(i,:)=gridsim%val
    id=id+1
end do
call cpu_time(finish)
timer=finish-start
end subroutine abre_sims

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

subroutine rand_int(a,b,ri)
integer,intent(in)::a,b
integer,intent(out)::ri
real::colheita
call random_number(colheita)
if (b<a) then
    ri=int(b+(a-b)*colheita+0.5)
else
    ri=int(a+(b-a)*colheita+0.5)
end if
end subroutine rand_int

subroutine rand_well(k,dg,px,py,wells,timer)
integer,intent(in)::k,dg(3)
integer,dimension(:),intent(in)::px,py
integer,allocatable,dimension(:,:,:),intent(out)::wells
real,intent(out)::timer
integer::fixo(2),i,j,maxp(2),minp(2),xx(2),yy(2),limx(2),limy(2),muv(2)
integer,allocatable,dimension(:,:)::difs
real::start,finish
call cpu_time(start)
! fixar um poco(x,y)
fixo(1)=px(1)
fixo(2)=py(1)
! diferencas(x,y) entre o fixo e todos
allocate(difs(size(px),2))
do i=1,size(px)
    difs(i,1)=fixo(1)-px(i)
    difs(i,2)=fixo(2)-py(i)
end do
! max(x,y) e min(x,y) entre todos
maxp(1)=maxval(px)
maxp(2)=maxval(py)
minp(1)=minval(px)
minp(2)=minval(py)
! quanto pode andar em x(direita,esquerda) e em y(cima,baixo)
xx(1)=dg(1)-maxp(1)
xx(2)=minp(1)-1
yy(1)=dg(2)-maxp(2)
yy(2)=minp(2)-1
! intervalo para gerar numeros aleatorios em x(min,max) e em y(min,max)
limx(1)=fixo(1)-xx(2)
limx(2)=fixo(1)+xx(1)
limy(1)=fixo(2)-yy(2)
limy(2)=fixo(2)+yy(1)
! gerar numeros aleatorios(x,y) e novas posicoes(x,y,k)
allocate(wells(size(px),2,k))
do i=1,k
    call rand_int(limx(1),limx(2),muv(1))
    call rand_int(limy(1),limy(2),muv(2))
    do j=1,size(px)
        wells(j,1,i)=muv(1)-difs(j,1)
        wells(j,2,i)=muv(2)-difs(j,2)
    end do
end do
call cpu_time(finish)
timer=finish-start
end subroutine rand_well

end module
