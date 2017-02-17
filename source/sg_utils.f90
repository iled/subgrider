! subgrider :: Julio Caineta, 2010
! sg_utils :: modulo de dependencias do interface (sg_main)

module sg_utils

use sg_griders
use sg_gridutils
implicit none

public::abre_pocos,abre_sims,header_ask,novo

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
if (batch>=0) then
    newfile='wells.prn'
    !newfile=ficheiro(1:len_trim(ficheiro)-4)//'_wells.prn'
else
    print *,"nome do ficheiro"
    read *,newfile
end if
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

end module
