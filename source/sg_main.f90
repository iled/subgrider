program subgrider

use qsort_c_module
implicit none

! declaracao de variaveis
integer::i,op,nvar,fid,batch,nr_sims,nr_p
integer,dimension(3)::dg,pa,pb,dg_2d
character(len=256)::ficheiro,output,pocos,base
character(len=1)::load,mp_p
logical::mp,header
real,allocatable::sims(:,:)
integer,allocatable::px(:),py(:)
real::P_x,P1,timer,nd,p,q
! tipo de dados grid
type::grid
    integer::dx,dy,dz
    real::nd
    real,allocatable::val(:)
end type grid
type(grid)::res
fid=10
! input inicial
print *,"----|| s u b g r i d e r  v0.19 ||----"
print *,""
print *,"carregar o ficheiro anterior? (s/n)"
read *,load
if (load=="s" .or. load=="S") then
	call checkfile("subgrider.cfg")
    open (9,file='subgrider.cfg',action='read')
    read (9,*) ficheiro,header
    read (9,*) dg,nd
elseif (load=="n" .or. load=="N") then
    print *,"ficheiro a carregar (caminho/nome)"
    read *,ficheiro
    call checkfile(ficheiro)
    call header_ask("r",header,nvar)
	print *,"dimensoes da grid: x y z"
	read (*,*) dg
	print *,"valor correspondente a 'no data'"
	read *,nd
	open (9,file='subgrider.cfg',action='write')
	write (9,*) ficheiro,header
	write (9,*) dg,nd
	close (9)
else
    print *,"erro - nao sei ler isso"
    stop
end if
print *,"a carregar o ficheiro ",trim(ficheiro),"..."
call abre(ficheiro,header,dg,nd,res,fid,timer)
print *,"ficheiro carregado em ",tempo(timer)

!ciclo principal
do
batch=-1
print *,""
print *,"escolher uma opcao"
print *,"1 ----- poco vertical "
print *,"  1.1 - a partir das coordenadas (x,y)"
print *,"2 - obter n pocos verticais a partir de um ficheiro com coordenadas (x,y)"
print *,"3 - obter uma grid a partir um volume (cubo, paralelipipedo) da grid inicial"
print *,"4 - criar uma mascara (data=1, no data=0)"
print *,"5 - calcular verosimilhanca"
print *,"6 - calcular percentil"
print *,"7 - corrigir quantil"
print *,"0 - sairs"
read *,op
if (op==1) then
    print *,"opcao 1"
    print *,"coordenadas do poco (x,y)"
    read *,pa(1),pa(2)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a furar o poco..."
    call pp(pa(1),pa(2),res,fid,output,header,batch,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==2) then
    print *,"opcao 2"
    print *,"ficheiro com coordenadas (caminho/nome)"
    read *,pocos
	call checkfile(pocos)
    print *,"a ler ficheiro ",trim(pocos),"..."
    call abre_pocos(pocos,px,py,timer)
    print *,trim(pocos)," carregado em",tempo(timer)
    call header_ask("w",header,nvar)
    batch=0
	call novo(header,nvar,fid,output,batch,timer)
	print *,"a furar os pocos..."
    do i=1,size(px)
    	call pp(px(i),py(i),res,fid,output,header,batch,timer)
    end do
    print *,"operacao concluida em ",tempo(timer) !tempo devolvido é do último poço
    batch=-1
elseif (op==3) then
    print *,"opcao 3"
    print *,"coordenadas x (min,max)"
    read *,pa(1),pb(1)
    print *,"coordenadas y (min,max)"
    read *,pa(2),pb(2)
    print *,"coordenadas z (min,max)"
    read *,pa(3),pb(3)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a criar subgrid..."
    call subgrid(pa,pb,res,fid,output,header,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==4) then
	print *,"opcao 4"
	call novo(header,nvar,fid,output,batch,timer)
	print *,"a criar mascara..."
	call mask(1,0,header,res,fid,output,timer)
	print *,"operacao concluida em ",tempo(timer)
elseif (op==5) then
    print *,"opcao 5"
    print *,"numero de simulacoes realizadas (ficheiros *_simX.OUT)"
    read *,nr_sims
    base=ficheiro(1:len_trim(ficheiro)-4)//'_sim'
    call header_ask("r",header,nvar)
    print *,"a carregar grids simuladas..."
    call abre_sims(nr_sims,sims,dg,base,'.OUT',header,timer)
    print *,nr_sims,"grids carregadas em ",tempo(timer)
    print *,"produzir grid com mapa de probabilidades? (s/n)"
    read *,mp_p
    if (mp_p=="s" .or. mp_p=="S") then
        mp=.true.
    elseif (mp_p=="n" .or. mp_p=="N") then
        mp=.false.
    else
        print *,"erro - input invalido."
        stop
    end if
    if (mp) then
    	call header_ask("w",header,nvar)
    	call novo(header,nvar,fid,output,batch,timer)
	end if
	print *,"inserir valor x (P[X=x])"
	read *,P_x
    print *,"a calcular probabilidades..."
    call likely(res,sims,P_x,P1,mp,fid,output,header,timer)
    print *,"a probabilidade e' ",P1
    print *,"operacao concluida em ",tempo(timer)
elseif (op==6) then
	print *,"opcao 6"
	print *,"inserir percentil a calcular"
	read *,p
	call percentil(res,p,q,timer)
	print '("o percentil ",f4.2," e ",f)',p,q
	print *,"operacao concluida em ",tempo(timer)
elseif (op==7) then
	print *,"opcao 7"
	print *,"inserir valor"
	read *,q
	call percentil_updt(res,q,p,timer)
	print *,"o valor ",q," corresponde ao percentil ",p
	print *,"operacao concluida em ",tempo(timer)
elseif (op==0) then
    print *,"programa terminado"
    stop
end if
end do

! funcoes e subrotinas

contains

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

subroutine header_ask(tipo,header,nvar)
character(len=1),intent(in)::tipo
logical,intent(out)::header
integer,intent(out)::nvar
character(len=1)::resp
do
	if (tipo=="r") print *,"ficheiro de input tem cabecario? (s/n)"
	if (tipo=="w") print *,"ficheiro de output tem cabecario? (s/n)"
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

subroutine header_skip(id)
integer,intent(in)::id
integer::nvar
read (id,*)
read (id,*) nvar
do i=1,nvar
	read (id,*)
end do
end subroutine header_skip

subroutine checkfile(file)
character(len=*),intent(in)::file
logical::file_exists
inquire(file=file,exist=file_exists)
if (.not.file_exists) then
    print *,"ficheiro ",trim(file)," nao encontrado."
    stop
end if
end subroutine checkfile

! carrega o ficheiro de dados para um vector
subroutine abre(ficheiro,header,dg,nd,grida,id,timer)
character(len=256), intent(in) :: ficheiro
integer, intent(in) :: dg(3),id
logical,intent(in)::header
type(grid),intent(out)::grida
real,intent(in)::nd
real,intent(out)::timer
real::start,finish
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
print '("carregadas as coordenadas de ",I3," pocos.")',n
print *,"ficheiro ",trim(pocos)," nao encontrado."
stop
call cpu_time(finish)
timer=finish-start
end subroutine abre_pocos

! carrega ficheiros de simulacoes para uma matriz
subroutine abre_sims(nr_sims,sims,dg,base,ext,header,timer)
integer,intent(in)::nr_sims,dg(3)
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
allocate(sims(nr_sims,product(dg)))
id=20
do i=1,nr_sims
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

! cria um novo ficheiro
subroutine novo(header,nvar,id,newfile,batch,timer)
logical,intent(in)::header
integer,intent(in)::nvar,batch
integer,intent(inout)::id
character(len=256),intent(out)::newfile
character(len=50)::nv,novonome
real,intent(out)::timer
real::start,finish
call cpu_time(start)
id=id+1
if (batch>=0) then
	newfile=ficheiro(1:len_trim(ficheiro)-4)//'_wells.prn'
else
	print *,"nome do ficheiro"
	read *,newfile
end if
open (id,file=newfile,action='write')
if (header) then
	if (batch>=0) then
		novonome=newfile(1:len_trim(newfile)-4)
    	write(id,*) novonome
    	write (id,*) nvar
    	nv='var'
    	do i=1,nvar
        	write (id,*) nv
    	end do
	else
	    print *,"nome do conjunto de dados"
    	read *,novonome
    	write(id,*) novonome
    	write (id,*) nvar
    	print *,"nomes das variaveis (espacados com enter)"
    	do i=1,nvar
        	read *,nv
        	write (id,*) nv
    	end do
	end if
end if
close(id)
call cpu_time(finish)
print *,"ai"
timer=finish-start
end subroutine novo

! papa um poco a partir de uma grid e devolve um point set
subroutine pp(xp,yp,res,id,output,header,batch,timer)
integer,intent(in)::xp,yp,id
integer,intent(inout)::batch
type(grid),intent(in)::res
character(len=256),intent(in)::output
logical,intent(in)::header
real::start,finish
real,intent(out)::timer
integer::z,p,m
call cpu_time(start)
open (id,file=output)
if (header) call header_skip(id)
if (batch>=0) then
	do m=1,batch
		read (id,*)
	end do
	batch=batch+93
end if
do z=1,res%dz
	p=xp+res%dx*(yp-1)+res%dx*res%dy*(z-1)
	write(id,*) xp,yp,z,res%val(p)
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine pp

! cria uma grid a partir de outra existente
subroutine subgrid(pa,pb,res,id,output,header,timer)
integer,dimension(3),intent(in)::pa,pb
integer,intent(in)::id
type(grid),intent(in)::res
character(len=256),intent(in)::output
logical,intent(in)::header
integer::pi,pf,d,p(3),j,c,v,nvar
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (id,file=output)
if (header) call header_skip(id)
p=pa
c=res%dx*res%dy
j=pa(3)-1
do while (j<pb(3))
    v=j*c+res%dx*(pa(2)-1)+pa(1)
    do while (v<=(pb(2)-1)*res%dx+pa(1)+j*c)
        i=0
        do while (i<=(pb(1)-pa(1)))
            write (id,*) res%val(v+i)
            i=i+1
        end do
        v=v+res%dx
    end do
    j=j+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine subgrid

! mascara (troca data por m1 e no data por m2)
subroutine mask(m1,m2,header,res,id,output,timer)
integer,intent(in)::id,m1,m2
type(grid),intent(in)::res
logical,intent(in)::header
character(len=256),intent(in)::output
real,intent(out)::timer
integer::j,i
real::start,finish
call cpu_time(start)
open (id,file=output)
if (header) call header_skip(id)
j=1
do while (j<=size(res%val))
	if (res%val(j)/=res%nd) then
		write (id,"(I2)") m1
	else
		write (id,"(I2)") m2
	end if
	j=j+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine mask

! calcula a verosimilhanca entre a grid inicial e as simuladas
subroutine likely(res,sims,k,p,mp,id,output,header,timer)
type(grid),intent(in)::res
real,dimension(:,:),intent(in)::sims
real,intent(in)::k
real,intent(out)::p
real,intent(out)::timer
integer,intent(in)::id
logical,intent(in)::mp,header
character(len=256),intent(in)::output
integer::j,l,a,b
real::start,finish
call cpu_time(start)
if (mp) then
	open(id,file=output)
	if (header) call header_skip(id)
end if
a=0
p=0
! open(33,file='likely.dbg',action='write') ! dd
! write (33,*) size(res%val), size(sims,1) ! dd
do j=1,size(res%val)
    if (res%val(j)>=k) then
        a=a+1
        if (mp) b=0
        do l=1,size(sims,1)
            if (sims(l,j)>=k) then
                p=p+1
                if (mp) b=b+1
            end if
!            write (33,*) p,b ! dd
        end do
        if (mp) write (id,*) b/size(sims,1)
    end if
!    write (33,*) j,a ! dd
end do
p=p/(size(sims,1)*a)
if (mp) close(id)
call cpu_time(finish)
timer=finish-start
end subroutine likely

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

end program
