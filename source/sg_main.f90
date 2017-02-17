program subgrider

use qsort_c_module
implicit none

! declaracao de variaveis
integer::i,op,nvar,fid,batch,nr_sims,nr_p,zef
integer,dimension(3)::dg,pa,pb,dg_2d
character(len=256)::ficheiro,output,pocos,base,mp_bend
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
type(grid)::res,hor
fid=10
! input inicial
print *,"----|| s u b g r i d e r  v0.19.3 ||----"
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
print *,"1 - obter um poco vertical a partir das coordenadas (x,y)"
!print *,"1 ----- poco vertical "
!print *,"  1.1 - a partir das coordenadas (x,y)"
print *,"2 - obter n pocos verticais a partir de um ficheiro com coordenadas (x,y)"
print *,"3 - obter uma grid a partir um volume (cubo, paralelipipedo) da grid inicial"
print *,"4 - criar uma mascara (data=1, no data=0)"
print *,"5 - calcular verosimilhanca"
print *,"6 - calcular percentil"
print *,"7 - corrigir quantil"
print *,"8 - bender"
print *,"0 - sair"
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
    print *,trim(pocos)," carregado em ",tempo(timer)
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
! <alt>
    print *,"base"
    read *,base
! <\alt>
!    base=ficheiro(1:len_trim(ficheiro)-4)//'_sim'
    call header_ask("r",header,nvar)
    print *,"a carregar grids simuladas..."
    call abre_sims(nr_sims,sims,dg,base,'.OUT',header,timer)
    print *,nr_sims,"grids carregadas em ",tempo(timer)
	do
	    print *,"produzir grid com mapa de probabilidades? (s/n)"
    	read *,mp_p
	    if (mp_p=="s" .or. mp_p=="S") then
    	    mp=.true.
    	    exit
    	elseif (mp_p=="n" .or. mp_p=="N") then
        	mp=.false.
        	exit
    	else
        	print *,"erro - input invalido."
		end if
	end do
    if (mp) then
    	call header_ask("w",header,nvar)
    	call novo(header,nvar,fid,output,batch,timer)
	end if
	print *,"inserir valor x (P[X=x])"
	read *,P_x
    print *,"a calcular probabilidades..."
    call likely(res,sims,P_x,P1,mp,fid,output,header,timer)
    print *,"a probabilidade e' ",trim(rtochar(P1))
    print *,"operacao concluida em ",tempo(timer)
    ! colocado aqui temporariamente, incorporar como batch do percentil
    print *,""
    print *,"percentis das simulacoes"
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"inserir percentil"
    read *,p
    print *,"a calcular percentis..."
    call percentil_sims(sims,p,P_x,fid,header,output,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==6) then
	print *,"opcao 6"
	print *,"inserir percentil a calcular"
	read *,p
	call percentil(res,p,q,timer)
	print '("o percentil ",f4.2," e ",a)',p,trim(rtochar(q))
	print *,"operacao concluida em ",tempo(timer)
elseif (op==7) then
	print *,"opcao 7"
	print *,"inserir valor"
	read *,q
	call percentil_updt(res,q,p,timer)
	print *,"o valor ",trim(rtochar(q))," corresponde ao percentil ",trim(rtochar(p))
	print *,"operacao concluida em ",tempo(timer)
elseif (op==8) then
    print *,"opcao 8"
    print *,"atencao: a grid inicial (alocada) vai ser alterada"
    print *,"dimensao maxima em z a escrever"
    read *,zef
    print *,"mapa de planificacao a abrir: caminho/nome"
    read *,mp_bend
    call checkfile(mp_bend)
	call header_ask("r",header,nvar)
    print *,"dimensoes do mapa: x y"
    read (*,*) dg_2d(1),dg_2d(2)
    dg_2d(3)=1
    print *,"a ler mapa de planificacao..."
    call abre(mp_bend,header,dg_2d,nd,hor,fid+20,timer) !nd
    print *,"mapa ",trim(mp_bend)," lido em ",tempo(timer)
	call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a planificar grid..."
	call bender(dg(1),dg(2),dg(3),-999.25,fid,output,header,zef,hor,res,timer)
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
var=sum(data**2)-sum(data)**2
end function variancia

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
	batch=batch+res%dz
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
integer::j,l,a
real::start,finish,b
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
    else
        if (mp) write (id,*) res%nd
    end if
!    write (33,*) j,a ! dd
end do
p=p/(size(sims,1)*a)
if (mp) close(id)
call cpu_time(finish)
timer=finish-start
end subroutine likely

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
real,intent(in)::p
real,intent(in)::q
integer,intent(in)::id
logical,intent(in)::header
character(len=256)::output
real,intent(out)::timer
type(grid)::aux
integer::i
real,allocatable,dimension(:)::ps,qs
real::start,finish,t,psim,qsim,pp
call cpu_time(start)
open (id,file=output)
if (header) call  header_skip(id)
allocate(ps(size(sims,1)))
allocate(qs(size(sims,1)))
allocate(aux%val(size(sims,2)))
pp=p
do i=1,size(sims,1)
    aux%val=sims(i,:)
    call percentil(aux,pp,qsim,t)
    call percentil_updt(aux,q,psim,t)
    ps(i)=qsim
    qs(i)=psim
    write (id,*) qsim,psim
end do
write (id,*) "media ",media(ps),media(qs)
write (id,*) "variancia ",variancia(ps),variancia(qs)
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

! planifica uma grid a partir de um mapa de planificacao previamente carregado
subroutine bender(x,y,z,nd,id,output,header,zef,hor,res,timer)
integer::yc,xc,zc,k
integer,intent(in)::x,y,z,zef,id
real,intent(in)::nd
type(grid),intent(inout)::res
type(grid),intent(in)::hor
real,intent(out)::timer
character(len=256),intent(in)::output
logical,intent(in)::header
real::start,finish
integer::h
call cpu_time(start)
yc=0
do while(yc<y)
    xc=0
    do while(xc<x)
        zc=0
        do while(zc<z)
        	h=hor%val(x*yc+xc+1) !tentar com nint(real)
            if (zc<z-h) then
                res%val(x*y*zc+x*yc+xc+1)=res%val(1+x*y*(h+zc)+x*yc+xc)
            else
                res%val(1+x*y*zc+x*yc+xc)=nd
            end if
            zc=zc+1
        end do
        xc=xc+1
    end do
    yc=yc+1
end do
k=0
open (id,file=output,action='write')
if (header) call header_skip(id)
do while(k<x*y*zef)
    write (id,*) res%val(k+1)
    k=k+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine bender

end program
