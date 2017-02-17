program subgrider

! declaracao de variaveis
integer::i,op,nvar,fid,batch,nr_sims,zef
integer,dimension(3)::dg,pa,pb,dg_2d
character(len=256)::ficheiro,output,pocos,base,mp_bend
character(len=1)::sgems,load,mp_p
logical::file_exists,mp
real,allocatable::sims(:,:)
integer,allocatable::px(:),py(:)
real::Q3,P1,timer,nd
! tipo de dados grid
type::grid
    integer::dx,dy,dz
    real,allocatable::val(:)
end type grid
type(grid)::res,hor
fid=10
! input inicial
print *,"----|| s u b g r i d e r  v0.16_bender.2 ||----"
print *,""
print *,"carregar o ficheiro anterior? (s/n)"
read *,load
if (load=="s" .or. load=="S") then
    inquire(file='subgrider.cfg',exist=file_exists)
    if (file_exists) then
	    open (9,file='subgrider.cfg',action='read')
	    read (9,*) ficheiro,sgems,nvar
	    read (9,*) dg
    else
        print *,"ficheiro de configuracao nao encontrado."
        stop
    end if
elseif (load=="n" .or. load=="N") then
    print *,"ficheiro a abrir: caminho/nome"
    read *,ficheiro
    print *,"ficheiro SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        nvar=0
	else
	    print *,"erro - nao sei ler isso"
	    stop
	end if
	print *,"dimensoes da grid: x y z"
	read (*,*) dg
	open (9,file='subgrider.cfg',action='write')
	write (9,*) ficheiro,sgems,nvar
	write (9,*) dg
	close (9)
else
    print *,"erro - nao sei ler isso"
    stop
end if
print *,"a carregar o ficheiro ",trim(ficheiro),"..."
call abre(ficheiro,nvar,dg,res,fid,timer)
print *,"ficheiro carregado em ",tempo(timer)

!ciclo principal
do
batch=-1
print *,""
print *,"escolher uma opcao"
print *,"1 - obter um poco vertical a partir das coordenadas (x,y)"
print *,"2 - obter n pocos verticais a partir do ficheiro pocos.dat"
print *,"3 - obter uma grid a partir um volume (cubo, paralelipipedo) da grid inicial"
print *,"4 - criar uma mascara (data=1, no data=0)"
print *,"5 - calcular verosimilhanca"
print *,"6 - bender"
print *,"0 - sair"
read *,op
if (op==1) then
    print *,"opcao 1"
    print *,"coordenadas do poco (x,y)"
    read *,pa(1),pa(2)
    print *,"output com cabecario SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        sgems="s"
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
        nvar=0
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    call novo(sgems,nvar+3,fid,output,batch,timer)
    print *,"a furar o poco..."
    call pp(dg,pa(1),pa(2),res,fid,output,sgems,nvar,batch,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==2) then
    print *,"opcao 2"
    write (pocos,*) "pocos.dat"
    print *,"a ler ficheiro ",trim(pocos),"..."
    call abre_pocos(pocos,px,py,timer)
    print *,"output com cabecario SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        sgems="s"
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
        nvar=0
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    batch=0
	call novo(sgems,nvar+3,fid,output,batch,timer)
	k=1
	print *,"a furar os pocos..."
    do k=1,10
    	call pp(dg,px(k),py(k),res,fid,output,sgems,nvar,batch,timer)
    end do
    print *,"operacao concluida em ",tempo(timer)
    batch=-1
elseif (op==3) then
    print *,"opcao 3"
    print *,"coordenadas x (min,max)"
    read *,pa(1),pb(1)
    print *,"coordenadas y (min,max)"
    read *,pa(2),pb(2)
    print *,"coordenadas z (min,max)"
    read *,pa(3),pb(3)
    print *,"output com cabecario SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        sgems="s"
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
        nvar=0
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    call novo(sgems,nvar,fid,output,batch,timer)
    print *,"a criar a subgrid..."
    call subgrid(pa,pb,res,fid,output,sgems,nvar,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==4) then
	print *,"opcao 4"
	print *,"inserir valor de no data"
	read *,nd
	call novo(sgems,nvar,fid,output,batch,timer)
	print *,"a criar mascara..."
	call mask(dg,1,0,nd,sgems,nvar,res,fid,output,timer)
	print *,"operacao concluida em ",tempo(timer)
elseif (op==5) then
    print *,"opcao 5"
    print *,"numero de simulacoes realizadas (ficheiros *_simX.OUT)"
    read *,nr_sims
    base=ficheiro(1:len_trim(ficheiro)-4)//'_sim'
    print *,"a carregar grids simuladas..."
    call abre_sims(nr_sims,sims,dg,base,'.OUT',nvar,timer)
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
    if (mp) call novo(sgems,nvar,fid,output,batch,timer)
    print *,"a calcular probabilidades..."
    Q3=28.0
    call likely(res,sims,Q3,P1,mp,fid,output,sgems,timer)
    print *,P1
    print *,"operacao concluida em ",tempo(timer)
elseif (op==6) then
    print *,"opcao 6"
    print *,"atencao: a grid inicial vai ser alterada"
    print *,"dimensao maxima em z a escrever"
    read *,zef
    print *,"mapa de planificacao a abrir: caminho/nome"
    read *,mp_bend
    print *,"ficheiro SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        nvar=0
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    print *,"dimensoes do mapa: x y"
    read (*,*) dg_2d(1),dg_2d(2)
    dg_2d(3)=1
    print *,"a ler mapa de planificacao"
    call abre(mp_bend,nvar,dg_2d,hor,id+20,timer)
    print *,"mapa lido em ",tempo(timer)
    print *,"output com cabecario SGeMS? (s/n)"
    read *,sgems
    if (sgems=="s" .or. sgems=="S") then
        sgems="s"
        nvar=1
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
        nvar=0
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    call novo(sgems,nvar,id,output,batch,timer)
    print *,"a planificar grid"
	call bender(dg(1),dg(2),dg(3),-999.25,id,output,zef,hor,res,timer)
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

! carrega o ficheiro de dados para um vector
subroutine abre(ficheiro,nvar,dg,grida,id,timer)
character(len=256), intent(in) :: ficheiro
integer, intent(in) :: dg(3),id
integer, intent(inout) :: nvar
type(grid),intent(out)::grida
real,intent(out)::timer
real::start,finish
call cpu_time(start)
grida%dx=dg(1)
grida%dy=dg(2)
grida%dz=dg(3)
allocate(grida%val(product(dg))) ! o vector-grid tem dimensao=produto(3dimensoes.grid)
open (id, file=ficheiro, action='read')
if (nvar>0) then
    read (id,*)
    read (id,*) nvar
    do i=1,nvar
        read (id,*)
    end do
end if
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
integer::n,i
logical::file_exists
real,intent(out)::timer
real::start,finish
call cpu_time(start)
inquire (file=pocos,exist=file_exists)
if (file_exists) then
    open (19,file=pocos,action='read')
    read (19,*) n
    allocate(px(n),py(n))
    do i=1,n
        read (19,*) px(i),py(i)
    end do
    close(19)
    print '("carregadas as coordenadas de ",I3," pocos.")',n
else
    print *,"ficheiro ",trim(pocos)," nao encontrado."
    stop
end if
call cpu_time(finish)
timer=finish-start
end subroutine abre_pocos

! carrega ficheiros de simulacoes para uma matriz
subroutine abre_sims(nr_sims,sims,dg,base,ext,nvar,timer)
integer,intent(in)::nr_sims,dg(3)
integer,intent(inout)::nvar
real,allocatable,intent(out)::sims(:,:)
character(len=256),intent(in)::base
character(len=4),intent(in)::ext
character(len=256)::sim
character(len=4)::simN,format
type(grid)::gridsim
logical::file_exists
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
    inquire (file=sim,exist=file_exists)
    if (file_exists) then
        call abre(sim,nvar,dg,gridsim,id,t)
        sims(i,:)=gridsim%val
        id=id+1
    else
        print *,"ficheiro ",trim(sim)," nao encontrado."
        stop
    end if
end do
call cpu_time(finish)
timer=finish-start
end subroutine abre_sims

! cria um novo ficheiro de texto com cabecario sgems (ou nao)
subroutine novo(sgems,nvar,id,newfile,batch,timer)
character(len=1),intent(in)::sgems
integer,intent(in)::nvar,batch
integer,intent(inout)::id
character(len=256),intent(out)::newfile
character(len=50)::nv,novonome
real,intent(out),optional::timer
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
if (sgems=="s") then
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
timer=finish-start
end subroutine novo

! papa um poco a partir de uma grid e devolve um point set
subroutine pp(dg,xp,yp,res,id,output,sgems,nvar,batch,timer)
integer,intent(in)::dg(3),xp,yp,id,nvar
integer,intent(inout)::batch
type(grid),intent(in)::res
character(len=256),intent(in)::output
character(len=1),intent(in)::sgems
real::start,finish
real,intent(out)::timer
integer::z,p,m
call cpu_time(start)
open (id,file=output)
if (sgems=="s") then
    do i=1,nvar+5
        read (id,*)
    end do
end if
if (batch>=0) then
	do m=1,batch
		read (id,*)
	end do
	batch=batch+93
end if
do z=1,dg(3)
	p=xp+dg(1)*(yp-1)+dg(1)*dg(2)*(z-1)
	write(id,*) xp,yp,z,res%val(p)
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine pp

! cria uma grid a partir de outra existente
subroutine subgrid(pa,pb,res,id,output,sgems,nvar,timer)
integer,dimension(3),intent(in)::pa,pb
integer,intent(in)::id,nvar
type(grid),intent(in)::res
character(len=256),intent(in)::output
character(len=1),intent(in)::sgems
integer::pi,pf,d,p(3),j,c,v
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (id,file=output)
if (sgems=="s") then
    do i=1,nvar+2
        read (id,*)
    end do
end if
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
subroutine mask(dg,m1,m2,nd,sgems,nvar,res,id,output,timer)
integer,dimension(3),intent(in)::dg
integer,intent(in)::id,nvar,m1,m2
type(grid),intent(in)::res
real,intent(in)::nd
character(len=256),intent(in)::output
character(len=1),intent(in)::sgems
integer::j,i
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (id,file=output)
if (sgems=="s") then
    do i=1,nvar+2
        read (id,*)
    end do
end if
j=1
do while (j<=product(dg))
	if (res%val(j)/=nd) then
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
subroutine likely(res,sims,k,p,mp,id,output,sgems,timer)
type(grid),intent(in)::res
real,dimension(:,:),intent(in)::sims
real,intent(in)::k
real,intent(out)::p,timer
integer,intent(in)::id
logical,intent(in)::mp
character(len=256),intent(in)::output
character(len=1),intent(in)::sgems
integer::j,l,a,b
real::start,finish
call cpu_time(start)
if (mp) open(id,file=output)
a=0
p=0
! open(33,file='likely.dbg',action='write') ! dd
! write (33,*) size(res%val), size(sims,1) ! dd
if (sgems=="s" .and. mp) then
    do i=1,nvar+2
        read (id,*)
    end do
end if
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

subroutine bender(x,y,z,nd,id,output,zef,hor,res,timer)
integer::yc,xc,zc,k
integer,intent(in)::x,y,z,zef
real,intent(in)::nd
type(grid),intent(inout)::res
type(grid),intent(in)::hor
real,optional,intent(out)::timer
character(len=256),intent(in)::output
real::start,finish
call cpu_time(start)
yc=0
do while(yc<y)
    xc=0
    do while(xc<x)
        zc=0
        do while(zc<z)
            if (hor%val(x*yc+xc+1)-zc>=0) then
                res%val(x*y*zc+x*yc+xc+1)=res%val(1+x*y*(hor%val(x*yc+xc+1)-zc)+x*yc+xc)
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
do while(k<x*y*zef)
    write (id,*) res%val(k+1)
    k=k+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine bender

end program
