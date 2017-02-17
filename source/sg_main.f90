program subgrider

! declaracao de variaveis
integer::i,op,nvar,fid
integer,dimension(3)::dg,pa,pb
character(len=256)::ficheiro,output !tentar len=* com f03 ou flibs
character(len=1)::sgems,load
logical::file_exists
! tipo de dados grid
type::grid
    integer::dx,dy,dz
    real,allocatable::val(:)
end type grid
type(grid)::res
fid=10
! input inicial
print *,"----|| s u b g r i d e r  v0.11 ||----"
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
print *,"a carregar o ficheiro..."
call abre(ficheiro,nvar,dg,res,fid)
print *,"ficheiro carregado"

!ciclo principal
do
print *,"escolher uma opcao"
print *,"1 - obter um poco vertical a partir das coordenadas (x,y)"
print *,"2 - obter outra grid a partir de um cubo ou paralelipipedo"
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
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    call novo(sgems,nvar+3,fid,output)
    call pp(dg,pa(1),pa(2),res,fid,output,sgems,nvar)
    print *,"operacao concluida"
elseif (op==2) then
    print *,"opcao 2"
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
    elseif (sgems=="n" .or. sgems=="N") then
        sgems="n"
    else
        print *,"erro - nao sei ler isso"
        stop
    end if
    call novo(sgems,nvar,fid,output)
    call subgrid(pa,pb,res,fid,output,sgems,nvar)
    print *,"operacao concluida"
elseif (op==0) then
    print *,"programa terminado"
    stop
end if
end do

! funcoes e subrotinas

contains

! carrega o ficheiro de dados para um vector
subroutine abre(ficheiro,nvar,dg,grida,id)
character(len=256), intent(in) :: ficheiro !tentar len=* com f03 ou flibs
integer, intent(in) :: dg(3),id
integer, intent(inout) :: nvar
type(grid),intent(out)::grida
real::start,finish ! timer
call cpu_time(start) ! timer
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
call cpu_time(finish) ! timer
print '("tempo = ",f6.3," segundos.")',finish-start ! timer
end subroutine abre

! cria um novo ficheiro de texto com cabecario sgems (ou nao)
subroutine novo(sgems,nvar,id,newfile)
character(len=1),intent(in)::sgems
integer,intent(in)::nvar
integer,intent(inout)::id
character(len=256),intent(out)::newfile !tentar len=* com f03 ou flibs
character(len=50)::nv,novonome  !tentar len=* com f03 ou flibs
id=id+1
print *,"nome do ficheiro"
read *,newfile
open (id,file=newfile,action='write')
if (sgems=="s") then
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
close(id)
end subroutine novo

! papa um poco a partir de uma grid e devolve um point set
subroutine pp(dg,xp,yp,res,id,output,sgems,nvar)
integer,intent(in)::dg(3),xp,yp,id,nvar
type(grid),intent(in)::res
character(len=256),intent(in)::output !tentar len=* com f03 ou flibs
character(len=1),intent(in)::sgems
real::start,finish ! timer
integer::z,p
print *,"a furar o poco..."
call cpu_time(start) ! timer
open (id,file=output)
if (sgems=="s") then
    do i=1,nvar+5
        read (id,*)
    end do
end if
do z=1,dg(3)
	p=xp+dg(1)*(yp-1)+dg(1)*dg(2)*(z-1)
	write(id,*) xp,yp,z,res%val(p)
end do
close(id)
call cpu_time(finish) ! timer
print '("tempo = ",f6.3," segundos.")',finish-start ! timer
end subroutine pp

! cria uma grid a partir de outra existente
subroutine subgrid(pa,pb,res,id,output,sgems,nvar)
integer,dimension(3),intent(in)::pa,pb
integer,intent(in)::id,nvar
type(grid),intent(in)::res
character(len=256),intent(in)::output !tentar len=* com f03 ou flibs
character(len=1),intent(in)::sgems
integer::pi,pf,d,p(3),j,c,v
real::start,finish
print *,"a criar a subgrid..."
call cpu_time(start) ! timer
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
call cpu_time(finish) ! timer
print '("tempo = ",f6.3," segundos.")',finish-start ! timer
close(id)
end subroutine subgrid

end program
