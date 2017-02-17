module sg_griders

use sg_gridutils
implicit none

public::abre,bender,likely,mask,pp,subgrid

contains

! carrega o ficheiro de dados para um vector
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
integer::i
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
real,intent(in)::k(2)
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
    if (res%val(j)>=k(1).and.res%val(j)<k(2)) then
        a=a+1
        if (mp) b=0
        do l=1,size(sims,1)
            if (sims(l,j)>=k(1).and.sims(l,j)<k(2)) then
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

end module
