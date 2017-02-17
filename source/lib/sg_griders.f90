! subgrider :: Julio Caineta, 2010
! sg_griders :: modulo de manipulacao de grids

module sg_griders

use sg_gridutils
use pointers
implicit none

public::abre,bender,likely,mask,pp,subgrid,simmedvar,unbender,blocking_dim,blocking_num,bgeost_ps,refin

private

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
integer::j,c,v,i
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (id,file=output)
if (header) call header_skip(id)
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

! discretiza uma grid em blocos e gera output para BGeost
! recebe as dimensoes pretendidas do bloco
subroutine blocking_dim(pt_size,bl_dim,bl_e,res,id,output,bl_n,timer)
integer,dimension(3),intent(in)::bl_dim
real,intent(in)::bl_e,pt_size(3)
integer,intent(in)::id
type(grid),intent(in)::res
character(len=256),intent(in)::output
integer,intent(out)::bl_n
integer::lyr,i_x,i_y,i_z,bl_pi(3),bl_pf(3),i_list,status,i_pset
character(len=4)::format,n_b
character(len=64)::copy
real,intent(out)::timer
real::start,finish,bl_m
type(lista)::temp_list
real,allocatable,dimension(:)::temp_arr
real,allocatable,dimension(:,:)::temp_pset
call cpu_time(start)
open (id,status='scratch') !file='tmp.log')
write (id,*) trim(output)//"_blocks"
write (id,*)
bl_n=0
bl_pi=(/1,1,1/)
!i_blz=dg(3)-1
bl_pf(3)=bl_pi(3)+bl_dim(3)-1
do while (bl_pf(3)<=res%dz)
    bl_pi(2)=1
    bl_pf(2)=bl_pi(2)+bl_dim(2)-1
    do while (bl_pf(2)<=res%dy)
        bl_pi(1)=1
        bl_pf(1)=bl_pi(1)+bl_dim(1)-1
        do while (bl_pf(1)<=res%dx)
            i_pset=1
            bl_n=bl_n+1
            !allocate(temp_arr(bl_pf(1)-b))
            allocate(temp_pset(res%dx*res%dy*res%dz,3)) !aldrabado
		    if (bl_n<10) format="(I1)"
		    if (bl_n>=10 .and. bl_n<100) format="(I2)"
		    if (bl_n>=100 .and. bl_n<1000) format="(I3)"
		    if (bl_n>=1000 .and. bl_n<10000) format="(I4)"
		    write (n_b,format) bl_n
		    write (id,*) "block_"//trim(n_b)
		    temp_list=lista_nova()
		    ! subgrider
		    lyr=res%dx*res%dy
		    i_z=bl_pi(3)-1
		    do while (i_z<bl_pf(3))
		        i_y=i_z*lyr+res%dx*(bl_pi(2)-1)+bl_pi(1)
		        do while (i_y<=(bl_pf(2)-1)*res%dx+bl_pi(1)+i_z*lyr)
		            i_x=0
		            do while (i_x<=(bl_pf(1)-bl_pi(1)))
		                call lista_adf(temp_list,res%val(i_y+i_x))
		                temp_pset(i_pset,:)=get_3dcoord(0.0,0.0,0.0,pt_size,res%dx,res%dy,res%dz,i_y+i_x) !xi,yi,zi,bl_size
		                i_x=i_x+1
		                i_pset=i_pset+1
		            end do
		            i_y=i_y+res%dx
		        end do
		        i_z=i_z+1
		    end do
		    ! /subgrider
		    call lista_array(temp_list,temp_arr)
		    write (id,*) sum(temp_arr,mask=temp_arr.ne.res%nd)/real(count(temp_arr.ne.res%nd))
		    write (id,*) bl_e
		    do i_list=1,lista_comp(temp_list)
                write (id,*) temp_pset(i_list,:),temp_arr(i_list)
		    end do
            bl_pf(1)=bl_pf(1)+bl_dim(1)
            if (bl_pf(1)-res%dx>0 .and. bl_pf(1)-res%dx<bl_dim(1)) bl_pf(1)=res%dx
            bl_pi(1)=bl_pi(1)+bl_dim(1)
            deallocate(temp_pset)
        end do
        bl_pf(2)=bl_pf(2)+bl_dim(2)
        if (bl_pf(2)-res%dy>0 .and. bl_pf(2)-res%dy<bl_dim(2)) bl_pf(2)=res%dy
        bl_pi(2)=bl_pi(2)+bl_dim(2)
    end do
    bl_pf(3)=bl_pf(3)+bl_dim(3)
    if (bl_pf(3)-res%dz>0 .and. bl_pf(3)-res%dz<bl_dim(3)) bl_pf(3)=res%dz
    bl_pi(3)=bl_pi(3)+bl_dim(3)
end do
open (id+2,file=output)
rewind(id)
read (id,'(A)') copy
write (id+2,'(A)') copy(front_trim(copy):len_trim(copy))
read (id,*)
write (copy,*) bl_n
write (id+2,*) copy(front_trim(copy):len_trim(copy))
status=0
do
    read (id,'(A)',iostat=status) copy
    if (status<0) exit
    write (id+2,'(A)') copy(front_trim(copy):len_trim(copy))
end do
close(id)
close(id+2)
call cpu_time(finish)
timer=finish-start
end subroutine blocking_dim

! discretiza uma grid em blocos e gera output para BGeost
! recebe o numero de blocos pretendidos
subroutine blocking_num(bl_num,bl_e,res,id,output,bl_d,timer)
real,intent(in)::bl_e
integer,intent(in)::bl_num,id
type(grid),intent(in)::res
character(len=256),intent(in)::output
integer,intent(out),dimension(3)::bl_d
integer::lyr,i_x,i_y,i_z,bl_pi(3),bl_pf(3)
character(len=4)::format,n_b
real,intent(out)::timer
real::start,finish,bl_m
call cpu_time(start)
open (id,file=output)
write (id,*) trim(output)//"_blocks"
write (id,*) bl_num

bl_d=0 !temp

close(id)
call cpu_time(finish)
timer=finish-start
end subroutine blocking_num

! transforma um ficheiro tipo BGeost em point set
subroutine bgeost_ps(id,input,output,timer)
integer,intent(in)::id
character(len=*),intent(in)::input,output
integer::j,nb,st,pt,bl,pts
character(len=128)::temp
character(len=1)::b
character(len=32)::x,y,z,val,n
real,intent(out)::timer
real::start,finish
call cpu_time(start)
open (id,file=input,action='read')
open (id+3,file=output,action='write')
read (id,*)
read (id,*) nb
write (id+3,*) "block_to_pointset"
write (id+3,*) 6
write (id+3,*) "x"
write (id+3,*) "y"
write (id+3,*) "z"
write (id+3,*) "block"
write (id+3,*) "average"
write (id+3,*) "noise"
pt=1
bl=0
do
    b=' '
    do while (b==' ')
        read(id,'(A)',iostat=st,advance="no") b
        if (st<0) then
            if (nb.ne.bl) then
                print *,"aviso: numero de blocos inconsistente",bl," blocos encontrados"
            end if
            close(id)
			close(id+3)
			call cpu_time(finish)
			timer=finish-start
            return
        end if
    end do
    if (iachar(b)>57 .or. iachar(b)<44) then
        bl=bl+1
        pt=0
        read (id,*)
        read (id,*) val
        read (id,*) n
    else
        backspace(id)
        read (id,*) x,y,z
        write (id+3,*) trim(x)," ",trim(y)," ",trim(z)," ",bl," ",trim(val)," ",trim(n)
        pt=pt+1
        pts=pts+pt
    end if
end do
end subroutine bgeost_ps

!refinamento sem interpolacao de uma grid (blocos em pontos)
subroutine refin(bl_dim,header,res,id,output,timer)
integer,dimension(3),intent(in)::bl_dim
integer,intent(in)::id
type(grid),intent(in)::res
logical,intent(in)::header
character(len=256),intent(in)::output
integer::i,j,k,coord_o(3),coord_n(3),nc
character(len=4)::format,n_b
real,intent(out)::timer
real::start,finish
type(lista)::temp_list
real,allocatable,dimension(:)::temp_arr
real,allocatable,dimension(:,:)::temp_pset
call cpu_time(start)
open (id,file=output)
if (header) call header_skip(id)
k=0
do while (k<res%dz*bl_dim(3))
    j=0
    do while (j<res%dy*bl_dim(2))
        i=0
        do while (i<res%dx*bl_dim(1))
            coord_n=(/i+1,j+1,k+1/)
            do nc=1,3
                coord_o(nc)=int(coord_n(nc)/real(bl_dim(nc))+0.5)
            end do
            write(id,*) res%val(coord_o(1)+res%dx*(coord_o(2)-1)+res%dx*res%dy*(coord_o(3)-1))
            i=i+1
        end do
        j=j+1
    end do
    k=k+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine refin

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
integer::j,l,a,pp,b
real::start,finish
call cpu_time(start)
if (mp) then
    open(id,file=output)
    if (header) call header_skip(id)
end if
a=0
p=0.0
pp=0
! open(33,file='likely.dbg',action='write') ! dd
! write (33,*) size(res%val), size(sims,1) ! dd
do j=1,size(res%val)
    if (res%val(j)>=k(1).and.res%val(j)<k(2)) then
        a=a+1
        if (mp) b=0
        do l=1,size(sims,1)
            if (sims(l,j)>=k(1) .and. sims(l,j)<k(2)) then
                pp=pp+1
!               write (33,*) j, res%val(j), sims(l,j)
                if (mp) b=b+1
            end if
        end do
        if (mp) write (id,*) b/real(size(sims,1))
    else
        if (mp) write (id,*) res%nd
    end if
end do
p=pp/real((size(sims,1)*a))
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

subroutine unbender(x,y,z,nd,id,output,header,hor,res,timer)
integer::yc,xc,zc,k
integer,intent(in)::x,y,z,id
real,intent(in)::nd
type(grid),intent(in)::res
type(grid),intent(in)::hor
real,intent(out)::timer
character(len=256),intent(in)::output
logical,intent(in)::header
real::start,finish
real,allocatable,dimension(:)::un
integer::h,maxi
call cpu_time(start)
allocate(un(x*y*(z+int(maxval(hor%val)))))
yc=0
un(:)=nd
do while(yc<y)
    xc=0
    do while(xc<x)
        zc=0
        do while(zc<z)
            h=hor%val(x*yc+xc+1) !tentar com nint(real)
            un(x*y*(zc+h)+x*yc+xc+1)=res%val(1+x*y*(h+zc)+x*yc+xc)
            zc=zc+1
        end do
        xc=xc+1
    end do
    yc=yc+1
end do
k=0
print *,""
open (id,file=output,action='write')
if (header) call header_skip(id)
maxi=z+maxval(hor%val)
do while(k<x*y*maxi)
    write (id,*) un(k+1)
    k=k+1
end do
close(id)
call cpu_time(finish)
timer=finish-start
end subroutine unbender

subroutine simmedvar(do_med,do_var,sims,id,output,header,timer)
real,dimension(:,:),intent(in)::sims
integer,intent(in)::id
character(len=256),intent(in)::output
logical,intent(in)::header,do_med,do_var
real,intent(out)::timer
integer::i,j
real::start,finish
!real,dimension(:),allocatable::med,var
call cpu_time(start)
!if (do_med) allocate(med(size(sims,2)))
!if (do_var) allocate(var(size(sims,2)))
open (id,file=output)
if (header) call header_skip(id)
do j=1,size(sims,2)
!    if (do_med) med(j)=media(sims(:,j))
!    if (do_var) var(j)=variancia(sims(:,j))
    if (do_med .and. .not.do_var) write(id,*) media(sims(:,j))
    if (.not.do_med .and. do_var) write(id,*) variancia(sims(:,j))
    if (do_med .and. do_var) write(id,*) media(sims(:,j)),variancia(sims(:,j))
end do
call cpu_time(finish)
timer=finish-start
end subroutine simmedvar

end module
