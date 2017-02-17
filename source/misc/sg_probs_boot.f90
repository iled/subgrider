program sg_probs_boot

use sg_gridutils
use sg_utils

implicit none
integer::i_cen,nr_sims,i_sim,len_well,i_cl,wi,wf,n
character(len=256)::path_well,well
character(len=4)::n_sim,format
real::start,finish,timer,k(2),p
real,allocatable::w(:),wells(:)
!type(grid)::aux

nr_sims=50
!aux%nd=-999.0
call chdir("/Users/julio/Desktop/sims")
open(9,file='probs_new_Z.txt',action='write')

do i_cen=1,9
    if (i_cen==1) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w3_G"
        well='wells_w3_G_'
        n=3
    elseif (i_cen==2) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w3_M"
        well='wells_w3_M_'
        n=3
    elseif (i_cen==3) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w3_P"
        well='wells_w3_P_'
        n=3
    elseif (i_cen==4) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w5_G"
        well='wells_w5_G_'
        n=5
    elseif (i_cen==5) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w5_M"
        well='wells_w5_M_'
        n=5
    elseif (i_cen==6) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w5_P"
        well='wells_w5_P_'
        n=5
    elseif (i_cen==7) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w10_G"
        well='wells_w10_G_'
        n=10
    elseif (i_cen==8) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w10_M"
        well='wells_w10_M_'
        n=10
    elseif (i_cen==9) then
        path_well="/Users/julio/Desktop/sims/bootstrap/wells_w10_P"
        well='wells_w10_P_'
        n=10
    else
        close(9)
        stop
    end if

write(9,*) path_well(27:len_trim(path_well))

call cpu_time(start)
print *,trim(path_well)
print *,"a alocar pocos"
call chdir(path_well)
len_well=len_trim(well)
do i_sim=1,nr_sims
    if (i_sim<10) format="(I1)"
    if (i_sim>=10 .and. i_sim<100) format="(I2)"
    if (i_sim>=100 .and. i_sim<1000) format="(I3)"
    if (i_sim>=1000 .and. i_sim<10000) format="(I4)"
    write (n_sim,format) i_sim
    well=well(1:len_well)//trim(n_sim)//'.prn'
    call abre_pocos_val(well,n*93,w,timer)
    !print *,i_sim,nr_sims,size(w),nr_sims*size(w)
    if (allocated(wells)) then
        if (size(wells).ne.nr_sims*size(w)) then
            deallocate(wells)
            allocate(wells(nr_sims*size(w)))
            wi=1
            wf=size(w)
        end if
    else
        allocate(wells(nr_sims*size(w)))
        wi=1
        wf=size(w)
    end if
    wells(wi:wf)=w
    wi=wf+1
    wf=wf+size(w)
end do
call cpu_time(finish)
timer=finish-start
print *,wells(23),wells(334),wells(9),size(wells)
print *,"pocos alocados em ",tempo(timer)
call cpu_time(start)
print *,"a calcular e a escrever resultados"
do i_cl=1,3
    if (i_cl==1) k=(/ 25.7612, 40.0/)
    if (i_cl==2) k=(/ 15.3437, 25.7612 /)
    if (i_cl==3) k=(/ 0.0, 15.3437 /)
    call quantil_range(wells,k(1),k(2),p,timer)
    write(9,*) "classe: ",k
    write(9,*) "probabilidade: ",p
end do
call cpu_time(finish)
timer=finish-start
print *,"concluido em ",tempo(timer)
print *,""
write(9,*) ""
deallocate(wells)
end do

contains

subroutine abre_pocos_val(pocos,n,w,timer)
character(len=256),intent(in)::pocos
real,allocatable,intent(out)::w(:)
integer,intent(in)::n
integer::i
real,intent(out)::timer
real::start,finish,aux(3)
call cpu_time(start)
open (19,file=pocos,action='read')
allocate(w(n))
do i=1,n
    read (19,*) aux(1),aux(2),aux(3),w(i)
end do
close(19)
call cpu_time(finish)
timer=finish-start
end subroutine abre_pocos_val

subroutine quantil_range(data,a,b,q,timer)
real,allocatable,dimension(:),intent(in)::data
real,intent(in)::a,b
real,intent(out)::q,timer
real::start,finish
integer::i,n
call cpu_time(start)
n=0
do i=1,size(data)
    if (data(i)>=a .and. data(i)<b) n=n+1
end do
q=n/real(size(data))
call cpu_time(finish)
timer=finish-start
end subroutine quantil_range

end program sg_probs_boot
