program simhist

use sg_griders

implicit none

logical::header
integer::dg(3),i_cen,id,nr_sims,i_sim,i_var,cl,len_sim
real::nd,start,finish,timer,k(2),q
character(len=256)::path_sim,sim,hist,hist1,hist2,hist3
character(len=1)::n_var
character(len=4)::format,simN
type(grid)::simgrid

header=.FALSE.
dg=(/ 334, 134, 93 /)
nd=-999.0
nr_sims=50
id=20


do i_cen=1,9
    if (i_cen==1) then
        path_sim="/home/data/sims/dss_w3_G"
        sim='dss_w3_G'
    elseif (i_cen==2) then
        path_sim="/home/data/sims/dss_w3_M"
        sim='dss_w3_M'
    elseif (i_cen==3) then
        path_sim="/home/data/sims/dss_w3_P"
        sim='dss_w3_P'
    elseif (i_cen==4) then
        path_sim="/home/data/sims/dss_w5_G"
        sim='dss_w5_G'
    elseif (i_cen==5) then
        path_sim="/home/data/sims/dss_w5_M"
        sim='dss_w5_M'
    elseif (i_cen==6) then
        path_sim="/home/data/sims/dss_w5_P"
        sim='dss_w5_P'
    elseif (i_cen==7) then
        path_sim="/home/data/sims/dss_w10_G"
        sim='dss_w10_G'
    elseif (i_cen==8) then
        path_sim="/home/data/sims/dss_w10_M"
        sim='dss_w10_M'
    elseif (i_cen==9) then
        path_sim="/home/data/sims/dss_w10_P"
        sim='dss_w10_P'
    else
        close(9)
        stop
    end if
    print *,trim(path_sim)
    hist=sim
    len_sim=len_trim(sim)
    do i_var=1,3
        if (i_var==1) n_var='G'
        if (i_var==2) n_var='M'
        if (i_var==3) n_var='P'
        hist=hist(1:len_sim)//'_'//n_var
        hist1=trim(hist)//'_'//'3Q'//'.txt'
        hist2=trim(hist)//'_'//'1Q3Q'//'.txt'
        hist3=trim(hist)//'_'//'1Q'//'.txt'
        call chdir("/home/data/sims")
        open(11,file=hist1,action='write')
        open(12,file=hist2,action='write')
        open(13,file=hist3,action='write')
        call chdir(path_sim)
        do i_sim=1,nr_sims
            if (i_sim<10) format="(I1)"
            if (i_sim>=10 .and. i_sim<100) format="(I2)"
            if (i_sim>=100 .and. i_sim<1000) format="(I3)"
            write (simN,format) i_sim
            sim=sim(1:len_sim)//'_'//trim(simN)//'_'//n_var//'.out'
            call checkfile(sim)
            call abre(sim,header,dg,nd,simgrid,id,timer)
            id=id+1
            do cl=1,3
                if (cl==1) k=(/ 25.7612, 40.0/)
                if (cl==2) k=(/ 15.3437, 25.7612 /)
                if (cl==3) k=(/ 0.0, 15.3437 /)
                call quantil_range(simgrid%val,k(1),k(2),q,timer)
                write(10+cl,*) q
            end do
            deallocate(simgrid%val)
        end do
        close(11)
        close(12)
        close(13)
    end do
end do

contains

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

end program simhist
