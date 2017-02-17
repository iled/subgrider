program dssparbaker

use sg_utils
use qsort_c
use sg_gridutils

implicit none
integer::i_cen,i_wells,k,i_k,i_var,seed,var(3,3),n,i_skip,cut_sim
character(len=256)::cen,sim,path_sim,input,path_in,input_set,path_well
character(len=1)::n_cen,n_var
character(len=4)::n_k,format
character(len=70)::copy
real::start,finish,timer,minw,maxw
real,allocatable,dimension(:)::set,w

print *,"usar K sets de pocos"
read *,k

var=reshape((/ (/ 165, 65, 25 /), (/ 110, 40, 25 /), (/ 60, 25, 25 /) /),(/3, 3/))

call cpu_time(start)
call chdir("/Users/julio/Desktop/sims/pars_cen")
print *,"a criar PARs..."
do i_cen=1,9
	cen="sim"
	write (n_cen,"(I1)") i_cen
	cen=trim(cen)//n_cen//'.out'
    do i_wells=1,9
        sim="dss_"
        if (i_wells==1) then
            sim=trim(sim)//'w3_G'
            n=3
            path_well='wells_w3_G'
        elseif (i_wells==2) then
            sim=trim(sim)//'w3_M'
            n=3
            path_well='wells_w3_M'
        elseif (i_wells==3) then
            sim=trim(sim)//'w3_P'
            n=3
            path_well='wells_w3_P'
        elseif (i_wells==4) then
            sim=trim(sim)//'w5_G'
            n=5
            path_well='wells_w5_G'
        elseif (i_wells==5) then
            sim=trim(sim)//'w5_M'
            n=5
            path_well='wells_w5_M'
        elseif (i_wells==6) then
            sim=trim(sim)//'w5_P'
            n=5
            path_well='wells_w5_P'
        elseif (i_wells==7) then
            sim=trim(sim)//'w10_G'
            n=10
            path_well='wells_w10_G'
        elseif (i_wells==8) then
            sim=trim(sim)//'w10_M'
            n=10
            path_well='wells_w10_M'
        elseif (i_wells==9) then
            sim=trim(sim)//'w10_P'
            n=10
            path_well='wells_w10_P'
        end if
        path_sim=sim
        !call system("mkdir "//trim(path_sim))
        !print *,"-- set: ",trim(path_sim)
        path_in="bootstrap/wells"//sim(4:len_trim(sim))//"/"
        sim=trim(sim)//'_'
        do i_k=1,k
            if (i_k<10) format="(I1)"
	        if (i_k>=10 .and. i_k<100) format="(I2)"
	        if (i_k>=100 .and. i_k<1000) format="(I3)"
	        if (i_k>=1000 .and. i_k<10000) format="(I4)"
	        write (n_k,format) i_k
            sim=sim(1:cut_sim-2)//trim(n_k)//'_' !trim(sim)//trim(n_k)//'_'
            cut_sim=len_trim(sim)
            input_set=trim(path_well)//'_'//trim(n_k)//'.prn'
            input=trim(path_in)//trim(input_set) !trim(path_in)//'wells'//sim(4:cut_sim-1)//'.prn'
            do i_var=1,3
                if (i_var==1) n_var='G'
                if (i_var==2) n_var='M'
                if (i_var==3) n_var='P'
                sim=sim(1:cut_sim)//n_var//'.par'
                call seeder(seed)
                call chdir("/Users/julio/Desktop/sims/bootstrap")
                call chdir(path_well)
                call abre_pocos_val(input_set,n*93,w,timer)
                call del_sorted_nd(w,set,-999.0,timer)
                call chdir("/Users/julio/Desktop/sims/pars_cen")
                minw=minval(set)
                maxw=maxval(set)
                open(20,file=sim,action='write')
                open(21,file=cen,action='read')
                do i_skip=1,37
                    if (i_skip==3) then
                        read (21,*)
                        write (20,*) " "
                    elseif (i_skip==5) then
                        read (21,*)
                        write (20,*) trim(input)
                    elseif (i_skip==8) then
                        read (21,*)
                        write(20,*) minw,maxw
                    elseif (i_skip==11) then
                        read (21,*)
                        write(20,*) minw,maxw
                    elseif (i_skip==12) then
                        read (21,*)
                        write(20,*) 1,minw
                    elseif (i_skip==13) then
                        read (21,*)
                        write(20,*) 1,maxw
                    elseif (i_skip==16) then
                        read (21,*)
                        write(20,*) sim(1:len_trim(sim)-4)//'.out'
                    elseif (i_skip==23) then
                        read (21,*)
                        write(20,*) seed
                    elseif (i_skip==37) then
                        read (21,*)
                        write(20,*) var(:,i_var)
                    else
                        read (21,'(A)') copy
                        write(20,'(A)') copy
                    end if
                end do
                close(20)
                close(21)
                call system("mv "//trim(sim)//" "//trim(path_sim))
            end do
        end do
    end do
end do
call cpu_time(finish)
timer=finish-start
print *,k*27," PARs criados em ",tempo(timer)

contains

subroutine seeder(seed)
integer,intent(out)::seed
integer::a,b
a=100000
b=1000000
call rand_int(a,b,seed)
if ((seed .and. 1) .ne. 1) seed=seed+1
end subroutine seeder

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

end program dssparbaker
