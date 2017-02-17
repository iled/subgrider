program gs_to_ecl

use sg_gridutils
use qsort_c

implicit none
integer:: i_k,k,st
character(len=256)::outputx,outputy,outputz,gs_i,gs,dir,temp
character(len=4)::n_k,format
real:: start,finish,timer,copy

!numero de simulacoes
k=30

call cpu_time(start)
!call chdir("/Users/julio/Desktop/hm-teste/sims_0")
print *,"a converter..."
!base do ficheiro gslib
gs_i="dss_"
do i_k=1,k
    !indice da realizacao
    if (i_k<10) format="(I1)"
    if (i_k>=10 .and. i_k<100) format="(I2)"
    if (i_k>=100 .and. i_k<1000) format="(I3)"
    if (i_k>=1000 .and. i_k<10000) format="(I4)"
    write (n_k,format) i_k
    !ficheiro gslib
    gs=trim(gs_i)//trim(n_k)//'.out'
    dir=trim(gs_i)//trim(n_k)
    outputx='PermX.txt' !trim(gs_i)//trim(n_k)//'_x'//'.txt'
    outputy='PermY.txt' !trim(gs_i)//trim(n_k)//'_y'//'.txt'
    outputz='PermZ.txt' !trim(gs_i)//trim(n_k)//'_z'//'.txt'
    call chdir("/Users/julio/Desktop/hm-teste/sims_0")
    open(19,file=gs,action='read') !19=original
    call chdir("/Users/julio/Desktop/hm-teste/sims_0/ecl")
    call system('mkdir -p '//trim(dir))
    call chdir('/Users/julio/Desktop/hm-teste/sims_0/ecl/'//trim(dir))
    open(20,file=outputx,action='write') !20=novo x
    open(21,file=outputy,action='write') !21=novo y
    open(22,file=outputz,action='write') !22=novo z
    write (20,'(A)') "PERMX"
    write (21,'(A)') "PERMY"
    write (22,'(A)') "PERMZ"
    do
        read(19,*,iostat=st) copy
        if (st<0) then
            exit
        end if
        write (temp,*) copy
        write(20,'(A)') temp(front_trim(temp):len_trim(temp))
        write(21,'(A)') temp(front_trim(temp):len_trim(temp))
        write (temp,*) copy/10
        write(22,'(A)') temp(front_trim(temp):len_trim(temp))
    end do
    write (20,'(A)')
    write (20,'(A)') "/"
    write (21,'(A)')
    write (21,'(A)') "/"
    write (22,'(A)')
    write (22,'(A)') "/"
    close (20)
    close (21)
    close (22)
    close (19)
end do

call cpu_time(finish)
timer=finish-start

print *,"feito em",timer

end program gs_to_ecl
