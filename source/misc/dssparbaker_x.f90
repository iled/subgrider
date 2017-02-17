program dssparbaker_x

use sg_utils
use qsort_c
use sg_gridutils
use sg_griders

implicit none
integer:: i_k,seed,i_skip,k
character(len=256):: par_o,par_sim,output,par_n
character(len=4)::n_k,format
character(len=70)::copy
real:: start,finish,timer


!numero de realizacoes
k=30

call cpu_time(start)
call chdir("/Users/julio/Desktop/hm-teste")
print *,"a criar PARs..."
!ficheiro par original
par_n="dss"
par_o=trim(par_n)//'.par'
do i_k=1,k
    !indice da realizacao
    if (i_k<10) format="(I1)"
    if (i_k>=10 .and. i_k<100) format="(I2)"
    if (i_k>=100 .and. i_k<1000) format="(I3)"
    if (i_k>=1000 .and. i_k<10000) format="(I4)"
    write (n_k,format) i_k
    !output da sim e ficheiro par
    output=trim(par_n)//'_'//trim(n_k)//'.out'
    par_sim='dss_'//trim(n_k)//'.par'
    !nova semente
    call seeder(seed)
    !ciclo para editar ficheiro
    call chdir("/Users/julio/Desktop/hm-teste")
    open(21,file=par_o,action='read') !21=par original
    call chdir("/Users/julio/Desktop/hm-teste/pars")
    open(20,file=par_sim,action='write') !20=novo
    do i_skip=1,37
        if (i_skip==16) then !caminho e ficheiro de saida
            read (21,*)
            write(20,*) trim(output)
        elseif (i_skip==23) then !semente
            read (21,*)
            write(20,*) seed
        else !copiar o resto
            read (21,'(A)') copy
            write(20,'(A)') copy
        end if
    end do
    close(20)
    close(21)
end do
call cpu_time(finish)
timer=finish-start
print *,k," PARs criados em ",tempo(timer)

contains

subroutine seeder(seed)
integer,intent(out)::seed
integer::a,b
a=100000
b=1000000
call rand_int(a,b,seed)
if ((seed .and. 1) .ne. 1) seed=seed+1
end subroutine seeder

end program dssparbaker_x
