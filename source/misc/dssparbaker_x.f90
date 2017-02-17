program dssparbaker_x2

! problema a escrever a linha do output, nao guarda caminho

use sg_utils
use qsort_c
use sg_gridutils
use sg_griders

implicit none
integer:: i_k,seed,i_skip,k,dss_v,lin_out,lin_seed,status
character(len=256):: par_o,par_sim,output,wd,par_b,out_b,par_p
character(len=4)::n_k,format
character(len=70)::copy
real:: start,finish,timer
logical::file_exists

print *,""
print *," - - -  dssparbaker-x v1.0 | julio 2011  - - - "
print *,""
print *,"preparar ficheiros de parametros para batch de simulacoes"
print *,""

print *,"escolher versao da dss:"
print *,"    1. geoMS ou NewIntel"
print *,"    2. paralela (RFMN)"
read *,dss_v
do
if (dss_v==1) then
    lin_out=19
    lin_seed=26
    exit
elseif (dss_v==2) then
    lin_out=16
    lin_seed=23
    exit
else
    print *,"opcao invalida"
end if
end do

print *,""

print *,"nome do ficheiro original"
read *,par_o
inquire(file=par_o,exist=file_exists)
if (.not. file_exists) then
    print *,"ficheiro ",trim(par_o)," nao encontrado."
    stop
end if

print *,"numero de simulacoes (ex: 30)"
read *,k

print *,"base dos ficheiros de parametros a criar (ex: dss_)"
read *,par_b

print *,"base dos ficheiros de output da dss (ex:dss_)"
read *,out_b

print *,"pasta onde ficheiros serao guardados (ex: pars)"
read *,par_p
inquire(file=par_p,exist=file_exists)
if (.not. file_exists) call system('mkdir '//trim(par_p))

call cpu_time(start)
!call chdir("/Users/julio/Desktop/hm-teste")
call getcwd(wd)
print *,"a criar PARs..."

do i_k=1,k
    !indice da realizacao
    if (i_k<10) format="(I1)"
    if (i_k>=10 .and. i_k<100) format="(I2)"
    if (i_k>=100 .and. i_k<1000) format="(I3)"
    if (i_k>=1000 .and. i_k<10000) format="(I4)"
    write (n_k,format) i_k
    !output da sim e ficheiro par
    output=trim(out_b)//trim(n_k)//'.out'
    par_sim=trim(par_b)//trim(n_k)//'.par'
    !nova semente
    call seeder(seed)
    !ciclo para editar ficheiro
    call chdir(wd)
    open(21,file=par_o,action='read') !21=par original
    call chdir(par_p)
    open(20,file=par_sim,action='write') !20=novo
    i_skip=1
    do
        if (i_skip==lin_out) then !caminho e ficheiro de saida
            read (21,*)
            write(20,*) trim(output)
        elseif (i_skip==lin_seed) then !semente
            read (21,*)
            write(20,*) seed
        else !copiar o resto
            read (21,'(A)',iostat=status) copy
            if (status<0) exit
            write(20,'(A)') copy
        end if
        i_skip=i_skip+1
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

end program dssparbaker_x2
