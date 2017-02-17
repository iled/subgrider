! probs :: Julio Caineta, 2010
! probs :: calcula as probabilidades a priori, a posteriori, interseccoes e condicionadas

program sg_probs

use sg_griders
use sg_utils
use sg_gridutils

implicit none
character(len=256)::ficheiro,base,path
logical::header
logical,allocatable::mapa_z(:)
integer::dg(3),id,nr_sims(2),nr_cen,i,j,a,b,ab,s,cl
real::nd,timer,k(2),start,finish,pa,pb,pab
real,allocatable::sims(:,:)
type(grid)::ref

header=.FALSE.
dg=(/ 334, 134, 93 /)
nd=-999
nr_sims=(/ 1, 50 /)
nr_cen=3
id=10
call chdir("/Users/julio/Desktop/sims")
open(9,file='probs_s2.txt',action='write')
write(9,*) "a = numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j"
write(9,*) "b = numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j"
write(9,*) "ab = numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j"

do s=1,3
    if (s==1) then
        call chdir("/Users/julio/Desktop/sims/DSS_G_s2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_G_s2.OUT'
        base='DSS_G_s2_sim'
    elseif (s==2) then
        call chdir("/Users/julio/Desktop/sims/DSS_M_s2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_M_s2.OUT'
        base='DSS_M_s2_sim'
    elseif (s==3) then
        call chdir("/Users/julio/Desktop/sims/DSS_P_s2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_P_s2.OUT'
        base='DSS_P_s2_sim'
    elseif (s==4) then
        call chdir("/Users/julio/Desktop/sims/DSS_M2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_M2.OUT'
        base='DSS_M2w'
    elseif (s==5) then
        call chdir("/Users/julio/Desktop/sims/DSS_M3")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_M3.OUT'
        base='DSS_M3w'
    elseif (s==6) then
        call chdir("/Users/julio/Desktop/sims/DSS_M4")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_M4.OUT'
        base='DSS_M4w'
    elseif (s==7) then
        call chdir("/Users/julio/Desktop/sims/DSS_P2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_P2.OUT'
        base='DSS_P2w'
    elseif (s==8) then
        call chdir("/Users/julio/Desktop/sims/DSS_P3")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_P3.OUT'
        base='DSS_P3w'
    elseif (s==9) then
        call chdir("/Users/julio/Desktop/sims/DSS_P4")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_P4.OUT'
        base='DSS_P4w'
    else
        close(9)
        stop
    end if
!print *,"ficheiro da grid de referencia"
!read *,ficheiro
print *,"a abrir grid de referencia ",trim(ficheiro)
call abre(ficheiro,header,dg,nd,ref,id,timer)
allocate(mapa_z(size(ref%val)))
print *,"operacao concluida em ",tempo(timer)
!print *,"base dos ficheiros das grids a comparar"
!read *,base
print *,"a abrir grids de comparacao ",trim(base)
call abre_sims(nr_sims,sims,dg,base,'.OUT',header,timer)
print *,"operacao concluida em ",tempo(timer)
do cl=1,3
    if (cl==1) k=(/ 24.2037, 40 /)
    if (cl==2) k=(/ 13.7452, 24.2037 /)
    if (cl==3) k=(/ 0.0, 13.7452 /)

print *,"intervalo (i,f) ",k
!read *,k
print *,"a percorrer grid de referencia..."
call cpu_time(start)
a=0 ! numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j
do i=1,size(ref%val)
    if (ref%val(i)>=k(1) .and. ref%val(i)<k(2)) then
        a=a+1
        mapa_z(i)=.TRUE.
    else
        mapa_z(i)=.FALSE.
!        if (sims
    end if
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a percorrer grids de comparacao..."
call cpu_time(start)
b=0 ! numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j
ab=0 ! numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j
do i=1,size(sims,1)
    do j=1,size(sims,2)
        if (sims(i,j)>=k(1) .and. sims(i,j)<k(2)) then
            b=b+1
            if (mapa_z(j)) ab=ab+1
        end if
    end do
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a escrever resultados..."
write(9,*) trim(path),"/",trim(ficheiro)
write(9,*) "classe ",k
write(9,*) ""
write(9,*) "blocos a: ",a
write(9,*) "blocos b: ",b
write(9,*) "blocos ab: ",ab
pa=a/real(size(ref%val)*nr_cen)
write(9,*) "P(i<=Z<f , S=s) = ",pa
write(9,*) "P(i<=Z<f | S=s) = ",pa/real(1/real(nr_cen))
pb=b/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z*<f , S=s) = ",pb
write(9,*) "P(i<=Z*<f | S=s) = ",pb/real(1/real(nr_cen))
pab=ab/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z<f , i<=Z*<f , S=s) = ",pab
write(9,*) "P(i<=Z*<f | i<=Z<f , S=s) = ",pab/pa
write(9,*) ""
end do
deallocate(sims,mapa_z)
end do

end program sg_probs
