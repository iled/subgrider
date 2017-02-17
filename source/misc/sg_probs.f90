! probs :: Julio Caineta, 2010
! probs :: calcula as probabilidades a priori, a posteriori, interseccoes e condicionadas

program sg_probs

use sg_griders
use sg_utils
use sg_gridutils

implicit none
character(len=256)::ficheiro,base,path
logical::header
integer,allocatable::mapa_z(:)
integer::dg(3),id,nr_sims(2),nr_cen,i,j,a1,a2,a3,b1,b2,b3,ab(3,3),s,cl
real::nd,timer,k(3,2),start,finish,pa(3),pb(3),pab(3,3)
real,allocatable::sims(:,:)
type(grid)::ref

header=.FALSE.
dg=(/ 334, 134, 93 /)
nd=-999
nr_sims=(/ 1, 50 /)
nr_cen=3
id=10
call chdir("/Users/julio/Desktop/sims")
open(9,file='probs3_s1.txt',action='write')
write(9,*) "a = numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j"
write(9,*) "b = numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j"
write(9,*) "ab = numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j"

do s=1,9
!    if (s==1) then
!        call chdir("/Users/julio/Desktop/sims/DSS_G_s2")
!        call getcwd(path)
!        print *,trim(path)
!        ficheiro='DSS_G_s2.OUT'
!        base='DSS_G_s2_sim'
!    elseif (s==2) then
!        call chdir("/Users/julio/Desktop/sims/DSS_M_s2")
!        call getcwd(path)
!        print *,trim(path)
!        ficheiro='DSS_M_s2.OUT'
!        base='DSS_M_s2_sim'
!    elseif (s==3) then
!        call chdir("/Users/julio/Desktop/sims/DSS_P_s2")
!        call getcwd(path)
!        print *,trim(path)
!        ficheiro='DSS_P_s2.OUT'
!        base='DSS_P_s2_sim'
    if (s==1) then
        call chdir("/Users/julio/Desktop/sims/DSS_G2")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_G2.OUT'
        base='DSS_G2w'
    elseif (s==2) then
        call chdir("/Users/julio/Desktop/sims/DSS_G3")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_G3.OUT'
        base='DSS_G3w'
    elseif (s==3) then
        call chdir("/Users/julio/Desktop/sims/DSS_G4")
        call getcwd(path)
        print *,trim(path)
        ficheiro='DSS_G4.OUT'
        base='DSS_G4w'
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
k(1,:)=(/ 24.3387, 40 /)
k(2,:)=(/ 13.9542, 24.3387 /)
k(3,:)=(/ 0.0, 13.9542 /)
!k(1,:)=(/ 24.2037, 40 /)
!k(2,:)=(/ 13.7452, 24.2037 /)
!k(3,:)=(/ 0.0, 13.7452 /)
write(9,*) "classes: ",k
print *,"a percorrer grid de referencia..."
call cpu_time(start)
a1=0 ! numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j
a2=0
a3=0
do i=1,size(ref%val)
    if (ref%val(i)>=k(1,1) .and. ref%val(i)<k(1,2)) then
        a1=a1+1
        mapa_z(i)=1
    elseif (ref%val(i)>=k(2,1) .and. ref%val(i)<k(2,2)) then
        a2=a2+1
        mapa_z(i)=2
    elseif (ref%val(i)>=k(3,1) .and. ref%val(i)<k(3,2)) then
        a3=a3+1
        mapa_z(i)=3
    else
        print *,"ooops!"
        stop
    end if
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a percorrer grids de comparacao..."
call cpu_time(start)
b1=0 ! numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j
b2=0
b3=0
!ab=0 ! numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j
ab=0
do i=1,size(sims,1)
    do j=1,size(sims,2)
        if (sims(i,j)>=k(1,1) .and. sims(i,j)<k(1,2)) then
            b1=b1+1
            if (mapa_z(j)==1) ab(1,1)=ab(1,1)+1
            if (mapa_z(j)==2) ab(2,1)=ab(2,1)+1
            if (mapa_z(j)==3) ab(3,1)=ab(3,1)+1
        elseif (sims(i,j)>=k(2,1) .and. sims(i,j)<k(2,2)) then
            b2=b2+1
            if (mapa_z(j)==1) ab(1,2)=ab(1,2)+1
            if (mapa_z(j)==2) ab(2,2)=ab(2,2)+1
            if (mapa_z(j)==3) ab(3,2)=ab(3,2)+1
        elseif (sims(i,j)>=k(3,1) .and. sims(i,j)<k(3,2)) then
            b3=b3+1
            if (mapa_z(j)==1) ab(1,3)=ab(1,3)+1
            if (mapa_z(j)==2) ab(2,3)=ab(2,3)+1
            if (mapa_z(j)==3) ab(3,3)=ab(3,3)+1
        end if
    end do
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a escrever resultados..."
write(9,*) trim(path),"/",trim(ficheiro)
write(9,*) "blocos a: ",a1,a2,a3
write(9,*) "blocos b: ",b1,b2,b3
write(9,*) "blocos a_b1: ",ab(:,1)
write(9,*) "blocos a_b2: ",ab(:,2)
write(9,*) "blocos a_b2: ",ab(:,3)
pa(1)=a1/real(size(ref%val)*nr_cen)
pa(2)=a2/real(size(ref%val)*nr_cen)
pa(3)=a3/real(size(ref%val)*nr_cen)
write(9,*) "P(i<=Z<f , S=s) = ",pa
write(9,*) "P(i<=Z<f | S=s) = ",pa/real(1/real(nr_cen))
pb(1)=b1/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
pb(2)=b2/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
pb(3)=b3/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z*<f , S=s) = ",pb
write(9,*) "P(i<=Z*<f | S=s) = ",pb/real(1/real(nr_cen))
pab=ab/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z<f , i<=Z*<f , S=s) = ",pab
write(9,*) "P(i<=Z*<f | i<=Z<f , S=s) = ",pab(1,:)/pa(1)
write(9,*) "P(i<=Z*<f | i<=Z<f , S=s) = ",pab(2,:)/pa(2)
write(9,*) "P(i<=Z*<f | i<=Z<f , S=s) = ",pab(3,:)/pa(3)
write(9,*) ""
deallocate(sims,mapa_z)
end do

end program sg_probs
