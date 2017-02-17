! probs :: Julio Caineta, 2010
! probs :: calcula as probabilidades a priori, a posteriori, interseccoes e condicionadas

program sg_probs

use sg_griders
use sg_utils
use sg_gridutils

implicit none
character(len=256)::cen,sim,path_cen,path_dir,path_sim
character(len=1)::n_var
logical::header
logical,allocatable::mapa_z(:)
integer::dg(3),id,nr_sims(2),nr_cen,i,j,a,b(3),ab(3),s,cl,n,i_var
real::nd,timer,k(2),start,finish,pa,pb(3),pab(3)
real,allocatable::sims(:,:,:),w(:)
type(grid)::ref

header=.FALSE.
dg=(/ 334, 134, 93 /)
nd=-999
nr_sims=(/ 1, 50 /)
nr_cen=3
id=10
call chdir("/Users/julio/Desktop/sims")
open(9,file='probs_new.txt',action='write')
write(9,*) "a = numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j"
write(9,*) "b = numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j"
write(9,*) "ab = numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j"

path_cen="/Users/julio/Desktop/sims/cens"

do s=1,9
    if (s==1) then
        path_sim="/Users/julio/Desktop/sims/dss_w3_G"
        cen='dss_w3_G.OUT'
        sim='dss_w3_G'
    elseif (s==2) then
        path_sim="/Users/julio/Desktop/sims/dss_w3_M"
        cen='dss_w3_M.OUT'
        sim='dss_w3_M'
    elseif (s==3) then
        path_sim="/Users/julio/Desktop/sims/dss_w3_P"
        cen='dss_w3_P.OUT'
        sim='dss_w3_P'
    elseif (s==4) then
        path_sim="/Users/julio/Desktop/sims/dss_w5_G"
        cen='dss_w5_G.OUT'
        sim='dss_w5_G'
    elseif (s==5) then
        path_sim="/Users/julio/Desktop/sims/dss_w5_M"
        cen='dss_w5_M.OUT'
        sim='dss_w5_M'
    elseif (s==6) then
        path_sim="/Users/julio/Desktop/sims/dss_w5_P"
        cen='dss_w5_P.OUT'
        sim='dss_w5_P'
    elseif (s==7) then
        path_sim="/Users/julio/Desktop/sims/dss_w10_G"
        cen='dss_w10_G.OUT'
        sim='dss_w10_G'
    elseif (s==8) then
        path_sim="/Users/julio/Desktop/sims/dss_w10_M"
        cen='dss_w10_M.OUT'
        sim='dss_w10_M'
    elseif (s==9) then
        path_sim="/Users/julio/Desktop/sims/dss_w10_P"
        cen='dss_w10_P.OUT'
        sim='dss_w10_P'
    else
        close(9)
        stop
    end if
print *,"a abrir grid de referencia ",trim(cen)
call chdir(path_cen)
call abre(cen,header,dg,nd,ref,id,timer)
allocate(mapa_z(size(ref%val)))
print *,"operacao concluida em ",tempo(timer)

!carregar ficheiros de simulacao
print *,"a abrir grids de comparacao ",trim(sim)
call chdir(path_sim)
call abre_sims_amp(nr_sims,sims,dg,sim,'.OUT',header,timer)
print *,"operacao concluida em ",tempo(timer)
!definicao das classes para o calculo das probabilidades (valores do shah)
do cl=1,3
    if (cl==1) k=(/ 25.7612, 40.0/)
    if (cl==2) k=(/ 15.3437, 25.7612 /)
    if (cl==3) k=(/ 0.0, 15.3437 /)
print *,"intervalo (i,f) ",k

print *,"a percorrer grid de referencia..."
call cpu_time(start)
a=0 ! numero de blocos que simultaneamente cumprem Z in classe_i e S=cenario_j
do i=1,size(ref%val)
    if (ref%val(i)>=k(1) .and. ref%val(i)<k(2)) then
        a=a+1
        mapa_z(i)=.TRUE.
    else
        mapa_z(i)=.FALSE.
    end if
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a percorrer grids de comparacao..."
call cpu_time(start)
b=0 ! numero de blocos que simultaneamente cumprem Z* in classe_i e S=cenario_j
ab=0 ! numero de blocos que simultaneamente cumprem Z in classe_i, Z* in classe_i e S=cenario_j
do i_var=1,3
    if (i_var==1) n_var='G'
    if (i_var==2) n_var='M'
    if (i_var==3) n_var='P'
	do i=1,size(sims,1)
	    do j=1,size(sims,2)
	        if (sims(i,j,i_var)>=k(1) .and. sims(i,j,i_var)<k(2)) then
	            b(i_var)=b(i_var)+1
	            if (mapa_z(j)) ab(i_var)=ab(i_var)+1
	        end if
	    end do
	end do
end do
call cpu_time(finish)
timer=finish-start
print *,"operacao concluida em ",tempo(timer)
print *,"a escrever resultados..."
write(9,*) trim(path_cen),"/",trim(cen)
write(9,*) "simulacao com cenario ",n_var
write(9,*) "classe ",k
write(9,*) ""
write(9,*) "blocos a: ",a
write(9,*) "P(i<=Z<f , S=s) = ",pa
write(9,*) "P(i<=Z<f | S=s) = ",pa/real(1/real(nr_cen))
write(9,*) ""
do i_var=1,3
    if (i_var==1) n_var='G'
    if (i_var==2) n_var='M'
    if (i_var==3) n_var='P'
write(9,*) "blocos b(",n_var,"): ",b(i_var)
write(9,*) "blocos ab(",n_var,"): ",ab(i_var)
pa=a/real(size(ref%val)*nr_cen)
pb(i_var)=b(i_var)/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z*<f , S=",n_var,") = ",pb(i_var)
write(9,*) "P(i<=Z*<f | S=",n_var,") = ",pb(i_var)/real(1/real(nr_cen))
pab(i_var)=ab(i_var)/real(size(sims,2)*nr_cen*(nr_sims(2)-nr_sims(1)+1))
write(9,*) "P(i<=Z<f , i<=Z*<f , S=",n_var,") = ",pab(i_var)
write(9,*) "P(i<=Z*<f | i<=Z<f , S=",n_var,") = ",pab(i_var)/pa
write(9,*) ""
end do
end do
deallocate(sims,mapa_z)
end do

contains

! carrega ficheiros de simulacoes para uma matriz, adaptado para as novas simulacoes, e.g., dss_w3_G_1_G.out
subroutine abre_sims_amp(nr_sims,sims,dg,base,ext,header,timer)
integer,intent(in)::nr_sims(2),dg(3)
logical,intent(in)::header
real,allocatable,intent(out)::sims(:,:,:)
character(len=256),intent(in)::base
character(len=4),intent(in)::ext
character(len=256)::sim
character(len=4)::simN,format
character(len=1)::n_var
type(grid)::gridsim
integer::i,id,i_var
real,intent(out)::timer
real::start,finish,t
call cpu_time(start)
! linhas=numero de simulacoes
! colunas=numero de blocos em cada simulacao
! 3d=numero de variogramas
allocate(sims(nr_sims(2)-nr_sims(1)+1,product(dg),3))
id=20
do i_var=1,3
    if (i_var==1) n_var='G'
    if (i_var==2) n_var='M'
    if (i_var==3) n_var='P'
    do i=nr_sims(1),nr_sims(2)
	    if (i<10) format="(I1)"
	    if (i>=10 .and. i<100) format="(I2)"
	    if (i>=100 .and. i<1000) format="(I3)"
	    write (simN,format) i
	    sim=trim(base)//'_'//trim(simN)//'_'//n_var//trim(ext)
        call checkfile(sim)
        call abre(sim,header,dg,-999.0,gridsim,id,t)
        sims(i,:,i_var)=gridsim%val
        id=id+1
    end do
end do
call cpu_time(finish)
timer=finish-start
end subroutine abre_sims_amp

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

end program sg_probs
