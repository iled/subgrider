program subgrider

use sg_utils
use sg_griders
use sg_gridutils

implicit none

! declaracao de variaveis
integer::i,nvar,fid,batch,nr_sims,nr_p,zef
integer,dimension(3)::dg,pa,pb,dg_2d
character(len=256)::ficheiro,output,pocos,base,mp_bend
character(len=1)::load,mp_p
logical::mp,header
real,allocatable::sims(:,:)
integer,allocatable::px(:),py(:)
real::P_x(2),P1,timer,nd,p,q,op,p_sim(2)
type(grid)::res,hor

!ciclo principal
do
batch=-1
print *,""
print *,"----|| s u b g r i d e r  v0.2.1 <<jc 2010>> ||----"
print *,""
print *,"escolher uma opcao"
print *,""
print *,"0 ----- carregar dados"
print *,"  0.1 - abrir ficheiro"
print *,"  0.2 - abrir ficheiro anterior"
print *,"1 ----- poco vertical "
print *,"  1.1 - a partir das coordenadas (x,y)"
print *,"  1.2 - a partir de um ficheiro com coordenadas (x,y)"
print *,"2 - obter uma grid a partir um volume (cubo, paralelipipedo) da grid inicial"
print *,"3 - criar uma mascara (data=1, no data=0)"
print *,"4 ----- update bayesiano"
print *,"  4.1 - calcular percentil"
print *,"  4.2 - calcular quantil"
print *,"  4.3 - calcular verosimilhanca"
print *,"5 - planificacao"
print *,"  5.1 - planificar uma grid a partir de um mapa de planificacao"
print *,"9 - sair"
read *,op
if (op==0.1) then
    print *,"abrir ficheiro"
    print *,""
    fid=10
    print *,"ficheiro a carregar (caminho/nome)"
    read *,ficheiro
    call checkfile(ficheiro)
    call header_ask("r",header,nvar)
    print *,"dimensoes da grid: x y z"
    read (*,*) dg
    print *,"valor correspondente a 'no data'"
    read *,nd
    open (9,file='subgrider.cfg',action='write')
    write (9,*) ficheiro,header
    write (9,*) dg,nd
    close (9)
    print *,"a carregar o ficheiro ",trim(ficheiro),"..."
    call abre(ficheiro,header,dg,nd,res,fid,timer)
    print *,"ficheiro carregado em ",tempo(timer)
elseif (op==0.2) then
    print *,"abrir ficheiro anterior"
    print *,""
    fid=10
    call checkfile("subgrider.cfg")
    open (9,file='subgrider.cfg',action='read')
    read (9,*) ficheiro,header
    read (9,*) dg,nd
    print *,"a carregar o ficheiro ",trim(ficheiro),"..."
    call abre(ficheiro,header,dg,nd,res,fid,timer)
    print *,"ficheiro carregado em ",tempo(timer)
elseif (op==1.1) then
    print *,"extrair poco vertical a partir das coordenadas (x,y)"
    print *,""
    print *,"coordenadas do poco (x,y)"
    read *,pa(1),pa(2)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a furar o poco..."
    call pp(pa(1),pa(2),res,fid,output,header,batch,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==1.2) then
    print *,"extrair pocos verticais a partir de um ficheiro com coordenadas (x,y)"
    print *,""
    print *,"ficheiro com coordenadas (caminho/nome)"
    read *,pocos
    call checkfile(pocos)
    print *,"a ler ficheiro ",trim(pocos),"..."
    call abre_pocos(pocos,px,py,timer)
    print *,trim(pocos)," carregado em ",tempo(timer)
    call header_ask("w",header,nvar)
    batch=0
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a furar os pocos..."
    do i=1,size(px)
        call pp(px(i),py(i),res,fid,output,header,batch,timer)
    end do
    print *,"operacao concluida em ",tempo(timer) !tempo devolvido é do último poço
    batch=-1
elseif (op==2) then
    print *,"obter uma grid a partir um volume (cubo, paralelipipedo) da grid inicial"
    print *,""
    print *,"coordenadas x (min,max)"
    read *,pa(1),pb(1)
    print *,"coordenadas y (min,max)"
    read *,pa(2),pb(2)
    print *,"coordenadas z (min,max)"
    read *,pa(3),pb(3)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a criar subgrid..."
    call subgrid(pa,pb,res,fid,output,header,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==3) then
    print *,"criar uma mascara (data=1, no data=0)"
    print *,""
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a criar mascara..."
    call mask(1,0,header,res,fid,output,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==4.3) then
    print *,"update bayesiano - calcular verosimilhanca"
    print *,""
    print *,"numero de simulacoes realizadas (ficheiros *_simX.OUT)"
    read *,nr_sims
! <alt>
    print *,"base"
    read *,base
! <\alt>
!    base=ficheiro(1:len_trim(ficheiro)-4)//'_sim'
    call header_ask("r",header,nvar)
    print *,"a carregar grids simuladas..."
    call abre_sims(nr_sims,sims,dg,base,'.OUT',header,timer)
    print *,nr_sims,"grids carregadas em ",tempo(timer)
    do
        print *,"produzir grid com mapa de probabilidades? (s/n)"
        read *,mp_p
        if (mp_p=="s" .or. mp_p=="S") then
            mp=.true.
            exit
        elseif (mp_p=="n" .or. mp_p=="N") then
            mp=.false.
            exit
        else
            print *,"erro - input invalido."
        end if
    end do
    if (mp) then
        call header_ask("w",header,nvar)
        call novo(header,nvar,fid,output,batch,timer)
    end if
    !print *,"inserir valor x (P[X=x])"
    !read *,P_x
    print *,"inserir limites do intervalo (xi, xf) (P[xi <= X < xf])"
    read *,P_x
    print *,"a calcular probabilidades..."
    call likely(res,sims,P_x,P1,mp,fid,output,header,timer)
    print *,"a probabilidade e' ",trim(rtochar(P1))
    print *,"operacao concluida em ",tempo(timer)
    ! colocado aqui temporariamente, incorporar como batch do percentil
    print *,""
    print *,"percentis das simulacoes"
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"inserir intervalo de percentis (pi, pf)"
    read *,p_sim
    print *,"a calcular percentis..."
    call percentil_sims(sims,p_sim,P_x,fid,header,output,timer) !verificar P_x(1)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==4.1) then
    print *,"update bayesiano - calcular percentil"
    print *,""
    print *,"inserir percentil a calcular"
    read *,p
    call percentil(res,p,q,timer)
    print '("o percentil ",f4.2," e ",a)',p,trim(rtochar(q))
    print *,"operacao concluida em ",tempo(timer)
elseif (op==4.2) then
    print *,"update bayesiano - calcular quantil"
    print *,""
    print *,"inserir valor"
    read *,q
    call percentil_updt(res,q,p,timer)
    print *,"o valor ",trim(rtochar(q))," corresponde ao percentil ",trim(rtochar(p))
    print *,"operacao concluida em ",tempo(timer)
elseif (op==5.1) then
    print *,"planificar uma grid a partir de um mapa de planificacao"
    print *,"atencao: a grid inicial (alocada) vai ser alterada"
    print *,""
    print *,"dimensao maxima em z a escrever"
    read *,zef
    print *,"mapa de planificacao a abrir: caminho/nome"
    read *,mp_bend
    call checkfile(mp_bend)
    call header_ask("r",header,nvar)
    print *,"dimensoes do mapa: x y"
    read (*,*) dg_2d(1),dg_2d(2)
    dg_2d(3)=1
    print *,"a ler mapa de planificacao..."
    call abre(mp_bend,header,dg_2d,nd,hor,fid+20,timer) !nd
    print *,"mapa ",trim(mp_bend)," lido em ",tempo(timer)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a planificar grid..."
    call bender(dg(1),dg(2),dg(3),-999.25,fid,output,header,zef,hor,res,timer)
    print *,"operacao concluida em ",tempo(timer)
elseif (op==9) then
    print *,"programa terminado"
    stop
else
    print *,"erro - input invalido."
end if
end do

! funcoes e subrotinas

contains

! devolve o tempo, em string, gasto pela ultima subrotina executada
function tempo(timer) result(time)
real,intent(in)::timer
real::te
!character(10),allocatable,dimension(:)::time
character(32)::time
character(len=6)::format,uni
if (timer<1) then
    te=timer*1000
    format="(f7.3)"
    uni=" ms"
elseif (timer>=1 .and. timer<60) then
    te=timer
    format="(f6.3)"
    uni=" s"
elseif (timer>=60 .and. timer<3600) then
    te=timer/60
    format="(f6.3)"
    uni=" min"
else
    te=timer/3600
    format="(f7.3)"
    uni=" h"
end if
write (time,format) te
time=time(1:len_trim(time))//trim(uni)
end function tempo

end program
