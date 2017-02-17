! subgrider :: Julio Caineta, 2010
! sg_main :: interface

program subgrider

use sg_utils
use sg_griders
use sg_gridutils

implicit none

! declaracao de variaveis
integer::i,nvar,fid,batch,nr_sims(2),nr_p,zef,ask_simmedvar
integer,dimension(3)::dg,pa,pb,dg_2d
character(len=256)::ficheiro,output,pocos,base,mp_bend
character(len=1)::load,mp_p
logical::mp,header,do_med,do_var,did_load,did_loadsims
real,allocatable::sims(:,:)
integer,allocatable::px(:),py(:)
real::P_x(2),P1,timer,nd,p,q,op,p_sim(2)
type(grid)::res,hor

did_load=.FALSE.
did_loadsims=.FALSE.
!ciclo principal
do
batch=-1
print *,""
print *,"----|| s u b g r i d e r  v0.2.3 <<jc 2010>> ||----"
print *,""
print *,"escolher uma opcao"
print *,""
print *,"0 ----- carregar dados"
print *,"  0.1 - abrir ficheiro"
print *,"  0.2 - abrir ficheiro anterior"
print *,"  0.3 - abrir varios ficheiros"
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
print *,"6 - analise de simulacoes"
print *,"  6.1 - media e variancia das simulacoes"
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
    did_load=.TRUE.
    call wait()
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
    did_load=.TRUE.
    call wait()
elseif (op==0.3) then
    print *,"abrir varios ficheiros numerados de A a B"
    print *,"inserir intervalo (A,B)"
    read *,nr_sims
    print *,"base"
    read *,base
    call header_ask("r",header,nvar)
    if (.not. did_load) then
        print *,"dimensoes das grids: x y z"
        read (*,*) dg
    end if
    print *,"a carregar ficheiros..."
    call abre_sims(nr_sims,sims,dg,base,'.OUT',header,timer)
    print *,nr_sims(2)-nr_sims(1)+1,"ficheiros carregados em ",tempo(timer)
    did_loadsims=.TRUE.
    call wait()
elseif (op==1.1 .and. did_load) then
    print *,"extrair poco vertical a partir das coordenadas (x,y)"
    print *,""
    print *,"coordenadas do poco (x,y)"
    read *,pa(1),pa(2)
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a furar o poco..."
    call pp(pa(1),pa(2),res,fid,output,header,batch,timer)
    print *,"operacao concluida em ",tempo(timer)
    call wait()
elseif (op==1.2 .and. did_load) then
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
    print *,"operacao concluida em ",tempo(timer) !tempo devolvido do ultimo poco
    batch=-1
    call wait()
elseif (op==2 .and. did_load) then
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
    call wait()
elseif (op==3 .and. did_load) then
    print *,"criar uma mascara (data=1, no data=0)"
    print *,""
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a criar mascara..."
    call mask(1,0,header,res,fid,output,timer)
    print *,"operacao concluida em ",tempo(timer)
    call wait()
elseif (op==4.3 .and. did_load .and. did_loadsims) then
    print *,"update bayesiano - calcular verosimilhanca"
    print *,""
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
    do
	    print *,"produzir ficheiro com os percentis das simulacoes? (s/n)"
	    read *,mp_p
	    if (mp_p=="s" .or. mp_p=="S") then
		    call header_ask("w",header,nvar)
		    call novo(header,nvar,fid,output,batch,timer)
		    print *,"inserir intervalo de percentis (pi, pf)"
		    read *,p_sim
		    print *,"a calcular percentis..."
		    !call percentil_sims(sims,p_sim,P_x,fid,header,output,timer) !verificar P_x(1)
		    call percentil_sims(sims,p_sim,P_x,fid,header,output,timer) !verificar P_x(1)
		    print *,"operacao concluida em ",tempo(timer)
	    elseif (mp_p=="n" .or. mp_p=="N") then
	        exit
	    else
	        print *,"erro - input invalido."
	    end if
    end do
    call wait()
elseif (op==4.1 .and. did_load) then
    print *,"update bayesiano - calcular percentil"
    print *,""
    print *,"inserir percentil a calcular"
    read *,p
    call percentil(res,p,q,timer)
    print '("o percentil ",f4.2," e ",a)',p,trim(rtochar(q))
    print *,"operacao concluida em ",tempo(timer)
    call wait()
elseif (op==4.2 .and. did_load) then
    print *,"update bayesiano - calcular quantil"
    print *,""
    print *,"inserir valor"
    read *,q
    call percentil_updt(res,q,p,timer)
    print *,"o valor ",trim(rtochar(q))," corresponde ao percentil ",trim(rtochar(p))
    print *,"operacao concluida em ",tempo(timer)
    call wait()
elseif (op==5.1 .and. did_load) then
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
    call wait()
elseif (op==6.1 .and. did_loadsims) then
    print *,"analise de simulacoes"
    print *,"media e variancia das simulacoes"
    print *,""
    print *,"calcular media (1), variancia (2) ou ambas (3)?"
    read *,ask_simmedvar
    do
        if (ask_simmedvar==1) then
            do_med=.TRUE.
            do_var=.FALSE.
            exit
        elseif (ask_simmedvar==2) then
            do_med=.FALSE.
            do_var=.TRUE.
            exit
        elseif (ask_simmedvar==3) then
            do_med=.TRUE.
            do_var=.TRUE.
            exit
        else
            print *,"erro - input invalido"
        end if
    end do
    call header_ask("w",header,nvar)
    call novo(header,nvar,fid,output,batch,timer)
    print *,"a calcular..."
    call simmedvar(do_med,do_var,sims,fid,output,header,timer)
    print *,"operacao concluida em ",tempo(timer)
    call wait()
elseif (op==9) then
    print *,"programa terminado"
    stop
else
    print *,"erro - input invalido."
end if
end do

contains

subroutine wait()
print *,""
print *,"primir [ENTER] para voltar ao menu principal"
read (*,'()')
end subroutine wait

end program
