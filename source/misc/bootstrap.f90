program bootstrap

use sg_utils
use sg_griders
use sg_gridutils

implicit none
integer::k,dg(3),i,j,batch,fid,s,cut,did
integer,allocatable::px(:),py(:),wells(:,:,:)
real::timer,start,finish,nd
character(len=256)::pocos,output,ficheiro,debug
character(len=1)::amp
logical::header
type(grid)::res
character(len=4)::kn,format

header=.FALSE.
fid=10
dg=(/ 334, 134, 93 /)
ficheiro='dss_w'
nd=-999.0

print *,"0. reamostragem aleatoria de K pocos"
print *,""
read *,k
!print *,"ficheiro com coordenadas (caminho/nome)"
!read *,pocos
do s=1,9
    call chdir("/Users/julio/Desktop/sims")
    if (s==1) then
        pocos='pocos3.cfg'
        output='wells_w3_'
        call checkfile(pocos)
	    print *,"1. a ler pocos ",trim(pocos),fid
	    call abre_pocos(pocos,px,py,timer)
	    print *,"2. ",trim(pocos)," carregado em ",tempo(timer)
	    amp='G'
	    ficheiro=ficheiro(1:5)//'3_'//amp//'.out'
    elseif (s==2) then
        output='wells_w3_'
        amp='M'
        ficheiro=ficheiro(1:7)//amp//'.out'
    elseif (s==3) then
        output='wells_w3_'
        amp='P'
        ficheiro=ficheiro(1:7)//amp//'.out'
    elseif (s==4) then
        pocos='pocos5.cfg'
        output='wells_w5_'
        print *,"1. a ler pocos ",trim(pocos),fid
        call abre_pocos(pocos,px,py,timer)
        print *,"2. ",trim(pocos)," carregado em ",tempo(timer)
        amp='G'
        ficheiro=ficheiro(1:5)//'5_'//amp//'.out'
    elseif (s==5) then
        output='wells_w5_'
        amp='M'
        ficheiro=ficheiro(1:7)//amp//'.out'
    elseif (s==6) then
        output='wells_w5_'
        amp='P'
        ficheiro=ficheiro(1:7)//amp//'.out'
    elseif (s==7) then
        pocos='pocos10.cfg'
        output='wells_w10_'
        print *,"1. a ler pocos ",trim(pocos),fid
        call abre_pocos(pocos,px,py,timer)
        print *,"2. ",trim(pocos)," carregado em ",tempo(timer)
        amp='G'
        ficheiro=ficheiro(1:5)//'10_'//amp//'.out'
    elseif (s==8) then
        output='wells_w10_'
        amp='M'
        ficheiro=ficheiro(1:8)//amp//'.out'
    elseif (s==9) then
        output='wells_w10_'
        amp='P'
        ficheiro=ficheiro(1:8)//amp//'.out'
    end if
    call chdir("/Users/julio/Desktop/sims/cens")
    print *,"3. a carregar cenario ",trim(ficheiro),fid
    call abre(ficheiro,header,dg,nd,res,fid,timer)
    print *,"4. cenario carregado em ",tempo(timer)
    fid=fid+1
    !call header_ask("w",header,nvar)
    batch=0
    !call novo(header,nvar,fid,output,batch,timer)
    print *,"5. a gerar posicoes..."
    call rand_well(k,dg,px,py,wells,timer)
    print *,"6. posicoes geradas em ",tempo(timer)
    output=trim(output)//amp
    print *,"7. a furar os pocos ",trim(output)
    call cpu_time(start)
    call chdir("/Users/julio/Desktop/sims/bootstrap")
    call chdir(output)
    debug=trim(output)//'.txt'
    output=trim(output)//"_"
    cut=len_trim(output)
    open(fid,file=debug,action='write')
    did=fid
    fid=fid+1
    do j=1,k
        if (j<10) format="(I1)"
        if (j>=10 .and. j<100) format="(I2)"
        if (j>=100 .and. j<1000) format="(I3)"
        if (j>=1000 .and. j<10000) format="(I4)"
        write (kn,format) j
        output=output(1:cut)//trim(kn)//'.prn'
        open(fid,file=output,action='write')
        do i=1,size(px)
            call pp(wells(i,1,k),wells(i,2,k),res,fid,output,header,batch,timer)
        end do
        write (did,*) wells(1,1,j),wells(1,2,j)
        batch=0
        fid=fid+1
    end do
    call cpu_time(finish)
    timer=finish-start
    print *,"8. pocos furados em ",tempo(timer)
end do

end program bootstrap
