! merge :: Julio Caineta, 2010
! merge :: junta mapas de probabilidade geradas pelo modulo sg_griders

program sg_merge

use sg_griders
use sg_utils

implicit none

integer::nr_sims,dg(3),i,fid
character(len=256)::base,newfile
real::timer,nd,start,finish
real,allocatable::sims(:,:)
logical::header
!real,allocatable,dimension(:)::merged

fid=10
print *,"numero de ficheiros"
read *,nr_sims
print *,"base"
read *,base
print *,"dimensoes da grid"
read *,dg
print *,"a carregar ficheiros..."
header=.TRUE.
call abre_sims(nr_sims,sims,dg,base,'.prn',header,timer)
print *,"operacao concluida em ",tempo(timer)
call novo(header,1,fid,newfile,-1,timer)
nd=-999
print *,"a unir grids..."
!allocate(merged(size(sims,2))
call cpu_time(start)
open (fid,file=newfile)
do i=1,size(sims,2)
    if (sims(1,i).ne.nd .and. sims(2,i).eq.nd .and. sims(3,i).eq.nd) then
        !merged(i)=sims(1,i)
        write (fid,*) sims(1,i)
    elseif (sims(1,i).eq.nd .and. sims(2,i).ne.nd .and. sims(3,i).eq.nd) then
        !merged(i)=sims(2,i)
        write (fid,*) sims(2,i)
    elseif (sims(1,i).eq.nd .and. sims(2,i).eq.nd .and. sims(3,i).ne.nd) then
        !merged(i)=sims(3,i)
        write (fid,*) sims(3,i)
    else
        print *,"erro: dois valores possiveis para a mesma celula"
        stop
    end if
end do
call cpu_time(finish)
timer=start-finish
print *,"operacao concluida em ",tempo(timer)

end program sg_merge
