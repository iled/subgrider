program randwell

use sg_utils, only:abre_pocos,rand_well,tempo

implicit none

integer::k,dg(3),i,j
integer,allocatable::px(:),py(:),wells(:,:,:)
real::timer,start,finish
character(len=256)::file

file='pocos.txt'
!file='test.txt'
k=100
dg=(/ 334, 134, 93 /)
!dg=(/ 10, 20, 1 /)
call abre_pocos(file,px,py,timer)
print *,"pocos carregados em",tempo(timer)
print *,"p ",px,py
call rand_well(k,dg,px,py,wells,timer)
print *,"posicoes geradas em ",tempo(timer)
open (9,file='randwell.txt',action='write')
open (10,file='randwell_fixo.txt',action='write')
call cpu_time(start)
do j=1,k
    do i=1,size(px)
        write (9,*) wells(i,1,j),wells(i,2,j)
    end do
    write (10,*) wells(1,1,j),wells(1,2,j)
end do
call cpu_time(finish)
print *,"resultados escritos em ",tempo(timer)

end program randwell
