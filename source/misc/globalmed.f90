program globalmed

implicit none
real::med,start,finish
integer::dg(3),i
character(len=256)::output

print *,"dimensoes da grid (x, y, z)"
read *,dg
print *,"media global a impor"
read *,med
print *,"ficheiro de saida"
read *,output

call cpu_time(start)
open(9,file=output,action='write')
do i=1,product(dg)
    write(9,*) med
end do
call cpu_time(finish)
print *,finish-start

end program globalmed
