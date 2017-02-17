module pointers

implicit none

public :: lista_nova, lista_comp, lista_adf, lista_remf, lista_array

private

! tipo de dados lista bidireccional
type, public :: lista
    private
    type(no), pointer :: prim,ult
    integer :: compr
end type lista

type, private :: no
    integer :: pos
    real :: val
    type(no), pointer :: ant, seg
end type no

contains

! cria uma nova lista
function lista_nova() result(list)
type(lista) :: list
list%prim=>null()
list%ult=>null()
list%compr=0
end function lista_nova

! devolve o comprimento da lista
function lista_comp(list) result(c)
type(lista), intent(in) :: list
integer :: c
c=list%compr
end function lista_comp

! adiciona uma entrada ao final da lista
subroutine lista_adf(list,x)
type(lista), intent(inout) :: list
real, intent(in) :: x
type(no), pointer :: temp
allocate(temp)
temp%val=x
temp%seg=>null()
if (associated(list%ult)) then
    temp%pos=list%ult%pos+1
    temp%ant=>list%ult
    list%ult%seg=>temp
else
    temp%pos=1
    temp%ant=>null()
    list%prim=>temp
end if
list%ult=>temp
list%compr=list%compr+1
end subroutine lista_adf

! remove n entradas a partir do final da lista
subroutine lista_remf(list,n)
type(lista), intent(inout) :: list
integer, intent(in) :: n
type(no), pointer :: temp
integer :: i
do i=1,n
	if (associated(list%ult)) then
        allocate(temp)
	    temp=>list%ult
	    if (associated(list%ult%ant)) then
    	    list%ult=>list%ult%ant
            list%ult%seg=>null()
        else
            list%ult=>null()
            list%prim=>null()
        end if
	    deallocate(temp)
	    list%compr=list%compr-1
	else
	    exit
	end if
end do
end subroutine lista_remf

! converte a lista num vector
subroutine lista_array(list,arr)
type(lista), intent(inout) :: list
real, dimension(:), allocatable, intent(out) :: arr
type(no), pointer :: temp
integer :: i
real::t
if (associated(list%ult)) then
	t=list%ult%val
	allocate(temp)
	allocate(arr(list%compr))
	temp=>list%prim
	do i=1,list%compr
	    arr(i)=temp%val
	    temp=>temp%seg
	end do
else
    print *,"erro: lista vazia (resultado=0)"
    allocate(arr(1))
    arr=0.0
end if
end subroutine lista_array

end module pointers
