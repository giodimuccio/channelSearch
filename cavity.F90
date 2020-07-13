! ----------------------------------------------------------
! Read and Convert dx volmap into FORTRAN matrix
! ----------------------------------------------------------

module volmap
real, dimension(:,:,:), allocatable  :: map,map2
integer :: xn,yn,zn ! nn number of map2 white points
real :: x0,y0,z0,a1,a2,a3,b1,b2,b3,c1,c2,c3
end module volmap

! ----------------------------------------------------------

subroutine dxtomatrix(map_dx_name)

use volmap

implicit none

real, dimension(:,:), allocatable :: x
real, dimension(:), allocatable :: M
integer :: n,i,j,k,error,modu
character(len=150) :: map_dx_name
character(len=20) :: s

!write(*,*) map_dx_name

open (unit=99, file=map_dx_name, status='old', action='read')

n=0
DO
 READ(99,*, iostat = error)
 IF (error == -1) EXIT
 n = n + 1
END DO
!write(*,*) "n = ",n
close(99)
allocate(x(n-10,3))
allocate(M((n-10)*3))
open (unit=99, file=map_dx_name, status='old', action='read')

! ------------- LEGGI FILE DX ------------------
read(99,*)  !comment
read(99,*)  s,s,s,s,s,xn,yn,zn
read(99,*)  s,x0,y0,z0   ! origin
read(99,*)  s,a1,a2,a3   ! vector basis 1 a
read(99,*)  s,b1,b2,b3   ! vector basis 2 b
read(99,*)  s,c1,c2,c3   ! vector basis 3 c
read(99,*)  !comment
read(99,*)  !comment

allocate(map(xn,yn,zn))
allocate(map2(xn,yn,zn))

modu=mod(xn*yn*zn,3)
k=1
do i=1,n-8-5-1 !8 head, 5 tail comment line, - last line
   read(99,*) x(i,:)
   do j=1,3
     M(k)=x(i,j)
     k=k+1
!     write(*,*) k, M(k)
   enddo
enddo

! last line
if (modu==0) then
   read(99,*) x(i,:)
   do j=1,3
     M(k)=x(i,j)
     k=k+1
   enddo
else if (modu==1) then
   read(99,*) M(k)
else if (modu==2) then
   read(99,*) M(k),M(k+1)
endif

!write(*,*) k
close(99)
! -------- CONVERTI in MATRICE 3D ----------------
n=1
do i=1,xn
 do j=1,yn
  do k=1,zn
   map(i,j,k)=M(n)
   n=n+1
!   write(*,*) M(i)
  enddo
 enddo
enddo

end subroutine dxtomatrix

! ----------------------------------------------------------
! ----------------------------------------------------------
! ----------------------------------------------------------

! ---------------------------------------------------------------------
! Search cavity, starting from the point (x0,y0,top/bottom)
! eps, r0: filters' parameters
! ---------------------------------------------------------------------

subroutine cavity(x0,y0,eps,r0)

use volmap
 
implicit none
integer x,y,x0,y0,z0,diff,i
real :: map0(xn,yn,zn),map_su(xn,yn,zn),map_giu(xn,yn,zn),eps,r0
 
map0=map
i=0
diff=1

do while(abs(diff)>0)
 map=map0

! Map su
 z0=8
 x=x0
 y=y0
 call initial_point(x,y,z0,map)
 call search_cavity(map,map2,x,y,z0,eps,+1)
 map_su=map2
! Map giù
 z0=zn-7
 map=map0
 call search_cavity(map,map2,x0,y0,z0,eps,-1)
 map_giu=map2
! call intersezione(map_su,map_giu,map2)

! Applica Filtri
 if (r0>0.0) then
  call filtra_raggio(map2,r0,x0,y0)
 endif
 if (eps>0.0) then
  call filtra_densita(map2,eps)
 endif

 diff=int(sum(map0-map2))
 map0=map2
 i=i+1
! print *, diff,i
end do
end subroutine cavity

! ---------------------------------------------------
! ------  READ PARA ---------------------------------
! ---------------------------------------------------
subroutine read_para_cavity(x0,y0,eps,r0,sel_out,zh)
integer :: x0,y0,sel_out,zh
real    :: eps,r0

open (unit=678, file='cavity.para', status='old', action='read')

read(678,*) ! comment
read(678,*) x0,y0
read(678,*) ! comment
read(678,*) eps,r0
read(678,*) ! comment
read(678,*) sel_out,zh

end subroutine read_para_cavity

! ---------------------------------------------------
! ------ INITIAL POINT TEST -------------------------
! ---------------------------------------------------
subroutine initial_point(x0,y0,z0,m)

use volmap
 
implicit none
real, dimension(xn,yn,zn)  :: m
integer :: x0,y0,z0,x,y,i,k
integer, dimension(:) :: XX(9),YY(9)

x=x0
y=y0
i=1
k=1
XX(1:9)=(/ -1, 0, 1, -1, 0, 1, -1, 0, 1 /)
YY(1:9)=(/ -1, -1, -1, 0, 0, 0, 1, 1, 1 /) 

! check punto iniziale
do while(((x<xn-5).and.(y<yn-5).and.(x>5).and.(y>5)).and.(m(x,y,z0).eq.(0.0)))
 x=x+XX(i)
 y=y+YY(i)
 if (i>9) then
  k=k+1
  XX=XX*k
  YY=YY*k
  i=1  
 end if
end do

x0=x
y0=y

end subroutine initial_point 

! ---------------------------------------------------
! ------ INTERSEZIONE -------------------------------
! ---------------------------------------------------
subroutine intersezione(map_su,map_giu,map_uni) 
use volmap, only : xn,yn,zn,nn
implicit none
real, dimension(xn,yn,zn)  :: map_su,map_giu,map_uni
integer x,y,z

nn=0
do x=1,xn
 do y=1,yn
  do z=1,zn
   map_uni(x,y,z)=sqrt(map_su(x,y,z)*map_giu(x,y,z))
   if (map_uni(x,y,z)>0) then
     nn=nn+1
   endif
  end do
 end do
end do

end subroutine intersezione 

! ---------------------------------------------------
! ------  SEARCH CAVITY -----------------------------
! ---------------------------------------------------
subroutine search_cavity(m_in,m_out,x0,y0,z0,eps,dd)

use volmap, only : xn,yn,zn,nn

implicit none
real, intent(inout) :: m_in(xn,yn,zn)
real, intent(inout) :: m_out(xn,yn,zn)
!real, dimension(xn,yn,zn)  :: m_in,m_out
integer i,j,k,x,y,z,x0,y0,z0,dd
real eps

! crea struttura dati per costruire la lista dinamica
! node è l'elemento della lista
type node
 integer :: x,y,z
 type( node ), pointer :: next
end type node

! elementi testa e coda della lista
type (node), pointer :: head,tail,current,t

nn=0 ! nn = number of occupyied cells - for .XYZ format
do x=1,xn
 do y=1,yn
  do z=1,zn
   m_out(x,y,z)=0
  enddo
 enddo
enddo


! Inizializza lista, pulisci e crea primo nodo
nullify(head)
nullify(tail)
nullify(current)
allocate(head)
head%x=x0
head%y=y0
head%z=z0
nullify(head%next)
tail => head
current => head

do while (associated(head))
! estrai coord dal nodo di testa della lista
x=head%x
y=head%y
z=head%z
! Controlla celle adiacenti, in salita dd=1, in discesa dd=-1
 k=0
 do i=-1,1
  do j=-1,1
!  if (abs(i).ne.abs(j)) then
   ! Check coord interne
   if ((x<xn).and.(x>1).and.(y<yn).and.(y>1)) then
    ! Se valore cella != 0 salva in una nuova mappa map2
    if (m_in(x+i,y+j,z+k*dd)>eps) then
     ! crea nuovo nodo ed "aggiungilo" alla lista, in coda
     allocate(current%next)
     current => current%next
     current%x=x+i
     current%y=y+j
     current%z=z+k*dd
     nullify(current%next)
     ! map2(x+i,y+j,z+k)=map(x+i,y+j,z+k) ! salva il valore della cella in nuova mappa
     m_out(x+i,y+j,z+k*dd)=m_in(x+i,y+j,z+k*dd) ! salva la cella in nuova mappa
     m_in(x+i,y+j,z+k*dd)=0 ! "annerisci cella già controllata
     tail => current
     nn=nn+1 ! conta celle occupate in map_out, per .XYZ file
    endif
   endif
!  endif 
  end do
 end do
 k=1
   ! Check coord interne
   if ((z<zn).and.(z>1)) then
    ! Se valore cella != 0 salva in una nuova mappa map2
    if (m_in(x,y,z+k*dd)>0.0) then
     ! crea nuovo nodo ed "aggiungilo" alla lista, in coda
     allocate(current%next)
     current => current%next
     current%x=x
     current%y=y
     current%z=z+k*dd
     nullify(current%next)
     ! map2(x+i,y+j,z+k)=map(x+i,y+j,z+k) ! salva il valore della cella in nuova mappa
     m_out(x,y,z+k*dd)=m_in(x,y,z+k*dd) ! salva la cella in nuova mappa
     m_in(x,y,z+k*dd)=0 ! "annerisci cella già controllata
     tail => current
     nn=nn+1 ! conta celle occupate in map_out, per .XYZ file
    endif
   endif

! sposta la testa della lista all'elemento successivo, e pulisci memoria
t=>head
head=>head%next
deallocate(t)

end do ! end while (list is empty)

end subroutine search_cavity

! ---------------------------------------------------
! ------  FILTRA MAPPA -------------------------------
! ---------------------------------------------------
subroutine filtra_densita(m,eps)
use volmap, only : xn,yn,zn
implicit none
real, dimension(xn,yn,zn)  :: m
integer x,y,z
real :: eps

 do x=1,xn
  do y=1,yn
   do z=1,zn
    if (m(x,y,z)<eps) then
     m(x,y,z)=0
    endif
   end do !z
  end do !y
 end do !x

end subroutine filtra_densita

! ---------------------------------------------------

subroutine filtra_raggio(m,r0,x0,y0)
use volmap, only : xn,yn,zn
implicit none
real, dimension(xn,yn,zn)  :: m
integer x,y,z,x0,y0
real r0

 do x=1,xn
  do y=1,yn
   do z=1,50
    if ((x-x0)**2+(y-y0)**2 > (r0/2)**2) then
     m(x,y,z)=0
    endif
   end do !z
   do z=51,zn
    if ((x-x0)**2+(y-y0)**2 > r0**2) then
     m(x,y,z)=0
    endif
   end do !z
  end do !y
 end do !x

end subroutine filtra_raggio

! ---------------------------------------------------

subroutine filtra_ponti(m_in,m_out,x0,y0)

use volmap, only : xn,yn,zn,nn

implicit none
real, intent(inout) :: m_in(xn,yn,zn)
real, intent(inout) :: m_out(xn,yn,zn)
!real, dimension(xn,yn,zn)  :: m_in,m_out
integer i,j,x,y,z,x0,y0

! crea struttura dati per costruire la lista dinamica
! node è l'elemento della lista
type node
 integer :: x,y,z
 type( node ), pointer :: next
end type node

! elementi testa e coda della lista
type (node), pointer :: head,tail,current,t

nn=0 ! nn = number of occupyied cells - for .XYZ format
do x=1,xn
 do y=1,yn
  do z=1,zn
   m_out(x,y,z)=0
  enddo
 enddo
enddo

do z=3,zn-3 ! scorri su tutti i piani

x=x0
y=y0
! Controlla punto iniziale
!print *, x,y,z,m_in(x,y,z)
call initial_point(x,y,z,m_in)

! Inizializza lista, pulisci e crea primo nodo
nullify(head)
nullify(tail)
nullify(current)
allocate(head)
head%x=x
head%y=y
head%z=z
nullify(head%next)
tail => head
current => head


do while (associated(head))
! estrai coord dal nodo di testa della lista
x=head%x
y=head%y

! Controlla celle adiacenti, muovendoti solo sul piano xy
 do i=-1,1
  do j=-1,1
  if (abs(i).ne.abs(j)) then ! escludi celle adiacenti in diagonale
   ! Check coord interne
   if ((x<xn).and.(x>1).and.(y<yn).and.(y>1)) then
    ! Se valore cella != 0 salva in una nuova mappa map2
    if (m_in(x+i,y+j,z).ne.(0.0)) then
     ! crea nuovo nodo ed "aggiungilo" alla lista, in coda
     allocate(current%next)
     current => current%next
     current%x=x+i
     current%y=y+j
     current%z=z
     nullify(current%next)
     m_out(x+i,y+j,z)=m_in(x+i,y+j,z) ! salva la cella in nuova mappa
     m_in(x+i,y+j,z)=0 ! "annerisci cella già controllata
     tail => current
     nn=nn+1 ! conta celle occupate in map_out, per .XYZ file
    endif
   endif
  endif
  end do
 end do

! sposta la testa della lista all'elemento successivo, e pulisci memoria
t=>head
head=>head%next
deallocate(t)

end do ! end while (list is empty)

end do ! end do z, passa a piano superiore

end subroutine filtra_ponti

! ---------------------------------------------------
! ------  PRINT OUT ---------------------------------
! ---------------------------------------------------
subroutine print_output(sel,zh)

use volmap

implicit none

integer sel,ss,zh
integer x,y,z,i
real, dimension(:), allocatable :: R,A 

if (sel==1) then
! ---- Compute Average Density on a plane --------------
ss=0
do x=1,xn
 do y=1,yn
  ss=ss+int(map(x,y,zh))
 enddo
enddo
ss=ss/(xn*yn)
print *, "densità media",ss
! -------------------------------------------------------
else if (sel==2) then
! ---- Print z-Profiles: z, Resistance, Area, Radius -------------------
allocate(R(zn))
allocate(A(zn))
do z=1,zn
 A(z)=0.0
 R(z)=0.0
 do y=1,yn
  do x=1,xn
      A(z)=A(z)+map2(x,y,z)
  enddo
 enddo
 R(z)=1.0/(A(z))
enddo
write(*,"(A5,A15,A15,A15)") "z", "Resistance","Area","Radius"
do z=1,zn
 write(*,"(I5,F15.7,F15.7,F15.7)") z,R(z), A(z), sqrt(A(z)/3.14)
enddo
! -------------------------------------------------------
else if (sel==4) then
! call matrixtodx
! Print for Paraview
! call um_volume()
! ----------------------------------------------------------
endif
200 FORMAT (A6,I5,A5,A1,A3,A2,I4,A4,3F8.3,2F6.2)
end subroutine print_output 

! ----------------------------------------------------------
! ----------------------------------------------------------
! ----------------------------------------------------------


