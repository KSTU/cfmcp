!дома ночь 7-8 12
module generator    !модуль генератора
 integer(4) n1, n2, n3 !случаные числа
 integer(4) m1,a1,b1,a2,b2,m3,m2 !параметры генератора
 integer(4) max1, max2,max3 !максимумы генераторов
 real(8) randmass(500) !массив случайных чисел
 real(8) outrand
contains
subroutine randomn
!use generator
    implicit none

integer(4) i
real(8) sseed1,sseed2,sseed3
integer(4) seed(20)
character(20) a,b,c
real(8) n11,n21,n31
call date_and_time(a,b,c,seed)
sseed1=0.0
sseed2=0.0
sseed3=0.0
!print *,seed
do i=1,20,1 !здесь шарной бред
 sseed1=sseed1+seed(i)
 if (mod(i,2)==0) then
  sseed2=sseed2+seed(i)
 else
  sseed2=sseed2/2.0
 endif
 if (mod(i,3)==0) then
  sseed3=sseed3+seed(i)
 else
  sseed3=sseed3-seed(i)/2.2
 endif
enddo

!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
print *,'- Random seeds -', n11,n21,n31
!call random_seed(int(sseed))
!call random_number(n11)
!call random_number(n21)
!call random_number(n31)
m1=2147483563
a1=40014
b1=12345
m2=2147483399
a2=40692
b2=54321
max3=2147483647
max1=m1-1
max2=m2-1
n1=ceiling((n11*max1)-1)
n2=ceiling((n21*max2)-1)
n3=ceiling((n31*max3)-1)
!   open(15,file='rand.txt')
 do i=1,500
 call randstart()
 randmass(i)=outrand
 enddo
print *,'--- Generator Start ---'
end subroutine
!*******************************************************
!Рандомная функция начало
subroutine  randstart()
!use generator
    implicit none

n1=abs(mod(a1*n1+b1,m1))
n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
if (float(n3)/float(max3)<0.5) then
 outrand=float(n1)/float(max1)
else
 outrand=float(n2)/float(max2)
endif
end subroutine
!*******************************************************
!Рандомная функция
function getrand()
!use generator
    implicit none
real(8) getrand
integer(4) repick
!n1=abs(mod(a1*n1+b1,m1))
!n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
repick=ceiling((float(n3)/float(max3))*500)
getrand=randmass(repick)
call randstart()
randmass(repick)=outrand
end function

end module

module global
! global varaibles
integer(4) nv       !numbers of particle types
integer(4) nrow     !numbers of particle
integer(4) ptip     !potential type
integer(4) neq      !number of MCsteps for equilibration
integer(4) nsamp    !numbers of samlles
integer(4) nprod    !numbers of MCmove per samples
character(5),allocatable:: label(:)
real(8),allocatable:: pp1(:)
real(8),allocatable:: pp2(:)
real(8),allocatable:: pp3(:)
real(8),allocatable:: pp4(:)
real(8),allocatable:: pp5(:)
real(8),allocatable:: pp6(:)
real(8),allocatable:: p1(:,:)
real(8),allocatable:: p2(:,:)
real(8),allocatable:: p3(:,:)
real(8),allocatable:: p4(:,:)
real(8),allocatable:: p5(:,:)
real(8),allocatable:: p6(:,:)
integer(4),allocatable:: ni(:)
real(8) ro,temp       !number density, temperature
integer(4) n            !numbers of particle
real(8),allocatable:: x(:),y(:),z(:)
real(8) kubl,xgrm,xgrb,ygrm,ygrb,zgrm,zgrb
real(8) rcut
integer(4),allocatable:: tip(:)
real(8) totalen,totalden,totalten

contains    !potential functions
function pfunc(i,j,r)
    implicit none
real(8) rz
real(8) r,pfunc
integer(4) i,j

if (ptip==1) then
    rz=p1(tip(i),tip(j))/r
    rz=rz*rz
    rz=rz*rz*rz
    pfunc=4.0d0*p2(tip(i),tip(j))*(rz*rz-rz)
endif
end function

end module

program cfmcp
use generator
use global
    implicit none
print *, '---Program start---'
call randomn
call initfile()
call xyzst()
call mixrule()
call potenctest()
call totalenergy()

print *, '--- Sucsesfulli ends ---'
end program

subroutine initfile()
use global
    implicit none
integer(4) i,j,k,nom

character(80) caption
print *, 'from input file'
open(21,file='input')
    read(21,'(a)') caption
    read(21,'(i6)') nv
    print *, trim(adjustl(caption)), ': ', nv

    read(21,'(a)') caption
    read(21,'(i6)') nrow
    print *, trim(adjustl(caption)), ': ', nrow

    read(21,'(a)') caption
    read(21,'(i6)') ptip
    print *, trim(adjustl(caption)), ': ', ptip

    read(21,'(a)') caption
    read(21,'(i6)') neq
    print *, trim(adjustl(caption)), ': ', neq

    read(21,'(a)') caption
    read(21,'(i6)') nsamp
    print *, trim(adjustl(caption)), ': ', nsamp

    read(21,'(a)') caption
    read(21,'(i6)') nprod
    print *, trim(adjustl(caption)), ': ', nprod

    read(21,'(a)') caption
    read(21,'(f20.10)') ro
    print *, trim(adjustl(caption)), ': ', ro

    read(21,'(a)') caption
    read(21,'(f20.10)') temp
    print *, trim(adjustl(caption)), ': ', temp

    print *, 'Particle parameters'
    allocate(label(nv))
    allocate(pp1(nv))
    allocate(pp2(nv))
    allocate(pp3(nv))
    allocate(pp4(nv))
    allocate(pp5(nv))
    allocate(pp6(nv))
    allocate(ni(nv))
    do i=1,nv
        read(21,'(a5)') label(i)
        read(21,'(i6)') ni(i)
        read(21,'(f20.10)') pp1(i)
        read(21,'(f20.10)') pp2(i)
        read(21,'(f20.10)') pp3(i)
        read(21,'(f20.10)') pp4(i)
        read(21,'(f20.10)') pp5(i)
        read(21,'(f20.10)') pp6(i)
    enddo
close(21)
if (ptip==1) then
    print *, 'Potential Lennard-Jones'
    do i=1,nv
        print *, label(i), ' sigma: ', pp1(i), ' epsilon: ', pp2(i)
    enddo
endif

n=nrow*nrow*nrow
allocate(x(n))
allocate(y(n))
allocate(z(n))
allocate(tip(n))
kubl=(float(n)/ro)**(1.0/3.0)
rcut=kubl/2.0
print *, 'Cut radius', rcut
print *, 'Cube lenght', kubl
xgrb=kubl
xgrm=0.0
ygrb=kubl
ygrm=0.0
zgrb=kubl
zgrm=0.0
nom=0
do i=1,nrow
    do j=1,nrow
        do k=1,nrow
           nom=nom+1
           x(nom)=float(i)/float(nrow)*kubl
           y(nom)=float(j)/float(nrow)*kubl
           z(nom)=float(k)/float(nrow)*kubl
        enddo
    enddo
enddo

nom=0
do i=1,nv
    nom=nom+ni(i)
enddo
if (nom/=n) then
    print *, 'Partial particles number not equal Total'
    stop
else
    print *, 'Paticles numbers correct'
endif
nom=0
do i=1,nv
    do j=1,ni(i)
        nom=nom+1
        tip(nom)=i
    enddo
enddo

end subroutine

subroutine xyzst()
use global
    implicit none
integer(4) i

open(12,file='initial.xyz')
    write(12,'(i5)') N
    write(12,'(a)') 'comentline'
    do i=1,N
        write(12, '(a3,f10.6,f10.6,f10.6)') label(tip(i)),x(i),y(i),z(i)
    enddo
close(12)
end subroutine

subroutine mixrule()
use global
    implicit none
integer(4) i,j
allocate(p1(nv,nv))
allocate(p2(nv,nv))
allocate(p3(nv,nv))
allocate(p4(nv,nv))
allocate(p5(nv,nv))
allocate(p6(nv,nv))
print *, 'Lorenc-Bertlo'
do i=1,nv
    do j=1,nv
        p1(i,j)=(pp1(i)+pp1(j))/2.0
        p2(i,j)=dsqrt(pp2(i)*pp2(j))
        print *,label(i),label(j), 'sigma: ', p1(i,j), 'eps: ',p2(i,j)
    enddo
enddo
end subroutine

subroutine potenctest()
use global
    implicit none
real(8) r
integer(4) i

open(22,file='pottest')
do i=1,500
    r=float(i)/500.0*10.0
    write(22,'(2f20.10)') r, pfunc(1,1,r)
enddo

print *, ' --- Potenc test DONE --- '
close(22)
end subroutine

subroutine totalenergy()
use global
    implicit none
integer(4) i,j
real(8) rmin

totalen=0.0
totalden=0.0
totalten=0.0
do i=1,N
    do j=1,N
        if (j/=i) then
            call minobr(x(i),x(j),y(i),y(j),rmin)
            totanen=totalen+pfunc(i,j,rmin)
            totalden=totalden+dpfunc(i,j,rmin)
            call checknear(i)
            if (calctri==1) then
                do k=1,N
                totalten=totalten+0.0
                enddo
            endif
        endif
    enddo
enddo

end subroutine
