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
!print *, sseed1,sseed2,sseed3
!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
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
!    implicit none

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

integer(4) nv   !number of spices
integer(4) ptip !potential type
integer(4) N    !number of particles
integer(4) stor !
real(8) kubl
real(8) rcut
character(20) temps

integer(4),allocatable:: proc(:)
integer(4),allocatable:: nkol(:)
character(10),allocatable:: label(:)
real(8),allocatable:: pp1(:)
real(8),allocatable:: pp2(:)
real(8),allocatable:: pp3(:)
real(8),allocatable:: pp4(:)
real(8),allocatable:: pp5(:)

real(8),allocatable:: p1(:,:)
real(8),allocatable:: p2(:,:)
real(8),allocatable:: p3(:,:)
real(8),allocatable:: p4(:,:)
real(8),allocatable:: p5(:,:)

real(8) temp
real(8) ro
real(8) pdo,ppos

integer(8),allocatable:: atip(:)
real(8),allocatable:: x(:)
real(8),allocatable:: y(:)
real(8),allocatable:: z(:)

integer(4) moven,move,emoven,emove

real(8) xlm,xlb,ylm,ylb,zlm,zlb
real(8) xbox,ybox,zbox
real(8) xboxh,yboxh,zboxh

real(8),allocatable:: drast(:,:)
real(8),allocatable:: temrast(:)
integer(4) trcalc
real(8) totalen,totalvir,totaltri

integer(4) naccept, accept
real(8),allocatable:: nearm(:,:)
real(8),allocatable:: nearold1(:)
real(8),allocatable:: nearold2(:)
real(8),allocatable:: kchm(:)
integer(4) kchnum,rdfkol1,rdfkol2,rdfkol3,rdfnum
real(8),allocatable:: rdfm1(:,:),rdfid1(:)
real(8),allocatable:: rdfm2(:,:),rdfid2(:)
real(8),allocatable:: rdfm3(:,:),rdfid3(:)
real(8) delrdf1,delrdf2,delrdf3
real(8),allocatable:: sumrdfm1(:,:),sumrdfm2(:,:),sumrdfm3(:,:)
real(8),allocatable:: nearyes(:)
real(8) maxdl
real(8) PI
integer(4) randmove
real(8) xold,yold,zold
!================================================================

contains
function potencfunc(r,i,j)
    implicit none
real(8) r,potencfunc
real(8) u1,u2,ul,rz,rz2,rz6
integer(4) i,j
    if (ptip==1) then
        rz=p1(atip(i),atip(j))/r    !первая степень
        rz=rz*rz        !вторая степень
        rz=rz*rz*rz     !шестая степень
        u1=4.0*p2(atip(i),atip(j))*(rz*rz-rz)
        u2=0.5*p5(atip(i),atip(j))*(r-p3(atip(i),atip(j)))*(r-p3(atip(i),atip(j)))-p4(atip(i),atip(j))
        ul=0.5*(1.0+dtanh((r-0.75*p1(atip(i),atip(j)))/0.015/p1(atip(i),atip(j))))
        potencfunc=(1.0-ul)*u2+ul*u1
    endif
    if (ptip==2) then
        if (r>0.1*p1(atip(i),atip(j))) then
            rz=p1(atip(i),atip(j))/r
            rz2=rz*rz
            rz6=rz2*rz2*rz2
            potencfunc=4.0*p2(atip(i),atip(j))*(rz6*rz6-rz6)
        else
            potencfunc=99999.999
        endif
    endif
    if (ptip==3) then
        potencfunc=1.0
    endif
end function

!contains
function ptrifunc(r1,r2,r3,i,j,k)
    implicit none
real(8) ptrifunc,r1,r2,r3
integer(4) i,j,k
integer(4) pp
    ptrifunc=0.0
    pp=0
    if (r1<0.75*p1(atip(i),atip(j))) then
        pp=pp+1
    endif
    if (r2<0.75*p1(atip(j),atip(k))) then
        pp=pp+1
    endif
    if (r3<0.75*p1(atip(k),atip(i))) then
        pp=pp+1
    endif
    if (pp>1) then
        ptrifunc=99999999.0
        !print *, 'ololo'
    endif
    !if (ptip==1) then
    !ptrifunc=0.7872/r1/r1/r1/r2/r2/r2/r3/r3/r3
    !endif

end function

function dpotenc(r,i,j)
    implicit none
integer(4) i,j
real(8) r, dpotenc
real(8) rp,wp,lr,rz,rz2,rz6,sep,U1,U2

    if (ptip==1) then
        !copy start
        !rz=sqrt(RastSQR)
        !print *, p1(atip(i),atip(j)),p2(atip(i),atip(j)),p3(atip(i),atip(j)),&
        !&p4(atip(i),atip(j)),p5(atip(i),atip(j))
        !pause
        rp=0.75*p1(atip(i),atip(j)) !r_pot(Ntip(N1),Ntip(N2),i1,i2)
        !print *,rp
        wp=0.015*p1(atip(i),atip(j)) !w_pot(Ntip(N1),Ntip(N2),i1,i2)
        !print *,wp
        lr=0.5*(1.0+dtanh((r-rp)/wp))
        !print *,lr
        rz=p1(atip(i),atip(j))/r
        rz2=rz*rz
        rz6=rz2*rz2*rz2
        U1=4.0*p2(atip(i),atip(j))*(rz6*rz6-rz6)
        !print *, U1
        U2=0.5*p5(atip(i),atip(j))*(r-p3(atip(i),atip(j)))*(r-p3(atip(i),atip(j)))&
        &-p4(atip(i),atip(j))
        !print *,U2
        sep=0.5/wp/((dcosh((r-rp)/wp))*(dcosh((r-rp)/wp)))
        !print *,sep
        dpotenc=p5(atip(i),atip(j))*(r-p3(atip(i),atip(j)))*(1.0-lr)
        !print *, dpotenc
        dpotenc=dpotenc+4.0*p2(atip(i),atip(j))*(6.0*rz6-12.0*rz6*rz6)/r*lr
        !print *,dpotenc
        dpotenc=dpotenc-sep*U2
        !print *, dpotenc
        dpotenc=dpotenc+sep*U1
        !print *, dpotenc
        dpotenc=dpotenc*r
        !print *, dpotenc
    endif
    if (ptip==2) then
        if (r>0.1*p1(atip(i),atip(j))) then
            rz=p1(atip(i),atip(j))/r
            rz2=rz*rz
            rz6=rz2*rz2*rz2
            dpotenc=4.0*p2(atip(i),atip(j))*(6.0*rz6-12.0*rz6*rz6)
        else
            dpotenc=99999999999.0
        endif
    endif
    if (ptip==3) then
        dpotenc=0.0
    endif
end function
end module


program main
use generator
use global
    implicit none
    integer(4) i,j
!начало программы----------------------------------------
call randomn
call initialfiles()
call initlatice()
call initmix()
call pottest()
!call planartest()
call xyzsatic('somefile.xyz')
call totalcalc()
!do i=1,N
!    call checknear(i)
!enddo
print *,totalen,totalvir,totaltri
!pause
print *,'--- Equilibration start ---'
open(33,file='movie.xyz')
open(34,file='numbers.eq')
open(35,file='energy.eq')
!print *, getrand(), getrand(), getrand(), getrand()

print *,'--- xyz file is opened ---'
do emove=1,emoven
    !
    call mcmove()
    if (mod(emove,20)==0) then
        call xyzanim()
        call equlibout(emove)
        if (mod(emove,20000)==0) then
            print *, emove
        endif
    endif

    !print *,'stepok',emove
enddo
print *, 'Productation started'
do move=1,moven
    !
    call mcmove()
    if (mod(move,1000)==0) then
        if (mod(move,50000)==0) then
            print *,move
        endif
        call calcrdf()
        call resout()
    endif

enddo
close(34)
close(33)
close(35)
print *, ' --- Sucsesfully end --- '
!--------------------------------------------------------
end program main

subroutine initialfiles()
use global
    implicit none
integer(4) i,j
real(8) rVerh,rNiz


open(11,file='input.txt')
    read(11,'(a)') temps
    read(11,'(i5)') nv
    read(11,'(a)') temps
    read(11,'(i5)') ptip
    read(11,'(a)') temps
    read(11,'(i5)') stor
    read(11,'(a)') temps
    read(11,'(f10.5)') temp
    read(11,'(a)') temps
    read(11,'(f10.5)') ro
    read(11,'(a)') temps
    read(11,'(i10)') emoven
    read(11,'(a)') temps
    read(11,'(i10)') moven
    read(11,'(a)') temps
    read(11,'(i10)') randmove
    read(11,'(a)') temps
    read(11,'(i10)') trcalc
    allocate(label(nv+1))
    allocate(pp1(nv+1))
    allocate(pp2(nv+1))
    allocate(pp3(nv+1))
    allocate(pp4(nv+1))
    allocate(pp5(nv+1))
    allocate(p1(nv+1,nv+1))
    allocate(p2(nv+1,nv+1))
    allocate(p3(nv+1,nv+1))
    allocate(p4(nv+1,nv+1))
    allocate(p5(nv+1,nv+1))
    allocate(proc(nv))
    allocate(nkol(nv))
    do i=1,nv
        read(11,'(a)') label(i)
        read(11,'(i5)') proc(i)
        read(11,'(f10.5)') pp1(i)
        read(11,'(f10.5)') pp2(i)
        read(11,'(f10.5)') pp3(i)
        read(11,'(f10.5)') pp4(i)
        read(11,'(f10.5)') pp5(i)
    enddo
close(11)
!some
PI=3.14159265358979323846
N=stor*stor*stor
kubl=(float(N)/ro)**(1.0/3.0)
rcut=kubl/2.0
allocate(atip(N))
allocate(x(N))
allocate(y(N))
allocate(z(N))
allocate(nearm(N,N))
allocate(drast(N,N))
allocate(nearyes(N))
allocate(kchm(100))
rdfkol1=ceiling(kubl/2.0/0.05)
rdfkol2=ceiling(kubl/2.0/0.02)
rdfkol3=ceiling(kubl/2.0/0.005)
delrdf1=rcut/float(rdfkol1)
delrdf2=rcut/float(rdfkol2)
delrdf3=rcut/float(rdfkol3)
allocate(rdfm1(nv*nv,rdfkol1))
allocate(rdfm2(nv*nv,rdfkol2))
allocate(rdfm3(nv*nv,rdfkol3))
allocate(rdfid1(rdfkol1))
allocate(rdfid2(rdfkol2))
allocate(rdfid3(rdfkol3))
allocate(sumrdfm1(nv*nv,rdfkol1))
allocate(sumrdfm2(nv*nv,rdfkol2))
allocate(sumrdfm3(nv*nv,rdfkol3))
allocate(nearold1(N))
allocate(nearold2(N))

do i=1,rdfkol1  !ideal distribution
 rNiz=float((i-1))*rcut/float(rdfkol1)
 rVerh=float(i)*rcut/float(rdfkol1)
 rdfid1(i)=4.0/3.0*PI/ro*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
enddo

do i=1,rdfkol2  !ideal distribution 2
 rNiz=float((i-1))*rcut/float(rdfkol2)
 rVerh=float(i)*rcut/float(rdfkol2)
 rdfid2(i)=4.0/3.0*PI/ro*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
enddo

do i=1,rdfkol3  !ideal distribution 3
 rNiz=float((i-1))*rcut/float(rdfkol3)
 rVerh=float(i)*rcut/float(rdfkol3)
 rdfid3(i)=4.0/3.0*PI/ro*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
enddo

open(23,file='rdfid.out')  !out of ideal part
do i=1,rdfkol1
    write(23,'(f20.10,a,f20.10)') (float(i)-0.5)*rcut/float(rdfkol1), ' ', rdfid1(i)
enddo
close(23)
!stop

do i=1,100          !zero initial
    kchm(i)=0.0
enddo
do j=1,nv*nv
    do i=1,rdfkol1
        rdfm1(j,i)=0.0
    enddo
    do i=1,rdfkol2
        rdfm2(j,i)=0.0
    enddo
    do i=1,rdfkol3
        rdfm3(j,i)=0.0
    enddo
enddo
maxdl=1.0
print *, '--- Initail file loaded ---'
print *, 'Number of substances: ', nv
print *, 'Potential type: ', ptip
print *, 'Temperature: ', temp
print *, 'Number density: ', ro
print *, 'Number of particles: ', N
print *, 'Length of cube: ', kubl
print *, 'Cut radius: ', rcut
end subroutine


subroutine initlatice()
use global
    implicit none
integer(4) i,j,k,ni
real(8) dar
real(8) tempreal1,tempreal2,tempreal3,tempreal4

ni=0
dar=kubl/float(stor)
do i=1,stor
    do j=1,stor
        do k=1,stor
            ni=ni+1
            x(ni)=float(i)*dar-0.5*dar
            y(ni)=float(j)*dar-0.5*dar
            z(ni)=float(k)*dar-0.5*dar
        enddo
    enddo
enddo

ni=0
do i=1,nv
    nkol(i)=0
    do j=1,proc(i)
        ni=ni+1
        atip(ni)=i
        nkol(i)=nkol(i)+1
    enddo
    print *, 'Molecules of type: ',i,'  ', nkol(i)
enddo
if (sum(nkol)/=N) then
    print *,'Numbers of nkol no equal N'
    stop
endif
!print *,int(N*proc(1)) , nkol(1)
print *,'Initlatice DONE'

!box bounds
xlm=0.0     !minimal bounds
ylm=0.0
zlm=0.0

xlb=kubl    !maimal bounds
ylb=kubl
zlb=kubl

xbox=xlb-xlm    !box lenght
ybox=ylb-ylm
zbox=zlb-zlm

kchnum=0        !initial numbers of sampling
rdfnum=0

print *, ' X box : ',xlm, xbox,xlb
print *, ' Y box : ',ylm, ybox,ylb
print *, ' Z box : ',zlm, zbox,zlb

do i=1,N    !initial start dinstancess
    call pcalc(i,tempreal1,tempreal2,tempreal3)
enddo
print *,'cut radius ',rcut
end subroutine

subroutine xyzsatic(fname)
use global
    implicit none
character(20) fname
integer(4) i

open(21,file=fname)

    write(21,'(i5)') N
    write(21,'(a)') 'commentline'
    do i=1,N
        write(21,'(a,3f15.10)') label(atip(i)),x(i),y(i),z(i)
    enddo
close(21)
end subroutine

subroutine trans(Nmol,maxtrans)
use global
use generator
    implicit none
integer(4) Nmol
real(8) maxtrans,dx,dy,dz

dx=(getrand()-0.5)*2.0*maxtrans
dy=(getrand()-0.5)*2.0*maxtrans
dz=(getrand()-0.5)*2.0*maxtrans

x(Nmol)=x(Nmol)+dx
y(Nmol)=y(Nmol)+dy
z(Nmol)=z(Nmol)+dz

call checkbox(Nmol)

end subroutine

subroutine checkbox(Nmol)
use global
    implicit none
integer(4) Nmol

if (x(Nmol)>xlb) then
    x(Nmol)=x(Nmol)-xbox
endif
if (x(Nmol)<xlm) then
    x(Nmol)=x(Nmol)+xbox
endif
if (y(Nmol)>ylb) then
    y(Nmol)=y(Nmol)-ybox
endif
if (y(Nmol)<ylm) then
    y(Nmol)=y(Nmol)+ybox
endif

if (z(Nmol)>zlb) then
    z(Nmol)=z(Nmol)-zbox
endif
if (z(Nmol)<zlm) then
    z(Nmol)=z(Nmol)+zbox
endif

end subroutine

subroutine pcalc(Nmol,pout,triple,dpout)
use global
    implicit none
integer(4) i,j,Nmol
real(8) pout,triple,dpout
real(8) pdx,pdy,pdz
real(8) rz
integer(4) im
real(8) dpot,ddpot

pout=0.0
triple=0.0
dpout=0.0
do im=1,N
    if (im .ne. Nmol) then
    !print *, im
        call minobr(x(Nmol),x(im),y(Nmol),y(im),z(Nmol),z(im),rz)
        !print *, rz,rcut
        drast(Nmol,im)=rz
        drast(im,Nmol)=rz
        if (rz<rcut) then
            pout=pout+potencfunc(rz,Nmol,im)
            dpout=dpout+dpotenc(rz,Nmol,im)  !считать в другом месте?
        endif
    endif
enddo

triple=0.0
if (trcalc==1) then
    do i=1,N
        if (i/=Nmol) then
            do j=i+1,N
                if (j/=Nmol) then
                !tripx1=x(i)-x(Nmol)
                !tripx2=x(j)-x(Nmol)
                !tripx3=x(i)-x(j)

                !tripy1=y(i)-y(Nmol)
                !tripy2=y(j)-y(Nmol)
                !tripy3=y(i)-y(j)

                !tripz1=z(i)-z(Nmol)
                !tripz2=z(j)-z(Nmol)
                !tripz3=z(i)-z(j)

                    triple=triple+ptrifunc(drast(Nmol,i),drast(i,j),drast(j,Nmol),Nmol,i,j)
                endif
            enddo
        endif
    enddo
endif

end subroutine

subroutine totalcalc()
use global
    implicit none
real(8) tempen,tempvir,temptri
integer(4) i

totalen=0.0
totalvir=0.0
totaltri=0.0


do i=1,N
    temptri=0.0
    tempen=0.0
    tempvir=0.0
    call pcalc(i,tempen,temptri,tempvir)
    totalen=totalen+tempen
    totalvir=totalvir+tempvir
    totaltri=totaltri+temptri
enddo
totalen=totalen/2.0
totalvir=totalvir/2.0
totaltri=totaltri/3.0
end subroutine


subroutine initmix()
use global
    implicit none
!real(8)
integer(4) i,j
real(8) p5temp

print *,'init mix'
do i=1,nv
    do j=1,nv
        p1(i,j)=(pp1(i)+pp1(j))/2.0
        p2(i,j)=dsqrt(pp2(i)*pp2(j))
        p3(i,j)=(pp3(i)+pp3(j))/2.0
        p4(i,j)=dsqrt(pp4(i)*pp4(j))
        p5temp=(pp5(i)+pp5(j))/2.0
        p5(i,j)=p4(i,j)/p5temp/p5temp
        print *,'sigma',p1(i,j)
        print *,'epsi',p2(i,j)
        print *,'L',p3(i,j)
        print *,'De',p4(i,j)
        print *,'k',p5(i,j)
        print *,'DL', p5temp
    enddo
enddo
!pause
print *, 'init mix DONE'

end subroutine

subroutine pottest()
use global
    implicit none
integer(4) i,j,ir
real(8) r

open(18,file='poten.test')
do ir=1,10000
    r=float(ir)*0.001
    write(18,'(f30.15,$)') r
    do i=1,nv
        do j=1,nv
            write(18,'(3(a,f30.15),$)') ' ',potencfunc(r,i,j), ' ', dpotenc(r,i,j)&
            &, ' ', dpotenc(r,i,j)/r
        enddo
    enddo
    write(18,'(a)') ' '
enddo
close(18)
print *, 'Potential test DONE'
end subroutine

subroutine planartest()
use global
    implicit none
integer(4) i,j
real(8) maxdim
real(8) x1,y1,z1
real(8) x2,y2,z2
real(8) x3,y3,z3
real(8) r1,r2,r3
real(8) pot

open(18,file='planar.test')
maxdim=10.0
x1=p3(1,1)/2.0
x2=-p3(1,1)/2.0
print *, x1,x2
y1=0.0
y2=0.0
z1=0.0
z2=0.0
z3=0.0
do i=-150,150
    do j=-150,150
        if ((i/=0) .and. (j/=0)) then
            pot=0.0
            x3=float(i)*maxdim/150.0
            y3=float(j)*maxdim/150.0
            r1=dsqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))
            r2=dsqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))
            r3=dsqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            pot=pot+potencfunc(r1,1,1)
            pot=pot+potencfunc(r2,1,1)
            pot=pot+ptrifunc(r1,r2,r3,1,1,1)
            if (pot<10.0) then
                write (18,'(3f30.15)') x3,y3,pot
            endif
        endif
    enddo
enddo
close(18)

end subroutine

subroutine mcmove()
use global
use generator
    implicit none
integer(4) Nmol
real(8) pnew,tnew,dpnew
real(8) pold,told,dpold
integer(4) i
    Nmol=ceiling(getrand()*float(N))
    call pcalc(Nmol,pold,told,dpold)
    !pdo=pold

    !print *, Nmol
    xold=x(Nmol)
    yold=y(Nmol)
    zold=z(Nmol)
    do i=1,Nmol !old matrix of connection
        nearold1(i)=nearm(Nmol,i)
        nearold2(i)=nearm(i,Nmol)
    enddo
    !print *,pold,told,dpold
    call trans(Nmol,pp1(atip(Nmol)))

    call pcalc(Nmol,pnew,tnew,dpnew)
    !ppos=pnew
    !print *, pdo,ppos
    call checkmove(Nmol,pold,told,dpold,pnew,tnew,dpnew)
    call checknear(Nmol)


end subroutine

subroutine checkmove(Nmol,pds,pts,dds,pdn,ptn,ddn)
use global
use generator
    implicit none
integer(4) Nmol,i
real(8) pdn,ptn,ddn,pds,pts,dds
real(8) delp

delp=pdn+ptn-pds-pts
if (dexp(-(delp)/Temp)<getrand()) then
    !не принимаются
    x(Nmol)=xold
    y(Nmol)=yold
    z(Nmol)=zold
    do i=1,Nmol
        nearm(Nmol,i)=nearold1(i)
        nearm(i,Nmol)=nearold2(i)
    enddo
    naccept=naccept+1
else
    accept=accept+1
    totalen=totalen+pdn-pds
    !if (pdn-pds<-1000) then
    !    print *,totalen, pdn,pds,pdo,ppos
    !    pause
    !endif
    totalvir=totalvir+ddn-dds
    totaltri=totaltri+ptn-pts
    !print *,'proshel', totalen,totalvir,totaltri
    !pause
    !принимается
endif


end subroutine


subroutine checknear(Nmol)
use global
    implicit none
integer(4) Nmol,i
real(8) rz

do i=1,N
    nearm(Nmol,i)=0.0
    nearm(i,Nmol)=0.0
    nearyes(Nmol)=0
    if (i/=Nmol) then
        call minobr(x(Nmol),x(i),y(Nmol),y(i),z(Nmol),z(i),rz)
        if (rz<p1(atip(Nmol),atip(i))*0.75) then
            nearm(Nmol,i)=1.0
            nearm(i,Nmol)=1.0
            nearyes(Nmol)=1
            !print *,'nashel rastoyanie', chrast, p1(Nmol,i)*0.75, i ,Nmol
        endif
    endif
enddo

end subroutine

subroutine modtrans(Nmol)
use global
use generator
    implicit none
integer(4) Nmol,i
integer(4) nearmol
real(8) randalfa,randbeta,randgamma
real(8) sinox,cosox,sinoy,cosoy,sinoz,cosoz
real(8) xv,yv,zv,rast_vr

do i=1,N
    if (nearm(Nmol,i)>0.5) then
        nearmol=i
    endif
enddo
randalfa=(getrand()-0.5)*0.3
randbeta=(getrand()-0.5)*0.3
randgamma=(getrand()-0.5)*0.3
!перемещаем молекулы так чтобы вокруг которой была в центре
x(Nmol)=x(Nmol)-x(nearmol)
y(Nmol)=y(Nmol)-y(nearmol)
z(Nmol)=z(Nmol)-z(nearmol)
!print *,x(Nmol)
!print *,y(Nmol)
!print *,dsqrt(x(Nmol)*x(Nmol)+y(Nmol)*y(Nmol)+z(Nmol)*z(Nmol))

!находим растояние вокруг которого надо вращать
rast_vr=1.0+(getrand()-0.5)*0.005
!добавляем случайное изменение
x(Nmol)=x(Nmol)*rast_vr
y(Nmol)=y(Nmol)*rast_vr
z(Nmol)=z(Nmol)*rast_vr
!задаем углы поворотов
sinox=dsin(randalfa)  !0.4*(randalfa-0.5)*2.0
cosox=dcos(randalfa)  !dsqrt(1.0-sinox*sinox)
sinoy=dsin(randbeta)  !0.4*(randbeta-0.5)*2.0
cosoy=dcos(randbeta)  !dsqrt(1.0-sinoy*sinoy)
sinoz=dsin(randgamma)  !0.4*(randgamma-0.5)*2.0
cosoz=dcos(randgamma)  !dsqrt(1.0-sinoz*sinoz)
xv=x(Nmol)
yv=y(Nmol)
zv=z(Nmol)
!поворачиваем
x(Nmol)=cosox*cosoy*xv+(cosox*sinoy*sinoz-sinox*cosoz)*yv+(cosox*sinoy*cosoz+sinox*sinoz)*zv
y(Nmol)=sinox*cosoy*xv+(sinox*sinoy*sinoz+cosox*cosoz)*yv+(sinox*sinoy*cosoz-cosox*sinoz)*zv
z(Nmol)=-sinoy*xv+cosoy*sinoz*yv+cosoy*cosoz*zv
!вставляем моекулу назад
x(Nmol)=x(Nmol)+x(nearmol)
y(Nmol)=y(Nmol)+y(nearmol)
z(Nmol)=z(Nmol)+z(nearmol)

call checkbox(Nmol)

end subroutine

subroutine kch()
use global
    implicit none
integer(4) i,j
real(8) sumkch
integer(4) histkch
kchnum=kchnum+1
do i=1,N
    histkch=0
    do j=1,N
        if (i/=j) then
            if (nearm(i,j)>0.5) then
                sumkch=sumkch+1.0
                histkch=histkch+1
            endif
        endif
    enddo
    kchm(histkch)=kchm(histkch)+1.0
enddo

end subroutine

subroutine calcrdf()
use global
    implicit none
integer(4) i,j
real(8) pdx,pdy,pdz
integer(4) hist1,hist2,hist3
real(8) rz

rdfnum=rdfnum+1
do i=1,N
    do j=i,N
        if (i/=j) then
            call minobr(x(i),x(j),y(i),y(j),z(i),z(j),rz)
            if (rz<rcut) then
                !print *, rcut, rz, rz/rcut*float(rdfkol1)
                hist1=ceiling(rz/delrdf1)
                hist2=ceiling(rz/delrdf2)
                hist3=ceiling(rz/delrdf3)
                !print *,hist1,rdfkol1
                !pause
                rdfm1((atip(i)-1)*nv+atip(j),hist1)=rdfm1((atip(i)-1)*nv+atip(j),hist1)+1.0
                rdfm2((atip(i)-1)*nv+atip(j),hist2)=rdfm2((atip(i)-1)*nv+atip(j),hist2)+1.0
                rdfm3((atip(i)-1)*nv+atip(j),hist3)=rdfm3((atip(i)-1)*nv+atip(j),hist3)+1.0
            endif
        endif
    enddo
enddo

end subroutine


subroutine resout()
use global
    implicit none
integer(4) i,j
real(8) obshkch,obshkch2
!kch output
open(13,file='rdf1.out')
do i=1,rdfkol1
    write(13,'(f20.10,a,$)') (float(i)-0.5)*rcut/float(rdfkol1), ' '
    do j=1,nv*nv
        write(13,'(f20.10,a,$)') rdfm1(j,i)/float(rdfnum)/rdfid1(i), ' '
        write(13,'(f20.10,a,$)') rdfm1(j,i)/float(rdfnum), ' '
    enddo
    write(13,'(a)') ' '
enddo
close(13)
open(13,file='rdf2.out')
do i=1,rdfkol2
    write(13,'(f20.10,a,$)') (float(i)-0.5)*rcut/float(rdfkol2), ' '
    do j=1,nv*nv
        write(13,'(f20.10,a,$)') rdfm2(j,i)/float(rdfnum)/rdfid2(i), ' '
        write(13,'(f20.10,a,$)') rdfm2(j,i), ' '
    enddo
    write(13,'(a)') ' '
enddo
close(13)
open(13,file='rdf3.out')
do i=1,rdfkol3
    write(13,'(f20.10,a,$)') (float(i)-0.5)*rcut/float(rdfkol3), ' '
    do j=1,nv*nv
        write(13,'(f20.10,a,$)') rdfm3(j,i)/float(rdfnum)/rdfid3(i), ' '
    enddo
    write(13,'(a)') ' '
enddo
close(13)
!
open(13,file='densdist.out')
    !распределение по плотности
close(13)
open(13,file='kch.out')
    !распределение по координационному числу
    obshkch=0.0
    obshkch2=0.0
    do i=1,100
        obshkch=obshkch+kchm(i)
        obshkch2=obshkch2+kchm(i)/float(i)
    enddo
    do i=1,100
        write (13,'(f20.10,a,f20.10,a,f20.10)') float(i),' ', kchm(i)/obshkch, ' ',&
        & kchm(i)/float(i)/obshkch2
    enddo
close(13)
end subroutine

subroutine xyzanim()
use global
    implicit none
integer(4) i

write(33,'(i5)') N
write(33,'(a)') 'comentline'
do i=1,N
    write(33,'(a,3f15.10)') label(atip(i)),x(i),y(i),z(i)
enddo
end subroutine

subroutine changemaxdl()
use global
    implicit none

if (float(accept)/float(naccept)>0.6) then
    maxdl=maxdl*1.05
endif

if (float(accept)/float(naccept)<0.4) then
    maxdl=maxdl*0.95
endif
accept=0
naccept=0

end subroutine

subroutine equlibout(nmove)
use global
    implicit none
integer(4) histkch
integer(4) kch0,kch1,kch2,kch3,kch4,kch5
real(8) kchsum1,kchsum2
integer(4) kchob
integer(4) nmove,i,j
!текущее КЧ
kchob=0
kch0=0
kch1=0
kch2=0
kch3=0
kch4=0
kch5=0
kchsum1=0.0
kchsum2=0.0

do i=1,N
    histkch=0
    do j=1,N
        if (i/=j) then
            if (nearm(i,j)>0.5) then
                kchob=kchob+1
                histkch=histkch+1
                !print *, nearm(i,j),histkch
            endif
        endif
    enddo
    !pause
    !print *,i,histkch
    if (histkch==0) then
        kch0=kch0+1
    endif
    if (histkch==1) then
        kch1=kch1+1
    endif
    if (histkch==2) then
        kch2=kch2+1
    endif
    if (histkch==3) then
        kch3=kch3+1
    endif
    if (histkch==4) then
        kch4=kch4+1
    endif
    if (histkch==5) then
        kch5=kch5+1
    endif
enddo
kchsum1=float(kch0+kch1+kch2+kch3+kch4+kch5)
kchsum2=float(kch0)+float(kch1)/2.0+float(kch2)/3.0+float(kch3)/4.0+float(kch4)/5.0&
        &+float(kch5)/6.0
!print *, kchob, kch0,kch1,kch2
!print *,nearm
!pause

write(34,'(i8,13f20.15)') nmove,float(kch0)/kchsum1,float(kch1)/kchsum1,float(kch2)/kchsum1,&
    &float(kch3)/kchsum1,float(kch4)/kchsum1,float(kch5)/kchsum1,float(kch0)/kchsum2,&
    &float(kch1)/kchsum2/2.0&
    &,float(kch2)/3.0/kchsum2,&
    float(kch3)/4.0/kchsum2,float(kch4)/5.0/kchsum2,float(kch5)/6.0/kchsum2,float(kchob)/float(N*N)
write(35,'(i8,3f40.15)') nmove,totalen,totaltri,totalvir
end subroutine

subroutine minobr(x1i,x2i,y1i,y2i,z1i,z2i,r)
use global
    implicit none
integer(4) i
real(8) pdx,pdy,pdz
real(8) x1i,x2i,y1i,y2i,z1i,z2i
real(8) r,xv,yv,zv

if (x1i-x2i>xboxh) then
    xv=x2i+xbox
else if (x2i-x1i>xboxh) then
    xv=x2i-xbox
else
    xv=x2i
endif
pdx=xv-x1i


if (y1i-y2i>yboxh) then
    yv=y2i+ybox
else if (y2i-y1i>yboxh) then
    yv=y2i-ybox
else
    yv=y2i
endif
pdy=yv-y1i

if (z1i-z2i>zboxh) then
    zv=z2i+zbox
else if (z2i-z1i>zboxh) then
    zv=z2i-zbox
else
    zv=z2i
endif
pdz=zv-z1i

r=dsqrt(pdx*pdx+pdy*pdy+pdz*pdz)

end subroutine
