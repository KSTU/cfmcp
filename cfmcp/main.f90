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

real(8),allocatable:: proc(:)
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

integer(8),allocatable:: atip(:)
real(8),allocatable:: x(:)
real(8),allocatable:: y(:)
real(8),allocatable:: z(:)

integer(4) moven,move,emoven,emove

real(8) xlm,xlb,ylm,ylb,zlm,zlb
real(8) xbox,ybox,zbox

real(8),allocatable:: drast(:,:)
real(8),allocatable:: temrast(:)
integer(4) trcalc
real(8) totalen,totalvir,totaltri

integer(4) naccept, accept
real(4),allocatable:: nearm(:,:)
real(4),allocatable:: kchm(:)
integer(4) kchnum,rdfkol,rdfnum
real(8),allocatable:: rdfm(:),rdfid(:)
real(8),allocatable:: nearyes(:)
real(8) maxdl
real(8) PI
integer(4) randmove
!================================================================

contains
function potencfunc(r,i,j)
real(8) r,potencfunc
real(8) u1,u2,ul,rz
integer(4) i,j
    if (ptip==1) then
        rz=p1(i,j)/r    !первая степень
        rz=rz*rz        !вторая степень
        rz=rz*rz*rz     !шестая степень
        u1=4.0*p2(i,j)*(rz*rz-rz)
        u2=p5(i,j)*(r-p3(i,j))*(r-p3(i,j))-p4(i,j)
        ul=0.5*(1.0+dtanh((r-0.75*p1(i,j))/0.015/p1(i,j)))
        potencfunc=(1.0-ul)*u2+ul*u1
    endif
end function

!contains
function ptrifunc(r1,r2,r3,i,j,k)
real(8) ptrifunc,r1,r2,r3
integer(4) i,j,k
integer(4) pp
    ptrifunc=0.0
    pp=0
    if (r1<0.75*p1(i,j)) then
        pp=pp+1
    endif
    if (r2<0.75*p1(j,k)) then
        pp=pp+1
    endif
    if (r3<0.75*p1(k,i)) then
        pp=pp+1
    endif
    if (pp>1) then
        ptrifunc=999999999.0
    endif
    !if (ptip==1) then
    !ptrifunc=0.7872/r1/r1/r1/r2/r2/r2/r3/r3/r3
    !endif

end function

function dpotenc(r,i,j)
real(8) r, dpotenc
    if (ptip==1) then
        dpotenc=0.0
    endif
end function
end module


program main
use generator
use global
    implicit none
    integer(4) i
!начало программы----------------------------------------
call randomn
print *,'--- Generator Start ---'
call initialfiles()
call initlatice()
call initmix()
call pottest()
call planartest()
call xyzsatic('somefile.xyz')
trcalc=0
call totalcalc()
do i=1,N
    call checknear(i)
enddo
print *,totalen,totalvir,totaltri
!pause
print *,'--- Equilibration start ---'
open(33,file='movie.xyz')
open(34,file='equlibration.out')
print *,'--- xyz file is opened ---'
do emove=1,emoven
    !
    call mcmove()
    call xyzanim()
    call equlibout(emove)
    !print *,'stepok',emove
enddo
do move=1,moven
    !

enddo
close(34)
close(33)
print *, ' --- Sucsesfully end --- '
!--------------------------------------------------------
end program main

subroutine initialfiles()
use global
    implicit none
integer(4) i
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
        read(11,'(f10.5)') proc(i)
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
allocate(drast(N,N))
allocate(atip(N))
allocate(x(N))
allocate(y(N))
allocate(z(N))
allocate(nearm(N,N))
allocate(nearyes(N))
allocate(kchm(100))
rdfkol=ceiling(kubl/2.0/0.1)
allocate(rdfm(rdfkol))
allocate(rdfid(rdfkol))
do i=1,rdfkol
 rNiz=float((i-1))*rcut/float(rdfkol)
 rVerh=float(i)*rcut/float(rdfkol)
 rdfid(i)=4.0/3.0*PI/ro*(rVerh*rVerh*rVerh-rNiz*rNiz*rNiz)
enddo
do i=1,100
    kchm(i)=0.0
enddo
do i=1,rdfkol
    rdfm(i)=0.0
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

ni=0
dar=kubl/float(stor)
do i=1,stor
    do j=1,stor
        do k=1,stor
            ni=ni+1
            x(ni)=float(i)*dar
            y(ni)=float(j)*dar
            z(ni)=float(k)*dar
        enddo
    enddo
enddo

ni=0
do i=1,nv
    nkol(i)=0
    do j=1,int(N*proc(i))
    ni=ni+1
    atip(ni)=i
    nkol(i)=nkol(i)+1
    enddo
    print *, 'Molecules of type: ',i,'  ', nkol(i)
enddo
!print *,int(N*proc(1)) , nkol(1)
print *,'Initlatice DONE'

xlm=0.0
ylm=0.0
zlm=0.0

xlb=kubl
ylb=kubl
zlb=kubl

xbox=xlb-xlm
ybox=ylb-ylm
zbox=zlb-zlm

kchnum=0
rdfnum=0

print *, ' X box : ', xbox
print *, ' Y box : ', ybox
print *, ' Z box : ', zbox

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
real(8) maxtrans

x(Nmol)=x(Nmol)+(getrand()-0.5)*2.0*maxtrans
y(Nmol)=y(Nmol)+(getrand()-0.5)*2.0*maxtrans
z(Nmol)=z(Nmol)+(getrand()-0.5)*2.0*maxtrans

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

pout=0.0
do i=1,N
    if (i/=Nmol) then
        pdx=abs(x(i)-x(Nmol))
        if (pdx>xbox/2.0) then
            pdx=xbox-pdx
        endif
        pdy=abs(y(i)-y(Nmol))
        if (pdy>ybox/2.0) then
            pdy=ybox-pdy
        endif
        pdz=abs(z(i)-z(Nmol))
        if (pdz>ybox/2.0) then
            pdz=ybox-pdz
        endif
        rz=sqrt(pdx*pdx+pdy*pdy+pdz*pdz)
        drast(Nmol,i)=rz
        if (rz<rcut) then
            pout=pout+potencfunc(rz,Nmol,i)
            dpout=dpout+dpotenc(rz,Nmol,i)  !считать в другом месте?
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

temptri=0.0
tempen=0.0
tempvir=0.0
do i=1,N
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

do i=1,nv
    do j=1,nv
        p1(i,j)=(pp1(i)+pp1(j))/2.0
        p2(i,j)=dsqrt(pp2(i)*pp2(j))
        p3(i,j)=(pp3(i)+pp3(j))/2.0
        p4(i,j)=dsqrt(pp4(i)*pp4(j))
        p5temp=(pp5(i)+pp5(j))/2.0
        p5(i,j)=p4(i,j)/p5temp/p5temp
        print *,p5(i,j)
    enddo
enddo


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
            write(18,'(a,f30.15,$)') ' ',potencfunc(r,i,j)
        enddo
    enddo
    write(18,'(a)') ' '
enddo
close(18)

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
real(8) pold,told,dpold    !
real(8) xold,yold,zold
    Nmol=ceiling(getrand()*float(N))
    !print *, Nmol
    xold=x(Nmol)
    yold=y(Nmol)
    zold=z(Nmol)
    call pcalc(Nmol,pold,told,dpold)
    !print *,pold,told,dpold
    if (getrand()>0.5) then
        if (randmove==1) then
            if (getrand()<0.5) then
                call trans(Nmol,0.15*pp1(Nmol))
            else
                call trans(Nmol,pp1(Nmol))
            endif
        else
            call trans(Nmol,maxdl)
        endif
        !print *,'transok'
    else
        if (nearyes(Nmol)==1) then
            !call modtrans(Nmol)
            !print *, 'modtrans ok'
            call trans(Nmol,maxdl)
        else
            if (randmove==1) then
                if (getrand()<0.5) then
                    call trans(Nmol,0.15*pp1(Nmol))
                else
                    call trans(Nmol,pp1(Nmol))
                endif
            else
                call trans(Nmol,maxdl)
            endif
        endif
        !print *,'modtransok'
    endif
    call pcalc(Nmol,pnew,tnew,dpnew)
    call checkmove(Nmol,pold+told,pnew+tnew,xold,yold,zold)
    call checknear(Nmol)


end subroutine

subroutine checkmove(Nmol,pnews,polds,xold,yold,zold)
use global
use generator
    implicit none
integer(4) Nmol
real(8) pnews,polds
real(8) xold,yold,zold

if (dexp(-(pnews-polds)/Temp)<getrand()) then
    !не принимаются
    x(Nmol)=xold
    y(Nmol)=yold
    z(Nmol)=zold
    naccept=naccept+1
else
    accept=accept+1
    !принимается
endif


end subroutine


subroutine checknear(Nmol)
use global
    implicit none
integer(4) Nmol,i
real(8) chrast,chdx,chdy,chdz

do i=1,N
    nearm(Nmol,i)=0
    nearm(i,Nmol)=0
    nearyes(Nmol)=0
    if (i/=Nmol) then
        chdx=x(Nmol)-x(i)
        chdy=y(Nmol)-y(i)
        chdz=z(Nmol)-z(i)
        chrast=dsqrt(chdx*chdx+chdy*chdy+chdz*chdz)
        if (chrast<p1(Nmol,i)*0.75) then
            nearm(Nmol,i)=1
            nearm(i,Nmol)=1
            nearyes(Nmol)=1
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
    if (nearm(Nmol,i)==1) then
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
print *,dsqrt(x(Nmol)*x(Nmol)+y(Nmol)*y(Nmol)+z(Nmol)*z(Nmol))

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
            if (nearm(i,j)==1) then
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
integer(4) hist
real(8) rz

rdfnum=rdfnum+1
do i=1,N
    do j=1,N
        if (i/=j) then
        pdx=abs(x(i)-x(j))
        if (pdx>xbox/2.0) then
            pdx=xbox-pdx
        endif
        pdy=abs(y(i)-y(j))
        if (pdy>ybox/2.0) then
            pdy=ybox-pdy
        endif
        pdz=abs(z(i)-z(j))
        if (pdz>zbox) then
            pdz=zbox-pdz
        endif
        rz=sqrt(pdx*pdx+pdy*pdy+pdz*pdz)
        if (rz<rcut) then
            hist=ceiling(rz/rcut*float(rdfkol))
            rdfm(hist)=rdfm(hist)+1.0
        endif
        endif
    enddo
enddo

end subroutine


subroutine resout()
use global
    implicit none
integer(4) i
real(8) obshkch,obshkch2
!kch output
open(13,file='rdf.out')
do i=1,rdfkol
    write(13,'(f20.10,a,f20.10)') (float(i)-0.5)/rcut/float(rdfkol), ' ', rdfm(i)/rdfnum*rdfid(i)
enddo
close(13)
!
open(13,file='densdist.out')
    !распределение по плотности
close(13)
open(13,file='')
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
            if (nearm(i,j)==1) then
                kchob=kchob+1
                histkch=histkch+1
            endif
        endif
    enddo
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

write(34,'(i8,13f20.15)') nmove,float(kch0)/kchsum1,float(kch1)/kchsum1,float(kch2)/kchsum1,&
    &float(kch3)/kchsum1,float(kch4)/kchsum1,float(kch5)/kchsum1,float(kch0)/kchsum2,&
    &float(kch1)/kchsum2/2.0&
    &,float(kch2)/3.0/kchsum2,&
    float(kch3)/4.0/kchsum2,float(kch4)/5.0/kchsum2,float(kch5)/6.0/kchsum2,float(kchob)/float(N*N)

end subroutine