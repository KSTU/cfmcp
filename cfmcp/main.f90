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
integer(4) calctri
real(8) xbox,xboxh,ybox,yboxh,zbox,zboxh
integer(4),allocatable:: nearm(:,:)
real(8) maxdl
integer(4) accept,naccept
real(8),allocatable:: sumten(:)
real(8),allocatable:: sumtden(:)
real(8),allocatable:: sumtten(:)
integer(4), allocatable:: sumsamp(:)
real(8) pi,ecorr,pcorr
real(8),allocatable:: sumrdf1(:,:),sumrdf2(:,:),sumrdf3(:,:)
real(8),allocatable:: idrdf1(:,:),idrdf2(:,:),idrdf3(:,:)
real(8) drrdf1,drrdf2,drrdf3
integer(4) rdfnum
integer(4),allocatable:: connected(:)
real(8),allocatable:: skch10(:),skch11(:),skch12(:),skch13(:),skch14(:),skch15(:)
real(8),allocatable:: skch20(:),skch21(:),skch22(:),skch23(:),skch24(:),skch25(:)
integer(4),allocatable:: kchnum(:)
integer(4) calcpe
real(8),allocatable:: sumdens1(:),sumdens2(:),sumdens3(:)
integer(4) numdens
real(8) drdens1,drdens2,drdens3
integer(4) indens1,indens2,indens3
!---functions----------------------------------------------------
contains    !potential functions
function pfunc(i,j,r)
    implicit none
real(8) rz
real(8) r,pfunc
integer(4) i,j
real(8) u1,u2,ul
real(8) rp,wp

pfunc=0.0
if (ptip==1) then
    rz=p1(tip(i),tip(j))/r
    rz=rz*rz
    rz=rz*rz*rz
    pfunc=4.0d0*p2(tip(i),tip(j))*(rz*rz-rz)
endif
if (ptip==2) then
    rz=p1(tip(i),tip(j))/r    !rm/r
    rz=rz*rz                  !(rm/r)^2
    rz=rz*rz*rz               !(rm/r)^6
    rp=0.75*p1(tip(i),tip(j))
    wp=0.015*p1(tip(i),tip(j))
    u1=4.0*p2(tip(i),tip(j))*(rz*rz-rz)
    u2=0.5*p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(r-p3(tip(i),tip(j)))-p4(tip(i),tip(j))
    ul=0.5*(1.0+dtanh((r-rp)/wp))
    pfunc=(1.0-ul)*u2+ul*u1
endif
if (ptip==3) then
    if (i/=j) then
        rz=p1(tip(i),tip(j))/r    !rm/r
        rz=rz*rz                  !(rm/r)^2
        rz=rz*rz*rz               !(rm/r)^6
        rp=0.75*p1(tip(i),tip(j))
        wp=0.015*p1(tip(i),tip(j))
        u1=4.0*p2(tip(i),tip(j))*(rz*rz-rz)
        u2=0.5*p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(r-p3(tip(i),tip(j)))-p4(tip(i),tip(j))
        ul=0.5*(1.0+dtanh((r-rp)/wp))
        pfunc=(1.0-ul)*u2+ul*u1
    else
        rz=p1(tip(i),tip(j))/r
        rz=rz*rz
        rz=rz*rz*rz
        pfunc=4.0d0*p2(tip(i),tip(j))*(rz*rz-rz)
    endif
endif

end function

function dpfunc(i,j,r)
    implicit none
integer(4) i,j
real(8) r,dpfunc,rz,rz2,rz6
real(8) lr,rp,U1,U2,sep,wp

dpfunc=0.0
if (ptip==1) then
    rz=p1(tip(i),tip(j))/r
    rz=rz*rz
    rz=rz*rz*rz
    dpfunc=4.0d0*p2(tip(i),tip(j))*(6.0*rz-12.0*rz*rz)
endif
if (ptip==2) then
    rp=0.75*p1(tip(i),tip(j))
    wp=0.015*p1(tip(i),tip(j))
    lr=0.5*(1.0+dtanh((r-rp)/wp))
    rz=p1(tip(i),tip(j))/r
    rz2=rz*rz
    rz6=rz2*rz2*rz2
    U1=4.0*p2(tip(i),tip(j))*(rz6*rz6-rz6)
    U2=0.5*p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(r-p3(tip(i),tip(j)))&
    &-p4(tip(i),tip(j))
    sep=0.5/wp/((dcosh((r-rp)/wp))*(dcosh((r-rp)/wp)))
    dpfunc=p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(1.0-lr)
    dpfunc=dpfunc+4.0*p2(tip(i),tip(j))*(6.0*rz6-12.0*rz6*rz6)/r*lr
    dpfunc=dpfunc-sep*U2
    dpfunc=dpfunc+sep*U1
    dpfunc=dpfunc*r
endif
if (ptip==3) then
    if (i/=j) then
        rp=0.75*p1(tip(i),tip(j))
        wp=0.015*p1(tip(i),tip(j))
        lr=0.5*(1.0+dtanh((r-rp)/wp))
        rz=p1(tip(i),tip(j))/r
        rz2=rz*rz
        rz6=rz2*rz2*rz2
        U1=4.0*p2(tip(i),tip(j))*(rz6*rz6-rz6)
        U2=0.5*p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(r-p3(tip(i),tip(j)))&
        &-p4(tip(i),tip(j))
        sep=0.5/wp/((dcosh((r-rp)/wp))*(dcosh((r-rp)/wp)))
        dpfunc=p5(tip(i),tip(j))*(r-p3(tip(i),tip(j)))*(1.0-lr)
        dpfunc=dpfunc+4.0*p2(tip(i),tip(j))*(6.0*rz6-12.0*rz6*rz6)/r*lr
        dpfunc=dpfunc-sep*U2
        dpfunc=dpfunc+sep*U1
        dpfunc=dpfunc*r
    else
        rz=p1(tip(i),tip(j))/r
        rz=rz*rz
        rz=rz*rz*rz
        dpfunc=4.0d0*p2(tip(i),tip(j))*(6.0*rz-12.0*rz*rz)
    endif
endif
end function
function tfunc(i,j,k,r1,r2,r3)
    implicit none
real(8) tfunc
integer(4) i,j,k
real(8) r1,r2,r3
integer(4) check

tfunc=0.0
if (calctri==1) then
    check=0
    if (r1<0.75*p1(tip(i),tip(j))) then
        check=check+1
    endif
    if (r2<0.75*p1(tip(i),tip(k))) then
        check=check+1
    endif
    if (r3<0.75*p1(tip(j),tip(k))) then
        check=check+1
    endif
    if (check>1.5) then
        tfunc=99999999.0
        !print *,'ololo', r1,r2,r3
    else
        tfunc=0.0
    endif
endif
end function

end module

program cfmcp
use generator
use global
    implicit none
integer(4) move,cmol
real(8) pold,dpold,told
real(8) pnew,dpnew,tnew
real(8) xnew,ynew,znew
integer csample

print *, '---Program start---'
call randomn
call initfile()
call xyzst()
call mixrule()
call potenctest()
call totalenergy()
call tailcorr()
!-----open write files
open(31,file='equil.xyz')
print *, ' --- Start Equilibration ---'
open(56,file='equilib.out')
open(34,file='enumber.out')
do move=1,neq         !-----------equilibration
    cmol=ceiling(getrand()*n)   !get currwnt molecule
    call pcalc(cmol,pold,dpold,told,x(cmol),y(cmol),z(cmol))
    call trans(cmol,xnew,ynew,znew)
    call pcalc(cmol,pnew,dpnew,tnew,xnew,ynew,znew)
    call checkmove(cmol,pold,dpold,told,pnew,dpnew,tnew,xnew,ynew,znew)
    call checknear(cmol)
    !print *,told,tnew
    !print *, totalen, pold, pnew, float(accept)/float(accept+naccept)
    if (mod(move,40)==0) then
        call equilibout(move)
    endif
    if (mod(move,200)==0) then
        !call xyzanim()
        if (mod(move,10000)==0) then
            call checkdl()
        endif
    endif
    if (mod(move,ceiling(float(neq)/2.0))==0) then
        zgrb=kubl+kubl*float(calcpe)
        zgrm=-kubl*float(calcpe)
        zbox=zgrb-zgrm
        zboxh=zbox/2.0
        print *, 'z- from ', zgrm, ' to ', zgrb
        print *, kubl, zbox
    endif
enddo
close(56)
print *, ' --- Start Productation ---'
call densinit()
open(56,file='prod.out')
do csample=1,nsamp       !----Productation
    print *, 'Sample # ', csample, ' of ', nsamp, ' started'
    do move=1,nprod
        cmol=ceiling(getrand()*n)   !get currwnt molecule
        call pcalc(cmol,pold,dpold,told,x(cmol),y(cmol),z(cmol))
        call trans(cmol,xnew,ynew,znew)
        call pcalc(cmol,pnew,dpnew,tnew,xnew,ynew,znew)
        call checkmove(cmol,pold,dpold,told,pnew,dpnew,tnew,xnew,ynew,znew)
        call checknear(cmol)
        call dosample(csample)
        if (mod(move,40)==0) then
        !    call equilibout(move+(csample-1)*nsamp)
            call calcrdf()
            call calckch(csample)
            call calcdens()
        endif
        if (mod(move,400)==0) then
            !call xyzanim()
            if (mod(move,4000)==0) then
                call resout(csample)
            endif
        endif
    enddo
    print *, 'Potential energy: ',sumten(csample)/float(sumsamp(csample))/float(N)+ecorr !,&
    !&sumten(csample),sumsamp(csample)
    print *, 'Pressure: ', ro*Temp-1.0/3.0*ro/float(n)*sumtden(csample)/&
    &float(sumsamp(csample)) +pcorr !,sumtden(csample),sumsamp(csample)
    print *, '3-body energy: ', sumtten(csample)/float(sumsamp(csample))/float(N) !,&
   ! &sumtten(csample),sumsamp(csample)
enddo
close(56)
close(34)
print *, '--- Sucsesfulli ends ---'
!---close write files
close(31)
end program

subroutine initfile()
use global
    implicit none
integer(4) i,j,k,nom
character(80) caption
real(8) rv,rn

pi=3.141592653589793238462643383279
print *, 'from input file'
open(21,file='input')
    read(21,'(a)') caption
    read(21,'(i8)') nv
    print *, trim(adjustl(caption)), ': ', nv

    read(21,'(a)') caption
    read(21,'(i8)') nrow
    print *, trim(adjustl(caption)), ': ', nrow

    read(21,'(a)') caption
    read(21,'(i8)') ptip
    print *, trim(adjustl(caption)), ': ', ptip

    read(21,'(a)') caption
    read(21,'(i8)') neq
    print *, trim(adjustl(caption)), ': ', neq

    read(21,'(a)') caption
    read(21,'(i8)') nsamp
    print *, trim(adjustl(caption)), ': ', nsamp

    read(21,'(a)') caption
    read(21,'(i8)') nprod
    print *, trim(adjustl(caption)), ': ', nprod

    read(21,'(a)') caption
    read(21,'(f20.10)') ro
    print *, trim(adjustl(caption)), ': ', ro

    read(21,'(a)') caption
    read(21,'(f20.10)') temp
    print *, trim(adjustl(caption)), ': ', temp

    read(21,'(a)') caption
    read(21,'(i5)') calctri
    print *, trim(adjustl(caption)), ': ', calctri

    read(21,'(a)') caption
    read(21,'(i5)') calcpe
    print *, trim(adjustl(caption)), ': ', calcpe

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
if (ptip==2) then
    print *, 'Lennard-Jones with harmnic'
    do i=1,nv
        print *, label(i), ' sigma: ', pp1(i), ' epsilon: ', pp2(i), ' L: ', pp3(i), &
        & ' De: ', pp4(i), ' k: ', pp5(i)
    enddo
endif

n=nrow*nrow*nrow
allocate(x(n))
allocate(y(n))
allocate(z(n))
allocate(tip(n))
allocate(nearm(n,n))
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
xbox=xgrb-xgrm
xboxh=xbox/2.0
ybox=ygrb-ygrm
yboxh=ybox/2.0
zbox=zgrb-zgrm
zboxh=zbox/2.0
nom=0
do i=1,nrow
    do j=1,nrow
        do k=1,nrow
           nom=nom+1
           x(nom)=(float(i)-0.5)/float(nrow)*kubl
           y(nom)=(float(j)-0.5)/float(nrow)*kubl
           z(nom)=(float(k)-0.5)/float(nrow)*kubl
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
maxdl=1.0
!allocating sums for sample
allocate(sumten(nsamp))
allocate(sumtden(nsamp))
allocate(sumtten(nsamp))
allocate(sumsamp(nsamp))
allocate(skch10(nsamp))
allocate(skch11(nsamp))
allocate(skch12(nsamp))
allocate(skch13(nsamp))
allocate(skch14(nsamp))
allocate(skch15(nsamp))
allocate(skch20(nsamp))
allocate(skch21(nsamp))
allocate(skch22(nsamp))
allocate(skch23(nsamp))
allocate(skch24(nsamp))
allocate(skch25(nsamp))
allocate(kchnum(nsamp))
sumsamp=0
sumten=0.0
sumtden=0.0
sumtten=0.0
drrdf1=rcut/100.0
drrdf2=rcut/200.0
drrdf3=rcut/300.0
allocate(sumrdf1(ceiling(rcut/drrdf1),nv*nv))
allocate(sumrdf2(ceiling(rcut/drrdf2),nv*nv))
allocate(sumrdf3(ceiling(rcut/drrdf3),nv*nv))
rdfnum=0
!sumrdf1=0.0
!sumrdf2=0.0
!sumrdf3=0.0
allocate(idrdf1(ceiling(rcut/drrdf1),nv*nv))
allocate(idrdf2(ceiling(rcut/drrdf2),nv*nv))
allocate(idrdf3(ceiling(rcut/drrdf3),nv*nv))
open(55,file='idrdf1.out')
do i=1,nv
    do j=1,nv
        do k=1,ceiling(rcut/drrdf1)
            rv=drrdf1*float(k)
            rn=drrdf1*float(k-1)
            idrdf1(k,(i-1)*nv+j)=4.0/3.0*pi*ro*(rv*rv*rv-rn*rn*rn)*float(ni(i)*ni(j))/float(n*n)
            sumrdf1(k,(i-1)*nv+j)=0.0
            write(55,'(2f20.10)') (rv+rn)/2.0, idrdf1(k,(i-1)*nv+j)
        enddo
    enddo
enddo
close(55)
!open(55,file='idrdf2.out')
do i=1,nv
    do j=1,nv
        do k=1,ceiling(rcut/drrdf2)
            rv=drrdf2*float(k)
            rn=drrdf2*float(k-1)
            idrdf2(k,(i-1)*nv+j)=4.0/3.0*pi*ro*(rv*rv*rv-rn*rn*rn)*float(ni(i)*ni(j))/float(n*n)
            sumrdf2(k,(i-1)*nv+j)=0.0
!            write(55,'(2f20.10)') (rv+rn)/2.0, idrdf2(k,(i-1)*nv+j)
        enddo
    enddo
enddo
!close(55)
!open(55,file='idrdf3.out')
do i=1,nv
    do j=1,nv
        do k=1,ceiling(rcut/drrdf3)
            rv=drrdf3*float(k)
            rn=drrdf3*float(k-1)
            idrdf3(k,(i-1)*nv+j)=4.0/3.0*pi*ro*(rv*rv*rv-rn*rn*rn)*float(ni(i)*ni(j))/float(n*n)
            sumrdf3(k,(i-1)*nv+j)=0.0
!            write(55,'(2f20.10)') (rv+rn)/2.0, idrdf3(k,(i-1)*nv+j)
        enddo
    enddo
enddo
!close(55)
allocate(connected(n))
skch10=0.0
skch11=0.0
skch12=0.0
skch13=0.0
skch14=0.0
skch15=0.0
skch20=0.0
skch21=0.0
skch22=0.0
skch23=0.0
skch24=0.0
skch25=0.0
kchnum=0
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
if (ptip==1) then
    print *, 'Lorenc-Bertlo'
    do i=1,nv
        do j=1,nv
            p1(i,j)=(pp1(i)+pp1(j))/2.0
            p2(i,j)=dsqrt(pp2(i)*pp2(j))
            print *,label(i),label(j), 'sigma: ', p1(i,j), 'eps: ',p2(i,j)
        enddo
    enddo
endif
if ((ptip==2) .or. (ptip==3)) then
    print *, 'Lorenc-Berthlo+GE'
    do i=1,nv
        do j=1,nv
            p1(i,j)=(pp1(i)+pp1(j))/2.0
            p2(i,j)=dsqrt(pp2(i)*pp2(j))
            p3(i,j)=(pp3(i)+pp3(j))/2.0
            p4(i,j)=dsqrt(pp4(i)*pp4(j))
            p5(i,j)=(pp5(i)+pp5(j))/2.0 !4.0*p2(i,j)/(pp5(i)+pp5(j))/(pp5(i)+pp5(j))
        enddo
    enddo
endif

end subroutine

subroutine potenctest()
use global
    implicit none
real(8) r
integer(4) i

open(22,file='pottest')
do i=1,500
    r=float(i)/500.0*10.0
    write(22,'(3f20.10)') r, pfunc(1,1,r), dpfunc(1,1,r)
enddo

print *, ' --- Potenc test DONE --- '
close(22)
end subroutine

subroutine totalenergy()
use global
    implicit none
integer(4) i,j,k
real(8) rmin,r1,r2,r3

totalen=0.0
totalden=0.0
totalten=0.0
do i=1,N
    do j=1,N
        if (j/=i) then
            call minobr(x(i),x(j),y(i),y(j),z(i),z(j),rmin)
            call checknear(i)
            if (rmin<rcut) then
                totalen=totalen+pfunc(i,j,rmin)
                totalden=totalden+dpfunc(i,j,rmin)
                if (calctri/=0) then
                    do k=1,N
                        call minobr(x(i),x(j),y(i),y(j),z(i),z(j),r1)
                        call minobr(x(i),x(k),y(i),y(k),z(i),z(k),r2)
                        call minobr(x(k),x(j),y(k),y(j),z(k),z(j),r3)
                        totalten=totalten+tfunc(i,j,k,r1,r2,r3)
                    enddo
                endif
            endif
        endif
    enddo
enddo
totalen=totalen/2.0
totalden=totalden/2.0
totalten=totalten/3.0
print *,'total 2-body energy', totalen
print *,'total virial ', totalden
print *,'total 3-body energy', totalten
end subroutine

subroutine minobr(x1,x2,y1,y2,z1,z2,r)
use global
    implicit none
real(8) x1,x2,y1,y2,z1,z2
real(8) r
real(8) pdx,pdy,pdz

pdx=abs(x1-x2)
if (pdx>xboxh) then
    pdx=xbox-pdx
endif
pdy=abs(y1-y2)
if (pdy>yboxh) then
    pdy=ybox-pdy
endif
pdz=abs(z1-z2)
if (pdz>zboxh) then
    pdz=zbox-pdz
endif
r=dsqrt(pdx*pdx+pdy*pdy+pdz*pdz)

end subroutine

subroutine checknear(Nmol)
use global
    implicit none
integer(4) i,Nmol
real(8) r
connected(Nmol)=0
do i=1,n
    if (i/=Nmol) then
        call minobr(x(Nmol),x(i),y(Nmol),y(i),z(Nmol),z(i),r)
        if(r<p1(tip(i),tip(Nmol))*0.75) then
            nearm(i,Nmol)=1
            nearm(Nmol,i)=1
            connected(Nmol)=1
        else
            nearm(i,Nmol)=0
            nearm(Nmol,i)=0
        endif
    endif
enddo

end subroutine

subroutine pcalc(nmol,p,dp,t,xm,ym,zm)
use global
    implicit none
real(8) p,dp,t
real(8) xm,ym,zm
integer(4) nmol
integer(4) i,j
real(8) r,r1,r2,r3

p=0.0
dp=0.0
t=0.0
do i=1,N
    if (i/=nmol) then
        call minobr(x(i),xm,y(i),ym,z(i),zm,r)
        if (r<rcut) then
            p=p+pfunc(i,nmol,r)
            dp=dp+dpfunc(i,nmol,r)
            if (calctri/=0) then
                do j=1,n
                    if ((j/=i).and.(j/=nmol)) then
                        call minobr(x(i),x(j),y(i),y(j),z(i),z(j),r1)
                        call minobr(x(i),xm,y(i),ym,z(i),zm,r2)
                        call minobr(xm,x(j),ym,y(j),zm,z(j),r3)
                        t=t+tfunc(i,j,nmol,r1,r2,r3)
                    endif
                enddo
            endif
        end if
    endif
enddo

!if (t>0) then
!    print *, p,dp,t
!endif

end subroutine


subroutine trans(nmol,xn,yn,zn)
use generator
use global
    implicit none
integer(4) nmol
real(8) xn,yn,zn
real(8) dx,dy,dz

if ((ptip==2).or.(ptip==3)) then
    if (getrand()>0.5) then
        maxdl=1.0
    else
        maxdl=0.15
    endif
endif

dx=(getrand()-0.5)*2.0*maxdl
dy=(getrand()-0.5)*2.0*maxdl
dz=(getrand()-0.5)*2.0*maxdl

xn=x(nmol)+dx
yn=y(nmol)+dy
zn=z(nmol)+dz

if (xn>xgrb) then
    xn=xn-xbox
endif
if (yn>ygrb) then
    yn=yn-ybox
endif
if (zn>zgrb) then
    zn=zn-zbox
endif
if (xn<xgrm) then
    xn=xn+xbox
endif
if (yn<ygrm) then
    yn=yn+ybox
endif
if (zn<zgrm) then
    zn=zn+zbox
endif

end subroutine

subroutine checkmove(nmol,po,dpo,to,pn,dpn,tn,xn,yn,zn)
use generator
use global
    implicit none
real(8) po,dpo,to,pn,dpn,tn,xn,yn,zn
integer(4) nmol
real(8) deltap

deltap=pn+tn-po-to
if (dexp(-deltap/temp)<getrand()) then
    !not accept
    naccept=naccept+1
else
    !accept
    accept=accept+1
    x(nmol)=xn
    y(nmol)=yn
    z(nmol)=zn
    totalen=totalen+pn-po
    totalden=totalden+dpn-dpo
    totalten=totalten+tn-to
    if (tn>0) then
        print *,totalten,tn,to
    endif

endif
end subroutine

subroutine xyzanim()
    use global
    implicit none
integer(4) i

write(31,'(i5)') N
write(31,'(a)') 'comentline'
do i=1,N
    write(31, '(a3,f10.6,f10.6,f10.6)') label(tip(i)),x(i),y(i),z(i)
enddo

end subroutine

subroutine checkdl()
use global
    implicit none
!print *, '--------------------------------------'
!print *, ' Old max distance: ', maxdl
if (accept==0) then
    maxdl=maxdl*0.75
else if (naccept==0) then
    maxdl=maxdl*1.25
else if (float(accept)/float(naccept)<0.4) then
    maxdl=maxdl*0.75
else if (float(accept)/float(naccept)>0.6) then
    maxdl=maxdl*1.25
endif
if (ptip==2) then
    maxdl=1.0*p1(1,1)
endif

!print *, ' Accept: ', accept, ' Not Accept: ', naccept
!print *, ' New max distance: ', maxdl
!print *, '------------------------------------- '
accept=0
naccept=0


end subroutine


subroutine dosample(cs)
use global
    implicit none
integer(4) cs

sumsamp(cs)=sumsamp(cs)+1
sumten(cs)=sumten(cs)+totalen
sumtden(cs)=sumtden(cs)+totalden
sumtten(cs)=sumtten(cs)+totalten


end subroutine

subroutine resout(csamp)
use global
    implicit none
integer(4) i,csamp,j,k

open(52,file='out.txt')
write(52,'(a,f20.10)') 'Temperature: ', temp
write(52,'(a,f20.10)') 'Density: ', ro
do i=1,csamp
    write(52,'(a,i5)') 'sample #', i
    write(52,'(a,f20.10)') 'Potential Energy: ', sumten(i)/float(sumsamp(i))/float(N)+ecorr
    write(52,'(a,f20.10)') 'Virial pressure: ', ro*Temp-1.0/3.0*ro/float(n)*sumtden(i)/&
    &float(sumsamp(i))+pcorr
    write(52,'(3a)') 'Coordination number ', '  by atom ', '  by molecule  '
    write(52,'(a,2f20.10)') ' 0 ', skch10(i)/kchnum(i),  skch20(i)/kchnum(i)
    write(52,'(a,2f20.10)') ' 1 ', skch11(i)/kchnum(i),  skch21(i)/kchnum(i)
    write(52,'(a,2f20.10)') ' 2 ', skch12(i)/kchnum(i),  skch22(i)/kchnum(i)
    write(52,'(a,2f20.10)') ' 3 ', skch13(i)/kchnum(i),  skch23(i)/kchnum(i)
    write(52,'(a,2f20.10)') ' 4 ', skch14(i)/kchnum(i),  skch24(i)/kchnum(i)
    write(52,'(a,2f20.10)') ' 5 ', skch15(i)/kchnum(i),  skch25(i)/kchnum(i)

enddo

close(52)
!---------------------------------------rdf1
open(53,file='rdf1.out')
write(53,'(a30,$)') 'r  '
do i=1,nv
    do j=1,nv
        write(53,'(5a,$)') ' g(', trim(adjustl(label(i))),'-',trim(adjustl(label(j))),')  '
    enddo
enddo
write(53,'(a)') ' '
do i=1,ceiling(rcut/drrdf1)
    write(53,'(f20.10,$)') drrdf1*(float(i)-0.5)
    do j=1,nv
        do k=1,nv
            write(53,'(f20.10,$)') sumrdf1(i,(j-1)*nv+k)/float(rdfnum)&
            &/idrdf1(i,(j-1)*nv+k)/float(n)
        enddo
    enddo
    write(53,'(a)')
enddo
close(53)
!----------------------------------------rdf2
open(53,file='rdf2.out')
write(53,'(a30,$)') 'r  '
do i=1,nv
    do j=1,nv
        write(53,'(5a,$)') ' g(', trim(adjustl(label(i))),'-',trim(adjustl(label(j))),')  '
    enddo
enddo
write(53,'(a)') ' '
do i=1,ceiling(rcut/drrdf2)
    write(53,'(f20.10,$)') drrdf2*(float(i)-0.5)
    do j=1,nv
        do k=1,nv
            write(53,'(f20.10,$)') sumrdf2(i,(j-1)*nv+k)/float(rdfnum)&
            &/idrdf2(i,(j-1)*nv+k)/float(n)
        enddo
    enddo
    write(53,'(a)')
enddo
close(53)
!---------------------------------------rdf3
open(53,file='rdf3.out')
write(53,'(a30,$)') 'r  '
do i=1,nv
    do j=1,nv
        write(53,'(5a,$)') ' g(', trim(adjustl(label(i))),'-',trim(adjustl(label(j))),')  '
    enddo
enddo
write(53,'(a)') ' '
do i=1,ceiling(rcut/drrdf3)
    write(53,'(f20.10,$)') drrdf3*(float(i)-0.5)
    do j=1,nv
        do k=1,nv
            write(53,'(f20.10,$)') sumrdf3(i,(j-1)*nv+k)/float(rdfnum)&
            &/idrdf3(i,(j-1)*nv+k)/float(n)
        enddo
    enddo
    write(53,'(a)')
enddo
close(53)
!!!
open(53,file='densdistr1.out')
    do i=1,indens1
        write(53,'(2f20.10)') float(i)*drdens1,sumdens1(i)/float(numdens)&
        &/(kubl*kubl*drdens1)
    enddo
close(53)
open(53,file='densdistr2.out')
    do i=1,indens2
        write(53,'(2f20.10)') float(i)*drdens2,sumdens2(i)/float(numdens)&
        &/(kubl*kubl*drdens2)
    enddo
close(53)
open(53,file='densdistr3.out')
    do i=1,indens3
        write(53,'(2f20.10)') float(i)*drdens3,sumdens3(i)/float(numdens)&
        &/(kubl*kubl*drdens3)
    enddo
close(53)
end subroutine

subroutine equilibout(nmov)
use global
    implicit none
integer(4) nmov
integer(4) kchob,histkch
real(8) kch0,kch1,kch2,kch3,kch4,kch5
real(8) kchsum1,kchsum2
integer(4) i,j

write(56,'(i10,4f20.10)') nmov,totalen,totalden,totalten, maxdl

kch0=0.0
kch1=0.0
kch2=0.0
kch3=0.0
kch4=0.0
kch5=0.0
kchob=0
kchsum1=0.0
kchsum2=0.0
do i=1,n
    histkch=0
    do j=1,n
        if (i/=j) then
            if (nearm(i,j)==1) then
                kchob=kchob+1
                histkch=histkch+1
            endif
        endif
    enddo
    if (histkch==0) then
        kch0=kch0+1.0
    endif
    if (histkch==1) then
        kch1=kch1+1.0
    endif
    if (histkch==2) then
        kch2=kch2+1.0
    endif
    if (histkch==3) then
        kch3=kch3+1.0
    endif
    if (histkch==4) then
        kch4=kch4+1.0
    endif
    if (histkch==5) then
        kch5=kch5+1.0
    endif
enddo
kchsum1=kch0+kch1+kch2+kch3+kch4+kch5
kchsum2=kch0+kch1/2.0+kch2/3.0+kch3/4.0+kch4/5.0+kch5/6.0

write(34,'(i8,13f20.15)') nmov,kch0/kchsum1,kch1/kchsum1,kch2/kchsum1,&
    &kch3/kchsum1,kch4/kchsum1,kch5/kchsum1,kch0/kchsum2,kch1/kchsum2/2.0&
    &,kch2/3.0/kchsum2,kch3/4.0/kchsum2,kch4/5.0/kchsum2,kch5/6.0/kchsum2,kchob/float(n*n)
end subroutine

subroutine tailcorr()
use global
    implicit none
integer(4) i,j

if ((ptip==1).or.(ptip==2)) then
    do i=1,nv
        do j=1,nv
            ecorr=ecorr+float(ni(i))*float(ni(j))/float(n*n)*4.0*p1(i,j)*&
                &(p1(i,j)**12/9.0/rcut**9-p1(i,j)**6/3.0/rcut**3)*2.0*pi*ro
            pcorr=pcorr+float(ni(i))*float(ni(j))/float(n*n)*4.0*p1(i,j)*&
                &(4.0*p1(i,j)**12/3.0/rcut**9-2.0*p1(i,j)**6/rcut**3)*&
                &2.0/3.0*pi*ro*ro
        enddo
    enddo
endif
print *, '2-body energy correction: ', ecorr
print *, 'pressure correction: ', pcorr
end subroutine

subroutine calcrdf()
use global
    implicit none
integer(4) i,j
real(8) r
integer(4) hist1,hist2,hist3

rdfnum=rdfnum+1
do i=1,N
    do j=1,N
        if (i/=j) then
            call minobr(x(i),x(j),y(i),y(j),z(i),z(j),r)
            if (r<rcut) then
                hist1=ceiling(r/drrdf1)
                hist2=ceiling(r/drrdf2)
                hist3=ceiling(r/drrdf3)
                sumrdf1(hist1,(tip(i)-1)*nv+tip(j))=sumrdf1(hist1,(tip(i)-1)*nv+tip(j))+1.0
                sumrdf2(hist2,(tip(i)-1)*nv+tip(j))=sumrdf2(hist2,(tip(i)-1)*nv+tip(j))+1.0
                sumrdf3(hist3,(tip(i)-1)*nv+tip(j))=sumrdf3(hist3,(tip(i)-1)*nv+tip(j))+1.0
            end if
        endif
    enddo
enddo
end subroutine

subroutine calckch(csamp)
use global
    implicit none
real(8) kch0,kch1,kch2,kch3,kch4,kch5
real(8) kchsum1,kchsum2
integer(4) histkch,i,j,kchob
integer(4) csamp



kch0=0.0
kch1=0.0
kch2=0.0
kch3=0.0
kch4=0.0
kch5=0.0
kchob=0
kchsum1=0.0
kchsum2=0.0
do i=1,n
    histkch=0
    do j=1,n
        if (i/=j) then
            if (nearm(i,j)==1) then
                kchob=kchob+1
                histkch=histkch+1
            endif
        endif
    enddo
    if (histkch==0) then
        kch0=kch0+1.0
    endif
    if (histkch==1) then
        kch1=kch1+1.0
    endif
    if (histkch==2) then
        kch2=kch2+1.0
    endif
    if (histkch==3) then
        kch3=kch3+1.0
    endif
    if (histkch==4) then
        kch4=kch4+1.0
    endif
    if (histkch==5) then
        kch5=kch5+1.0
    endif
enddo
kchsum1=kch0+kch1+kch2+kch3+kch4+kch5
kchsum2=kch0+kch1/2.0+kch2/3.0+kch3/4.0+kch4/5.0+kch5/6.0

kchnum(csamp)=kchnum(csamp)+1
skch10(csamp)=skch10(csamp)+kch0/kchsum1
skch11(csamp)=skch11(csamp)+kch1/kchsum1
skch12(csamp)=skch12(csamp)+kch2/kchsum1
skch13(csamp)=skch13(csamp)+kch3/kchsum1
skch14(csamp)=skch14(csamp)+kch4/kchsum1
skch15(csamp)=skch15(csamp)+kch5/kchsum1

skch20(csamp)=skch20(csamp)+kch0/kchsum2
skch21(csamp)=skch21(csamp)+kch1/kchsum2/2.0
skch22(csamp)=skch22(csamp)+kch2/kchsum2/3.0
skch23(csamp)=skch23(csamp)+kch3/kchsum2/4.0
skch24(csamp)=skch24(csamp)+kch4/kchsum2/5.0
skch25(csamp)=skch25(csamp)+kch5/kchsum2/6.0

end subroutine

subroutine mod1trans(Nmol,xn,yn,zn)
use global
use generator
    implicit none
integer(4) Nmol,i
real(8) xn,yn,zn
real(8) randalfa,randbeta,randgamma
integer(4) nearmol
real(8) xtemp,ytemp,ztemp
real(8) xv,yv,zv
real(8) cosox,cosoy,cosoz
real(8) sinox,sinoy,sinoz
real(8) rast_vr
real(8) dx,dy,dz

if (connected(Nmol)==1) then
do i=1,N
    if (nearm(Nmol,i)>0.5) then
        nearmol=i
    endif
enddo
randalfa=(getrand()-0.5)*0.3
randbeta=(getrand()-0.5)*0.3
randgamma=(getrand()-0.5)*0.3
!перемещаем молекулы так чтобы вокруг которой была в центре
xtemp=x(Nmol)-x(nearmol)
ytemp=y(Nmol)-y(nearmol)
ztemp=z(Nmol)-z(nearmol)
!находим растояние вокруг которого надо вращать
rast_vr=1.0+(getrand()-0.5)*0.005
!добавляем случайное изменение
xtemp=xtemp*rast_vr
ytemp=ytemp*rast_vr
ztemp=ztemp*rast_vr
!задаем углы поворотов
sinox=dsin(randalfa)  !0.4*(randalfa-0.5)*2.0
cosox=dcos(randalfa)  !dsqrt(1.0-sinox*sinox)
sinoy=dsin(randbeta)  !0.4*(randbeta-0.5)*2.0
cosoy=dcos(randbeta)  !dsqrt(1.0-sinoy*sinoy)
sinoz=dsin(randgamma)  !0.4*(randgamma-0.5)*2.0
cosoz=dcos(randgamma)  !dsqrt(1.0-sinoz*sinoz)
xv=xtemp
yv=ytemp
zv=ztemp
!поворачиваем
xtemp=cosox*cosoy*xv+(cosox*sinoy*sinoz-sinox*cosoz)*yv+(cosox*sinoy*cosoz+sinox*sinoz)*zv
ytemp=sinox*cosoy*xv+(sinox*sinoy*sinoz+cosox*cosoz)*yv+(sinox*sinoy*cosoz-cosox*sinoz)*zv
ztemp=-sinoy*xv+cosoy*sinoz*yv+cosoy*cosoz*zv
!вставляем моекулу назад
xn=xtemp+x(nearmol)
yn=ytemp+y(nearmol)
zn=ztemp+z(nearmol)
else
    if (ptip==2) then
        if (getrand()>0.5) then
            maxdl=1.0
        else
            maxdl=0.15
        endif
    endif
    dx=(getrand()-0.5)*maxdl
    dy=(getrand()-0.5)*maxdl
    dz=(getrand()-0.5)*maxdl
    xn=x(Nmol)+dx
    yn=y(Nmol)+dy
    zn=z(Nmol)+dz
endif

if (xn>xgrb) then
    xn=xn-xbox
endif
if (yn>ygrb) then
    yn=yn-ybox
endif
if (zn>zgrb) then
    zn=zn-zbox
endif
if (xn<xgrm) then
    xn=xn+xbox
endif
if (yn<ygrm) then
    yn=yn+ybox
endif
if (zn<zgrm) then
    zn=zn+zbox
endif

end subroutine

subroutine calcdens()
use global
    implicit none
integer(4) i
integer(4) hist1,hist2,hist3

do i=1,n
    hist1=ceiling((z(i)-zgrm)/drdens1)
    hist2=ceiling((z(i)-zgrm)/drdens2)
    hist3=ceiling((z(i)-zgrm)/drdens3)
    sumdens1(hist1)=sumdens1(hist1)+1.0
    sumdens2(hist2)=sumdens2(hist2)+1.0
    sumdens3(hist3)=sumdens3(hist3)+1.0
enddo
numdens=numdens+1

end subroutine

subroutine densinit()
use global
    implicit none
!integer(4) i

indens1=30
indens2=60
indens3=90

allocate(sumdens1(indens1))
allocate(sumdens2(indens2))
allocate(sumdens3(indens3))

sumdens1=0.0
sumdens2=0.0
sumdens3=0.0

drdens1=zbox/float(indens1)
drdens2=zbox/float(indens2)
drdens3=zbox/float(indens3)
numdens=0

end subroutine
