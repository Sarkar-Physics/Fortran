! program to find the eigen function and eigen energ of a particle in infinite square well
program infinite_square_well
real*8::psip,psif,dpsip,sf,h,sp,E,E_min,E_max,psi1,psi2,norm,pi
real*8,allocatable::psi(:),dpsi(:),s(:)
integer::n,node
! The value of the wave function at s=0 , psip and s=1, psif  is 0 and the derivative of psi is set as 1.0
! step size is taken as 0.01
psip=0.0;psif=0.0;dpsip=1.0;sp=0.0;sf=1.0;h=1e-6;node=0;pi=acos(-1.0)
!initial guess to determine energy eigen value
E_min=100.0;E_max=150.0
! no of discrete data points
n=(sf-sp)/h
allocate(psi(0:n),dpsi(0:n),s(0:n))
! setting the initial value of the first element of array
psi(0)=0.0;dpsi(0)=1.0;s(0)=0.0
write(*,*)"Minimum Energy:",E_min,"Maximum Energy:",E_max
open(1,file='infinite_square_well.dat',status='unknown')
!Checking the value of psi at s=1 for Minimum Energy
do i=1,n
    call rk4(psip,dpsip,sp,sf,E_min,h)
    psi(i)=psip
    psi1=psip
    s(i)=sp
    end do
! Counting the number of Nodes for initial trial Energy
do i=1,n-1
    if(psi(i)/=sign(psi(i),psi(i-1))) node=node+1
end do
psip=0.0;psif=0.0;dpsip=1.0;sp=0.0;sf=1.0
! checking the value of psi at s=1 for Maximum Energy
do i=1,n
    call rk4(psip,dpsip,sp,sf,E_max,h)
    psi2=psip
    end do
! If the product of value of Psi for two different values of energy is positive then the bisection method wont work
if(psi1*psi2>0) then
    write(*,*)"!!! Wrong initial guess value of Energy given !!!"
    stop
end if
! Calling the subroutine energy_eigen_value to find the approximate energy using bisection method
call energy_eigen_value(E_min,E_max,E,node)
psi(0)=0.0;dpsi(0)=1.0;s(0)=0.0
psip=0.0;psif=0.0;dpsip=1.0;sp=0.0;sf=1.0;node=0
!Calculating the Wave function using the energy found from bisection method
do i=1,n
    call rk4(psip,dpsip,sp,sf,E,h)
    psi(i)=psip
    s(i)=sp
end do
! Finding the nodes and the quantum number of state
do i=1,n-1
    if(psi(i)/=sign(psi(i),psi(i-1))) node=node+1
end do
write(*,*)"Nodes:",node,"Energy Eigen State:",node+1
write(*,*)"Eigen Energy:",E,"Theoritical Value of Eigen Energy:",(node+1)**2*pi**2/2.0
! Calling the subroutine normalization to find the normalized wave function using simpson 1/3 rd rule
call normalization(E,norm)
! Writing the wave function into a file for plotting
do i=0,n
    write(1,*)s(i),psi(i)/sqrt(norm),sqrt(2.0)*sin((node+1)*pi*s(i))
end do
close(1)
end program
! Defining the Function for solving the Schrodinger Equation
real*8 function f(counter,E,sp,dpsip,psip)
integer::counter
real*8::psip,dpsip,sp,E,V
V=0.0
select case(counter)
case(1)
    f=dpsip         !f=d/ds(psi)
case(2)
    f=2*(V-E)*psip     !f=d/ds(dpsi)
end select
end function
! Subroutine rk4 to calculate the Wave function using Runge-Kutta 4th order method
subroutine rk4(psip,dpsip,sp,sf,E,h)
real*8::psip,dpsip,sp,sf,h,E
real*8::k1_psi,k2_psi,k3_psi,k4_psi,k1_dpsi,k2_dpsi,k3_dpsi,k4_dpsi,f
    ! solving psi using runge kutta 4th order algorithm
    k1_psi=h*f(1,E,sp,dpsip,psip)
    k2_psi=h*f(1,E,sp+0.5*h,dpsip,psip+0.5*k1_psi)
    k3_psi=h*f(1,E,sp+0.5*h,dpsip,psip+0.5*k2_psi)
    k4_psi=h*f(1,E,sp+h,dpsip,psip+k3_psi)
    psip=psip+(1/6.0)*(k1_psi+k4_psi+2*(k2_psi+k3_psi))
    ! solving for dpsi
    k1_dpsi=h*f(2,E,sp,dpsip,psip)
    k2_dpsi=h*f(2,E,sp+0.5*h,dpsip+0.5*k1_dpsi,psip)
    k3_dpsi=h*f(2,E,sp+0.5*h,dpsip+0.5*k2_dpsi,psip)
    k4_dpsi=h*f(2,E,sp+h,dpsip+k3_dpsi,psip)
    dpsip=dpsip+(1/6.0)*(k1_dpsi+k4_dpsi+2*(k2_dpsi+k3_dpsi))
    sp=sp+h
return
end subroutine
! Subroutine energy_eigen_value to find the approximate energy eigen value using bisection method
subroutine energy_eigen_value(E_min,E_max,E,node)
real*8::psip,psif,dpsip,sp,sf,h,last,E,E_min,E_max
integer::n,node
psip=0.0;psif=0.0;dpsip=1.0;sp=0.0;sf=1.0;h=1e-6
E=E_min
n=(sf-sp)/h
    do i=1,n
        call rk4(psip,dpsip,sp,sf,E,h)
    end do
do while(dabs(psif-psip)>1e-10)
    E=(E_min+E_max)/2.0
    psip=0.0;psif=0.0;dpsip=1.0;s0=0.0;sf=1.0;h=1e-6;
    do j=1,n
    call rk4(psip,dpsip,sp,sf,E,h)
    end do
    if(modulo(node,2)==0) then
        if(psip>0) then
            E_min=E
        else
            E_max=E
        end if
    else
        if(psip<0) then
            E_min=E
        else
            E_max=E
        end if
    end if
end do
end subroutine
! Subroutine normalization to find the normalized wave function using Simpson's 1/3 rd Rule
subroutine normalization(E,norm)
    real*8::psip,psif,dpsip,sp,sf,h,t1,t2,norm,E
    real*8,allocatable::psi(:),dpsi(:),s(:)
    psip=0.0;psif=0.0;dpsip=1.0;sp=0.0;s0=0.0;sf=1.0;h=1e-6
    t1=0.0;t2=0.0
    n=(sf-sp)/h
    allocate(psi(0:n),dpsi(0:n),s(0:n))
    psi(0)=0.0;dpsi(0)=1.0;s(0)=0.0
    do i=1,n
        call rk4(psip,dpsip,sp,sf,E,h)
        psi(i)=psip
        s(i)=sp
    end do
    do i=1,n-1,2
        t1=t1+psi(i)**2
    end do
    do i=2,n-2,2
        t2=t2+psi(i)**2
    end do
    norm=(h/3.0)*(psi(0)**2+psi(n)**2+4.0*t1+2.0*t2)
end subroutine
