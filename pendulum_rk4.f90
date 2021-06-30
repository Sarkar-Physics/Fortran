program harmonic_oscillator
real::x,v,t,tf,h
x=0.5;v=0.0;t=0.0;tf=5.0;h=0.01;m=1.0;k=80.0
!x is the displacement at t=0
!v is the velocity at t=0
!tf is the final time
!h is the step size
!m is the mass of the object and k is the spring constant
open(1,file="pendulum.dat",status='unknown')
open(2,file='pendulum_analytical.dat',status='unknown')
call rk4(x,v,t,tf,h)
! calling the subroutine to solve the differential equation
close(1)
close(2)
end program

real function f(counter,t,v,x)
integer::counter
real::x,v,t,k=80.0,m=1.0
select case(counter)
case(1)
    f=v         !f=dx/dt
case(2)
    f=-(k/m)*x  !f=dv/dt
end select
end function

subroutine rk4(x,v,t,tf,h)
    real::k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,x,v,t,tf,h,f,k,m
    k=80.0;m=1.0
    ! solving x and v using runge kutta 4th order algorithm
    do while (t<=tf)
        write(1,*)t,x,m*v,0.5*m*v**2,0.5*k*x**2,0.5*(m*v**2+k*x**2)
        analy_x=0.5*cos(sqrt(k/m)*t)
        analy_v=-0.5*sqrt(k/m)*sin(sqrt(k/m)*t)
        write(2,*)t,analy_x,m*analy_v,0.5*m*analy_v**2,0.5*k*analy_x**2,0.5*(m*analy_v**2+k*analy_x**2)
        !solve for x
        k1x=h*f(1,t,v,x)
        k2x=h*f(1,t+0.5*h,v,x+0.5*k1x)
        k3x=h*f(1,t+0.5*h,v,x+0.5*k2x)
        k4x=h*f(1,t+h,v,x+k3x)
        x=x+(1/6.0)*(k1x+k4x+2*(k2x+k3x))
        ! solve for v
        k1v=h*f(2,t,v,x)
        k2v=h*f(2,t+0.5*h,v+0.5*k1v,x)
        k3v=h*f(2,t+0.5*h,v+0.5*k2v,x)
        k4v=h*f(2,t+h,v+k3v,x)
        v=v+(1/6.0)*(k1v+k4v+2*(k2v+k3v))
        t=t+h
        !print*,t,x,analy_x
    end do
    return
end subroutine
