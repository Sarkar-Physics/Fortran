program damped_oscillator
real::x,v,t,tf,h
x=0.5;v=0.0;t=0.0;tf=10.0;h=0.01;m=1.0
open(1,file="damped_oscillator_1.dat",status='unknown')
open(2,file="damped_oscillator_2.dat",status='unknown')
open(3,file="damped_oscillator_3.dat",status='unknown')
call rk4(x,v,t,tf,h)
close(1)
close(2)
close(3)
end program

real function f(counter,b,t,v,x)
integer::counter,i
real::x,v,t,k=80.0,m=1.0,b
select case(counter)
case(1)
    f=v                     !f=dx/dt
case(2)
    f=-2.0*b*v-(k/m)*x        !f=dv/dt
end select
end function f

subroutine rk4(x,v,t,tf,h)
    integer::i
    real::k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,x,v,t,h,f,k=80.0,m=1.0
    real::b(3)
    b(1)=0.0;b(2)=1.0;b(3)=abs(sqrt(k/m))
    do i=1,3
        t=0.0;x=0.5;v=0.0
        do while (t<=tf)
            write(i,*)t,x,m*v,0.5*m*v**2,0.5*80.0*x**2,0.5*(m*v**2+80.0*x**2)
            !print*,i,b(i),t,x
            !solve for v
            k1v=h*f(2,b(i),t,v,x)
            k2v=h*f(2,b(i),t+0.5*h,v+0.5*k1v,x)
            k3v=h*f(2,b(i),t+0.5*h,v+0.5*k2v,x)
            k4v=h*f(2,b(i),t+h,v+k3v,x)
            v=v+(1/6.0)*(k1v+k4v+2*(k2v+k3v))
            ! solve for x
            k1x=h*f(1,b(i),t,v,x)
            k2x=h*f(1,b(i),t+0.5*h,v,x+0.5*k1x)
            k3x=h*f(1,b(i),t+0.5*h,v,x+0.5*k2x)
            k4x=h*f(1,b(i),t+h,v,x+k3x)
            x=x+(1/6.0)*(k1x+k4x+2*(k2x+k3x))
            t=t+h
        end do
    end do
    return
end subroutine
