program forced_oscillator
real*8::x,v,t,tf,h
x=0.5;v=0.0;t=0.0;tf=10.0;h=0.01
open(1,file="try.dat",status='unknown')
call rk4(t,tf,v,x,h)
close(1)
end program

real*8 function f(counter,t,v,x)
integer::counter
real*8::x,v,t,k=80.0,b=1.0,m=1.0,w=80.0
select case(counter)
case(1)
    f=v     !f=dx/dt
case(2)
    f=-(k/m)*x-2*b*v+w*cos(2*sqrt((k/m)-b**2)*t)  !f=dv/dt
end select
end function f

subroutine rk4(t,tf,v,x,h)
    real*8::k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,x,v,t,tf,h,f
    do while (t<=tf)
        write(1,*)t,x
        !print*,t,x
        !solve for v
        k1v=h*f(2,t,v,x)
        k2v=h*f(2,t+0.5*h,v+0.5*k1v,x)
        k3v=h*f(2,t+0.5*h,v+0.5*k2v,x)
        k4v=h*f(2,t+h,v+k3v,x)
        v=v+(1/6.0)*(k1v+k4v+2*(k2v+k3v))
        ! solve for x
        k1x=h*f(1,t,v,x)
        k2x=h*f(1,t+0.5*h,v,x+0.5*k1x)
        k3x=h*f(1,t+0.5*h,v,x+0.5*k2x)
        k4x=h*f(1,t+h,v,x+k3x)
        x=x+(1/6.0)*(k1x+k4x+2*(k2x+k3x))
        t=t+h
        print*,t,x
    end do
    return
end subroutine
