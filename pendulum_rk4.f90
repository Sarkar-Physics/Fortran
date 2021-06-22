program harmonic_oscillator
real::x,v,t,tf,h,fx,fv,k1,k2,k3,k4,m
x=0.5;v=0.0;t=0.0;tf=5.0;h=0.01;m=1.0
open(1,file="pendulum_displacement.dat",status='unknown')
open(2,file="pendulum_phase_space.dat",status='unknown')
open(3,file="pendulum_K.E.dat",status='unknown')
open(4,file="pendulum_P.E.dat",status='unknown')
open(5,file="pendulum_T.E.dat",status='unknown')
do while (t<tf)
    !write(1,*)t,x
    print*,t,x
    !write(2,*)m*v,x
    !write(3,*)t,0.5*m*v**2
    !write(4,*)t,0.5*10.0*x**2
    !write(5,*)t,0.5*(m*v**2+10.0*x**2)
    call rk4(x,v,t,h)
end do

end program

real function fx(t,v)
real::x,v,t
fx=v
end function

real function fv(t,x)
real::x,v,t,k=10.0
fv=-k*x
end function

subroutine rk4(x,v,t,h)
    real::k1,k2,k3,k4,dx,dv,x,v,t,h
    k1=h*fv(t,x)
    k2=h*fv(t+0.5*h,x+0.5*k1)
    k3=h*fv(t+0.5*h,x+0.5*k2)
    k4=h*fv(t+h,x+k3)
    dv=(1/6.)*(k1+k4+2*(k2+k3))
    v=v+dv
    k1=h*fx(t,v)
    k2=h*fx(t+0.5*h,v+0.5*k1)
    k3=h*fx(t+0.5*h,v+0.5*k2)
    k4=h*fx(+h,v+k3)
    dx=(1/6.)*(k1+k4+2*(k2+k3))
    x=x+dx
    t=t+h
    return
end subroutine
