program forced_oscillator
real::x,v,t,tf,h,fx,fv,k1,k2,k3,k4
x=0.5;v=0.0;t=0.0;tf=15.0;h=0.01;m=1.0
open(1,file="forced_displacement.dat",status='unknown')
do while (t<=tf)
    write(1,*)t,x
    call rk4(x,v,t,h)
end do

end program

real function fx(t,v)
real::x,v,t
fx=v
return
end function

real function fv(t,x,v)
real::x,v,t,k=80.0,b=1.0,m=1.0,f=80.0
fv=-(k/m)*x-2*b*v+f*cos(2*sqrt((k/m)-b**2)*t)
return
end function

subroutine rk4(x,v,t,h)
    real::k1,k2,k3,k4,dx,dv,x,v,t,h
    k1=h*fv(t,x,v)
    k2=h*fv(t+0.5*h,x+0.5*k1,v+0.5*k1)
    k3=h*fv(t+0.5*h,x+0.5*k2,v+0.5*k2)
    k4=h*fv(t+h,x+k3,v+k3)
    dv=(1/6.)*(k1+k4+2*(k2+k3))
    v=v+dv
    k1=h*fx(t,v)
    k2=h*fx(t+0.5*h,v+0.5*k1)
    k3=h*fx(t+0.5*h,v+0.5*k2)
    k4=h*fx(t+h,v+k3)
    dx=(1/6.)*(k1+k4+2*(k2+k3))
    x=x+dx
    t=t+h
    print*,x
    return
end subroutine
