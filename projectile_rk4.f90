program projectile
real::x,vx,y,vy,t,h,pi,theta
pi=acos(-1.0)
theta=45*pi/180.0
x=0.0;y=0.0;vx=650.0*cos(theta);vy=650.0*sin(theta);t=0.0;h=0.01;m=1.0
open(1,file="proj.dat",status='unknown')
call rk4(x,y,vx,vy,t,h)
print*,"Range of projectile(in m):",x
end program

real function f(counter,t,vx,vy)
integer::counter
real::x,y,vx,vy,t,g,b
g=9.8;b=0.0000389
select case(counter)
case(1)
    f=vx                                      !f=dx/dt
case(2)
    f=-b*vx*sqrt(vx**2+vy**2)       !f=dvx/dt
case(3)
    f=vy                                  !f=dy/dt
case(4)
    f=-g-b*vy*sqrt(vx**2+vy**2) !f=dvy/dt
end select
end function f

subroutine rk4(x,y,vx,vy,t,h)
    real::k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1vx,k2vx,k3vx,k4vx,k1vy,k2vy,k3vy,k4vy,x,y,vx,vy,f
    do while (y>=0.0)
    write(1,*)x,y
    k1vx=h*f(2,t,vx,vy)
    k2vx=h*f(2,t+0.5*h,vx+0.5*k1vx,vy)
    k3vx=h*f(2,t+0.5*h,vx+0.5*k2vx,vy)
    k4vx=h*f(2,t+h,vx+k3vx,vy)
    vx=vx+(1/6.0)*(k1vx+k4vx+2*(k2vx+k3vx))
    k1x=h*f(1,t,vx,vy)
    k2x=h*f(1,t+0.5*h,vx+0.5*k1x,vy)
    k3x=h*f(1,t+0.5*h,vx+0.5*k2x,vy)
    k4x=h*f(1,t+h,vx+k3x,vy)
    x=x+(1/6.0)*(k1x+k4x+2*(k2x+k3x))
    k1vy=h*f(4,t,vx,vy)
    k2vy=h*f(4,t+0.5*h,vx,vy+0.5*k1vy)
    k3vy=h*f(4,t+0.5*h,vx,vy+0.5*k2vy)
    k4vy=h*f(4,t+h,vx,vy+k3vy)
    vy=vy+(1/6.0)*(k1vy+k4vy+2*(k2vy+k3vy))
    k1y=h*f(3,t,vx,vy)
    k2y=h*f(3,t+0.5*h,vx,vy+0.5*k1y)
    k3y=h*f(3,t+0.5*h,vx,vy+0.5*k2y)
    k4y=h*f(3,t+h,vx,vy+k3y)
    y=y+(1/6.0)*(k1y+k4y+2*(k2y+k3y))
    t=t+h
    end do
    return
end subroutine
