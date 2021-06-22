program projectile
real::x,vx,y,vy,t,h,pi,theta
pi=acos(-1.0)
theta=19.7678280*pi/180.0
x=0.0;y=0.0;vx=650*cos(theta);vy=650*sin(theta);t=0.0;h=0.01;m=1.0
open(1,file="projectile2_rk4.dat",status='unknown')
do while (y>=0.0)
    write(1,*)x,y
    call rk4(x,y,vx,vy,t,h)
end do
print*,"Range of projectile(in m):",x
end program

real function fx(t,vx)
real::x,vx,t
fx=vx
end function

real function fvx(t,vx,vy)
real::x,vx,vy,t
fvx=-0.0000389*vx*sqrt(vx**2+vy**2)
end function

real function fvy(t,vx,vy)
real::y,vx,vy,t
fvy=-9.8-0.0000389*vy*sqrt(vx**2+vy**2)
end function

real function fy(t,vy)
real::y,vy,t
fy=vy
end function

subroutine rk4(x,y,vx,vy,t,h)
    real::k1,k2,k3,k4,x,dx,y,dy,vx,dvx,vy,dvy
    k1=h*fvx(t,vx,vy)
    k2=h*fvx(t+0.5*h,vx+0.5*k1,vy+0.5*k1)
    k3=h*fvx(t+0.5*h,vx+0.5*k2,vy+0.5*k2)
    k4=h*fvx(t+h,vx+k3,vy+k3)
    dvx=(1/6.)*(k1+k4+2*(k2+k3))
    vx=vx+dvx
    k1=h*fx(t,vx)
    k2=h*fx(t+0.5*h,vx+0.5*k1)
    k3=h*fx(t+0.5*h,vx+0.5*k2)
    k4=h*fx(t+h,vx+k3)
    dx=(1/6.)*(k1+k4+2*(k2+k3))
    x=x+dx
    k1=h*fvy(t,vx,vy)
    k2=h*fvy(t+0.5*h,vx+0.5*k1,vy+0.5*k1)
    k3=h*fvy(t+0.5*h,vx+0.5*k2,vy+0.5*k2)
    k4=h*fvy(t+h,vx+k3,vy+k3)
    dvy=(1/6.)*(k1+k4+2*(k2+k3))
    vy=vy+dvy
    k1=h*fy(t,vy)
    k2=h*fy(t+0.5*h,vy+0.5*k1)
    k3=h*fy(t+0.5*h,vy+0.5*k2)
    k4=h*fy(t+h,vy+k3)
    dy=(1/6.)*(k1+k4+2*(k2+k3))
    y=y+dy
    t=t+h
    return
end subroutine
