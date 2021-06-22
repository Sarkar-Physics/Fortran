program projectile_shooting
real::x,y,xf,t,h,pi,theta_initial,theta_final,theta_mid
pi=acos(-1.0)
theta_initial=10.0*pi/180.0;theta_final=30.0*pi/180.0
x=0.0;y=0.0;xf=16769.0;t=0.0;h=0.01
call shooting(x,y,t,h,theta_initial,theta_final,theta_mid)
print*,"Required angle(in degree):",theta_mid*180.0/pi

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

subroutine shooting(x0,y0,t0,step,initial_parameter,final_parameter,theta_mid)
    real::x,y,vx,vy,xf,t,h,initial_parameter,final_parameter,theta_mid,theta1,theta2
    theta1=initial_parameter
    theta2=final_parameter
    theta_mid=theta1
    do i=0,50
        x=x0;y=y0;t=t0;h=step
        xf=16769.0
        vx=650.0*cos(theta_mid);vy=650.0*sin(theta_mid)
        do while (y>=0.0)
            call rk4(x,y,vx,vy,t,h)
        end do
        if(x>xf) then
            theta2=theta_mid
        elseif(x-xf<0) then
            theta1=theta_mid
         end if
        theta_mid=(theta1+theta2)/2.0
        print*,x,theta_mid*180/3.14
        if (abs(x-xf)<=0.1) then
            exit
        end if
    end do
    return
end subroutine
