program projectile_shooting
real::x,y,xf,t,h,pi,theta_initial,theta_final,theta_mid
pi=acos(-1.0)
theta_initial=10.0*pi/180.0;theta_final=30.0*pi/180.0
x=0.0;y=0.0;xf=16769.0;t=0.0;h=0.01
call shooting(x,y,t,h,theta_initial,theta_final,theta_mid)
print*,"Required angle(in degree):",theta_mid*180.0/pi

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

subroutine shooting(x0,y0,t0,step,initial_parameter,final_parameter,theta_mid)
    real::x,y,vx,vy,xf,t,h,initial_parameter,final_parameter,theta_mid,theta1,theta2
    theta1=initial_parameter
    theta2=final_parameter
    theta_mid=theta1
    do i=0,50
        x=x0;y=y0;t=t0;h=step
        xf=16769.0
        vx=650.0*cos(theta_mid);vy=650.0*sin(theta_mid)
        call rk4(x,y,vx,vy,t,h)
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
