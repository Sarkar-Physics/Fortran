program newton_cooling
    real::T_0,t,tf,h,temp,f,x,y
    T_0=1200.0
    t=0.0
    tf=8000.0
    h=80.0
    x=0
    y=T_0
    temp=T_0
    open(1,file="newton_cooling_3.dat",status='unknown')
    !open(2,file="analytical_cooling.dat",status='unknown')
    do while(t<=tf)
        call cooling(f,temp)
        !print *,t,temp
        write(1,*)t,temp
        !write(2,*)x,y
        t=t+h
        temp=temp+h*f
        !x=((1.8519*atan(0.00333*y)-0.92593*log((y-300)/(y+300))-2.9282)/(0.2267*10.**(-3)))
        !y=y-h
    end do
end program
subroutine cooling(f,temp)
    real::f,temp
    f=-2.2067*10.**(-12)*(temp**4-81*10.**8)
    return
end subroutine
