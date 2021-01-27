program airdrag_m1
    real::t,v,h,dvdt,tf,a,v_analy,g,m,z
    print *,"Enter initial time(sec):"
    read *,t
    print *,"Enter initial velocity(m/s):"
    read *,v
    print *,"Enter the incriment:"
    read *,h
    print *,"Enter final value of time(sec):"
    read *,tf
    print *,"Enter the value of g(m/s^2):"
    read *,g
    print *,"Enter mass of the object(kg):"
    read *,m
    print *, "Enter the drag coefficient(kg/m):"
    read *,a
    v_analy=v
    z=a/m
    !open(1,file='air_drag_m1.dat',status='unknown')
    !open(2,file='analytical_m1.dat',status='unknown')
    do while (t<=tf)
        call f(dvdt,g,v,z)
        print *,t,v
        print *,t,v_analy
        !write(1,*)t,v
        !write(2,*)t,v_analy
        t=t+h
        v=v+h*dvdt
        v_analy=((9.8/z)*(1-exp(-z*t)))
    end do
close(1)
close(2)

end program airdrag_m1
subroutine f(dvdt,g,v,z)
    dvdt=g-(z*v)
    return
end subroutine
