program airdrag_v2
    real::t,v,h,dvdt,tf,a,v_analy,g,m,z,v_0,T_0
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
    print *, "Enter the drag coefficient b(kg/m):"
    read *,b
    v_analy=v
    z=b/m
    v_0=sqrt((m*g)/b)
    T_0=sqrt(m/(g*b))
    print *,"Terminal velocity:",v_0
    open(1,file='air_drag_v2.dat',status='unknown')
    open(2,file='analytical_v2.dat',status='unknown')
    do while (t<=tf)
        call f(dvdt,g,v,z)
        print *,t,v
        print *,t,v_analy
        write(1,*)t,v
        write(2,*)t,v_analy
        t=t+h
        v=v+h*dvdt
        v_analy=v_0*tanh(t/T_0)
    end do
close(1)
close(2)

end program airdrag_v2
subroutine f(dvdt,g,v,z)
    dvdt=g-z*v**2
    return
end subroutine
