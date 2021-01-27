program airdrag_1
    real::t,v,h,dvdt,tf,v0,a,v_analy,g
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
    print *, "Enter the constant a:"
    read *,a
    v_analy=v
    !open(1,file='air_drag.dat',status='unknown')
    !open(2,file='analytical_k_1.dat',status='unknown')
    do while (t<=tf)
        call f(dvdt,g,a,v)
        print *,t,v
        print *,t,v_analy
        !write(1,*)t,v
        !write(2,*)t,v_analy
        t=t+h
        v=v+h*dvdt
        v_analy=((9.8/a)*(1-exp(-a*t)))
    end do
close(1)
close(2)

end program airdrag_1
subroutine f(dvdt,g,a,v)
    dvdt=g-a*v
    return
end subroutine
