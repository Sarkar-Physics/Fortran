program freefall
    real::t,v,h,f,tf
    print *,"Enter initial time(sec):"
    read(*,*)t
    print *,"Enter initial velocity(m/s):"
    read(*,*)v
    print *,"Enter the incriment:"
    read(*,*)h
    print *,"Enter final value of time(sec):"
    read(*,*)tf
    !open(1,file='freefall.dat',status='unknown')
    do while (t<=tf)
        call dvdt(f,t)
        print *,t,v
        !write(1,*)t,v
        t=t+h
        v=v+h*f
    end do
close(1)
end program freefall
subroutine dvdt(f,t)
    f=9.8
    return
end subroutine
