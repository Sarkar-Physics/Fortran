program parachute
    real::t,v,h,dvdt,tf,a,g,m,z,v_0,T_0,t_new,v_new
    print *,"Enter initial time(sec):"
    read *,t
    print *,"Enter initial velocity(m/s):"
    read *,v
    print *,"Enter final value of time(sec):"
    read *,tf
    print *,"Enter terminal velocity:"
    read *,v_0
    print *,"Enter the value of g(m/s^2):"
    read *,g
    print *,"Enter mass of the object(kg):"
    read *,m
    print *, "Enter the drag coefficient b(kg/m):"
    read *,b
    z=b/m
    open(1,file='parachute.dat',status='unknown')
    do while(v<=v_0)
          call airdrag1(dvdt,g,v)
          h=0.1
          t=t+h
          v=v+h*dvdt
          print*,t,v
          write(1,*)t,v
    end do
    t_new=t
    v_new=v
    do while (t_new<=tf)
        call airdrag2(f,v_new,g,z)
        h=0.001
        print *,t_new,v_new
        write(1,*)t_new,v_new
        t_new=t_new+h
        v_new=v_new+h*f
    end do
close(1)
close(2)

end program parachute
subroutine airdrag2(f,v,g,z)
    f=g-z*v**2
    return
end subroutine
subroutine airdrag1(dvdt,g,v)
    dvdt=g-0.22*v
    return
end subroutine

