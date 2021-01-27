program growth_decay
    real::t,tf,i_growth,i_decay,V,R,L,f,h,prev_i,new_i,t_new
    t=0.0
    tf=3.0
    h=0.1
    i_growth=0.0
    V=20.0
    R=2.0
    L=0.4
    t_new=0
    open(1,file="growth1.dat",status='unknown')
    open(2,file="decay1.dat",status='unknown')
    do while(t<=tf)
        call growth(f,i_growth,V)
        print *,t,i_growth
        write(1,*)t,i_growth
        t=t+h
        prev_i=i_growth
        i_growth=i_growth+h*f
        new_i=i_growth
        if (new_i-prev_i<=0.0001) then
        exit
        end if
    end do
    print *,"last value of current:", new_i, "at time:",t
    tf=t
    do while(t_new<tf)
        V=0
        print *,t_new,new_i
        write(2,*)t_new,new_i
        call growth(f,new_i,V)
        t_new=t_new+h
        new_i=new_i+h*f
        end do
end program
subroutine growth(f,i,V)
    real::f,i,V,R,L
    R=2.0
    L=0.4
    f=(V-i*R)/L
    return
end subroutine
