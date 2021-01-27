program euler_LR
    real::t,i,h,f,tf,V,i_analy,R,L,prev_i,new_i
    h=0.1
    i_analy=i
    t=0.0
    i=0.0
    tf=2.0
    open(1,file='growth_LR.dat',status='unknown')
    open(2,file='analytical_LR.dat',status='unknown')
    do while (t<=tf)
        call didt(f,i,V,R,L)
        print *,t,i,i_analy
        write(1,*)t,i
        write(2,*)t,i_analy
        prev_i=i
        i=i+h*f
        t=t+h
        new_i=i
        i_analy=((V/R)*(1-exp(-(R/L)*t)))
        if (new_i-prev_i<=0.0001) then
        exit
        end if
    end do
close(1)
close(2)
print *,"Final value of current:",i
end program euler_LR
subroutine didt(f,i,V,R,L)
    real::f,i,V,R,L
    V=20.0
    R=2.0
    L=0.4
    f=(V-R*i)/L
    return
end subroutine
