program decay_LR
    real::t,i,h,f,tf,i_analy,V,R,L,prev_i,new_i
    t=0.0
    tf=2.0
    i=10.0
    h=0.1
    i_analy=i
    open(1,file='decay.dat',status='unknown')
    open(2,file='analytical_decay.dat',status='unknown')
    do while (t<=tf)
        call didt(f,i,R,L)
        print *,t,i,i_analy
        write(1,*)t,i
        write(2,*)t,i_analy
        prev_i=i
        t=t+h
        i=i+h*f
        new_i=i
        i_analy=10.0*(exp(-5.0*t))
        if(prev_i-new_i<=0.0001) then
        exit
        end if
    end do
close(1)
close(2)

end program decay_LR
subroutine didt(f,i,R,L)
    real::f,R,L,i
    R=2.0
    L=0.4
    f=-(R/L)*i
    return
end subroutine
