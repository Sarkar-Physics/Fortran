program timeconstant
    real::t,i,h,f,tf,V,R,L,prev_i,new_i,e,time
    h=0.01
    t=0.0
    i=0.0
    tf=2.0
    e=2.71
    !open(1,file='timeconstant.dat',status='unknown')
    do while (t<=tf)
        call didt(f,i,V,R,L)
        !print *,t,i
        !write(1,*)t,i
        prev_i=i
        i=i+h*f
        t=t+h
        new_i=i
        if (new_i-(V/R*(1-1/e))<=0.01) then
            time=t
        end if
        !if (new_i-prev_i<=0.0001) then
       ! exit
        !end if
    end do

close(1)
close(2)
print *,"Final value of current:",i
print *,"Timeconstant of the circuit:",time
end program timeconstant
subroutine didt(f,i,V,R,L)
    real::f,i,V,R,L
    V=20.0
    R=2.0
    L=0.4
    f=(V-R*i)/L
    return
end subroutine
