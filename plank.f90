program plank
    real::h,c,k,T,x=0,y=0,f,r,lambda
    h=6.626e-34
    c=3e8
    k=1.38e-23
    T=5500
    open(1,file='blackbody_T=6000K.dat',status='unknown')
    open(2,file='lambda.dat',status='unknown')
    call root(r)
    print*,"The value of root using bisection method",r
    lambda=h*c/(r*k*T)
    !write(2,*)lambda,(2*h*c**2/lambda**5)*(1/(exp(h*c/(k*T*lambda))-1))
    print *, "Maximum wavelength:",lambda
    do while(x<=4e-6)
        y=(2*h*c**2/x**5)*(1/(exp(h*c/(k*T*x))-1))
        x=x+0.01e-6
        write(1,*)x,y
    end do
    close(1)
    close(2)
end program
function f(q)
    real::q
    f=q-5*(1-exp(-q))
    return
end function
subroutine root(c)
    real::a,b,tol,c,x1,xh
    a=4.5
    b=5.2
    tol=0.01
    x1=a
    xh=b
    do while(abs(x1-xh)>=tol)
        c=(x1+xh)/2
        if(f(c)*f(x1)<0) then
            xh=c
        else
            x1=c
       end if
    end do
    return
end subroutine
