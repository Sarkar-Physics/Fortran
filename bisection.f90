program bisection
    real::f,c,x=0,y1,y2
    open(1,file='tan.dat',status='unknown')
    open(2,file='y_x.dat',status='unknown')
    open(3,file='root.dat',status='unknown')
    call root(c)
    print*,"The value of root using bisection method",c
   ! write(3,*)c,f(c)+c
    do while(x<=8.0)
        y1=tan(x)
        y2=x
        x=x+0.01
        !write(1,*)x,y1
       ! write(2,*)x,y2
    end do
    close(1)
    close(2)
    close(3)
end program
function f(q)
    real::q
    f=tan(q)-q
    return
end function
subroutine root(c)
    real::a,b,tol,c,x1,xh
    a=1.6
    b=4.6
    tol=0.001
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
