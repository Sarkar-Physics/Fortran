program bisection
    real::f,c,r
    call root(c)
    print*,"The value of root using bisection method",c
    r=sqrt(2.123)
    print *,"The value using sqrt function",r
end program
function f(q)
    real::q
    f=q**2-2.123
    return
end function
subroutine root(c)
    real::a,b,tol,c,x1,xh
    a=1.0
    b=2.0
    tol=0.001
    x1=a
    xh=b
    do while(abs(x1-xh)>=tol)
        c=(x1+xh)/2
        if (abs(f(c))<=tol) then
            exit
        elseif(f(c)*f(x1)<0) then
            xh=c
        else
            x1=c
       end if
    end do
    return
end subroutine
