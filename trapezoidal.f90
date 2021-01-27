program trapizoidal
    real::a,b
    print*,"insert lower limit:"
    read *,a
    print*,"insert upper limit:"
    read*,b
    call trapizoid(a,b)
end program

real function f(x)
f=sqrt(1+x**2)
return
end function

subroutine trapizoid(a,b)
    real::integration,sum,h,p
    integer::n,i
    n=100
    h=(b-a)/n
    p=(h/2.0)*(f(a)+f(b))
    sum=0
    do i=1,n-1
        sum=sum+h*f(a+h*i)
    end do
    integration=p+sum
    print*,"The value of integration:",integration
end subroutine
