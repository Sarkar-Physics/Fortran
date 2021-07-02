program legendre_polynomials
    integer::i
    real::x
    open(1,file='legendre.dat',status='unknown')
    do i=1,6
        call legendre(i,x)
    end do
    close(1)
end program

real function fact(n)
integer::n
fact=1.0
do i=2,n
    fact=fact*i
end do
return
end function

subroutine legendre(n,leg)
    integer::n,k,u
    real::p,q,r,s,x,leg,fact
    x=-1.0
    if(modulo(n,2)==0) then
        u=n/2
    else
        u=(n-1)/2
    end if
    do while(x<=1.0)
    leg=0.0
        do k=0,u
        leg=leg+((((-1.0)**k)*fact(2*n-2*k))/((2.0**n)*fact(k)*fact(n-k)*fact(n-2*k)))*x**(n-2.0*k)
        end do
    write(1,*)x,leg
    x=x+0.01
    end do
    write(1,*)
    write(1,*)
    return
end subroutine

