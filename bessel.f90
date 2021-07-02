program bessel_funcions
    integer::i
    real*8::x
    open(1,file='bessel.dat',status='unknown')
    do i=0,5
        call bessel(i,x)
    end do
    close(1)
end program

real*8 function fact(n)
integer::n
fact=1.0
do i=2,n
    fact=fact*i
end do
return
end function


subroutine bessel(n,bes)
    integer::n,s
    real*8::p,q,x,bes,fact
    x=0.0
    do while(x<=20.0)
    bes=0.0
    do s=0,50
        bes=bes+((((-1.0)**s)/(fact(n+s)*fact(s))))*(x/2.)**(n+2.0*s)
        end do
        write(1,*)x,bes!,bessel_jn(0,x)
        x=x+0.001
    end do
    write(1,*)
    write(1,*)
    return
end subroutine
