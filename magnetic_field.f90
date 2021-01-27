!program to show the variation of magnetic field
program magnetic_field
real::B_z,a=0,b=2*3.14,z=-5.0
open(1,file='magnetic_field.dat',status='unknown')
do while(z<=5.0)
call trapizoid(a,b,z,B_z)
write(1,*)z,B_z
z=z+0.1
end do
close(1)
end program

real function f(z,x)
real::x,z
f=0.25/(4*3.14*(0.25+z**2)**1.5)*x**0
return
end function

subroutine trapizoid(a,b,z,integration)
real::integration,sum,h,p
integer::n,i
n=100
h=(b-a)/n
p=(h/2.0)*(f(z,a)+f(z,b))
sum=0
do i=1,n-1
sum=sum+h*f(z,a+h*i)
end do
integration=p+sum
return
end subroutine
