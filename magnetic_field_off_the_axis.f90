program magnetic_field
real::B_z,a=0,b=2*3.14
call trapizoid(a,b,B_z)
print*,B_z
close(1)
end program

real function f(x)
real::x,z=1.0,y=1.0
f=(0.5*z*sin(x))/(4*3.14*(0.25+y**2+z**2-2*y*0.5*sin(x))**1.5)
return
end function

subroutine trapizoid(a,b,integration)
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
return
end subroutine
