!program to calculate moment of inertia of a eliptical lamina
program MI
real::a=14,b=8.0,c,l=-8.0,m=1,d=1
call trapizoid(l,b,c)
print*,"Moment of inertia of the object about y axis:",(2*m*a*a*c)/(3*3.1415*b)
print*,"Moment of inertia of the object about y' axis:",m*d**2+(2*m*a*a*c)/(3*3.1415*b)
end program

real function f(x)
f=(1-x**2/8**2)**1.5
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
