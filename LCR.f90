program LCR
real::R,x=0
print *,"Enter initial value:"
read*,R
call new_rap(R)
open(1,file="LCR.dat",status="unknown")
do while(x<=440)
x=x+0.1
write(1,*)x,f(x)
end do
end program
function f(x)
real::x,L=5.0,C=0.0001,t=0.05,a,b
a=x/(2*L)
b=sqrt(1/(L*C)-a**2)
f=exp(-a*t)*cos(b*t)-0.01
return
end function

function g(x)
real::x,L=5.0,C=0.0001,t=0.05,a,b,d
a=x/(2*L)
b=sqrt(1/(L*C)-a**2)
d=1/(2*L)**2
g=exp(-a*t)*(-t/(2*L)*cos(b*t)+x*(t*d/b)*sin(b*t))
return
end function

subroutine new_rap(x)
real::h,tol,x
integer::iteration
tol=0.0001
iteration=0
h=-f(x)/g(x)
do while(abs(h)>tol)
if(iteration>100) then
print *,"root finding not posssible"
exit
end if
if(g(x)==0.0) then
print *,"root Finding not possible"
exit
end if
h=-f(x)/g(x)
x=x+h
iteration=iteration+1
!print*,iteration,h,x
end do
print*,"The required reisitance(in ohm):",x
print *, "No of Iteration is",iteration
return
end subroutine
