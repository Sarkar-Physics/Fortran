program projectile
real::x,y,iteration,z,p=0.0
c=10.0
call new_rap(c)
z=r(c)
print*,"The range of the projectile is(in m):",z
open(1,file="projectile.dat",status="unknown")
do i=1,1000
x=r(p)
y=f(p)
p=p+0.1
write(1,*)x,y
end do

end program

function f(x)
real::x
f=1600*(1-exp(-x/5))-160*x
return
end function

function g(x)
real::x
g=160*(2/5*exp(-x/5)-1)
return
end function

function r(x)
real::x
r=800*(1-exp(-x/5))
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
print*,"The object hits the ground at(in sec):"
print "(f10.3)",x
print *, "No of Iteration is",iteration
return
end subroutine
