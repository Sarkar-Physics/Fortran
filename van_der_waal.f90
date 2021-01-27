program van_der_waal_gas
real*8::v,f,g,x
v=0.5
print *,"Enter initial value:"
read*,x
open(1,file="gas.dat",status="unknown")
do while(v<=1.5)
v=v+0.001
write(1,*)v,f(v)
!print*,v,f(v)
end do
call new_rap(x)
end program

real*8 function f(v)
real*8::v,p=1.0,a=0.00874,b=0.0023,R=0.00368652,T=273
f=p*v**3-(p*b+R*T)*v**2+a*v-a*b
return
end function

real*8 function g(v)
real*8::v,p=1.0,a=0.00874,b=0.0023,R=0.00368652,T=273
g=3*p*v**2-2*(p*b+R*T)*v+a
return
end function

subroutine new_rap(x)
real*8::x,h,tol,f,g
integer::iteration
tol=1e-8
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
print*,"Required Molar volume:"
print"(f10.8)",x
print *, "No of Iteration is",iteration
return
end subroutine
