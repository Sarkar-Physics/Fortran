program pion_multiplicity
    integer::x(119),y(119),z(119),count,a,b=-1
    open(1,file='multso2_1',status='unknown')
    open(2,file='pion_multiplicity_frequency.dat',status='unknown')
! reading the data from the given file
    do i=1,119
    read(1,*)x(i),y(i)
    z(i)=0
    end do
! sorting the array in ascending order
   do i=1,119
    do j=i+1,119
        if(y(i)>y(j)) then
        a=y(i)
        y(i)=y(j)
        y(j)=a
        end if
    end do
   end do
! counting the frquency of data
    do i=1,119
        count=0
        do j=1,119
            if(y(i)==y(j)) then
            count=count+1
            z(j)=b
            end if
        end do
    if(z(j)/=b) then
    z(i)=count
    end if
    end do
! writing the data into the output file
    do i=1,119
        if(z(i)/=b) then
            write(2,*)y(i),z(i)
        end if
    end do
    close(1)
    close(2)
end program pion_multiplicity
