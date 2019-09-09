subroutine h1e(dm,bas2atm,hcore,natm,nao,energy_h1e)
implicit none

integer i,j,l
integer a,b
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(natm,nao,nao)
real*8 energy_h1e(natm)
real*8 tem

!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e

! -----------------------------------------------------------------------
! f2py -m h1e_new -c h1e_new.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

energy_h1e = 0.0d0
call system_clock(time1)

do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
    do l=1,natm
      tem = dm(i,j) * hcore(l,j,i)
      energy_h1e(l) = energy_h1e(l) + tem*1/2.0d0
      energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
      energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
    enddo
  enddo
enddo
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e

subroutine bgh1e(dm,bas2atm,hcore,natm,nao,energy_h1e)
implicit none

integer i,j,l
integer a,b
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(nao,nao)
real*8 energy_h1e(natm)
real*8 tem

!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e

! -----------------------------------------------------------------------
! f2py -m h1e_new -c h1e_new.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

energy_h1e = 0.0d0
call system_clock(time1)

do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
      tem = dm(i,j) * hcore(j,i)
      energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
      energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
  enddo
enddo
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine bgh1e

subroutine h1e_inter(dm,bas2atm,hcore,natm,nao,energy_h1e, e1_1, e1_2, e1_3)
implicit none

integer i,j,l
integer a,b
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(natm,nao,nao)
real*8 energy_h1e(natm)
real*8 e1_1(natm)
real*8 e1_2(natm,natm)
real*8 e1_3(natm,natm,natm)
real*8 tem

call system_clock(time1)
!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e,e1_1, e1_2, e1_3

energy_h1e = 0.0d0
e1_1=0.0d0 
e1_2=0.0d0
e1_3=0.0d0

!$OMP PARALLEL DO schedule(guided) &
!$omp default(private) &
!$omp shared(dm,bas2atm,hcore,nao,natm) &
!$omp reduction(+:energy_h1e,e1_1,e1_2,e1_3)
do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
    do l=1,natm
      tem = dm(i,j) * hcore(l,j,i)
      energy_h1e(l) = energy_h1e(l) + tem*1/2.0d0
      energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
      energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
      if (a==b) then
        if (a==l) then
          e1_1(a) = tem
        else
          e1_2(a,l) = tem
        endif
      else 
        if ((a==l) .or. (b==l)) then
          e1_2(a,b) = tem
        else
          e1_3(a,b,l) = tem
        endif
      endif
    enddo
  enddo
enddo
!$omp end parallel do
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e_inter
