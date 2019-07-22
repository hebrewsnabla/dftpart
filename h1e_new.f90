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
write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e

