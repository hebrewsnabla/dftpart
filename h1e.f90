subroutine h1e(dm,bas2atm,hcore,natm,nao,energy_h1e)
implicit none

integer i,j,l
integer a,b
integer ii,jj,kk,m
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(natm,nao,nao)
real*8 energy_h1e(natm+1)
real*8 tem

!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e

! -----------------------------------------------------------------------
! f2py -m h1e_old -c h1e_old.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

energy_h1e = 0.0d0
call system_clock(time1)

do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
    do l=1,natm
      tem = dm(i,j) * hcore(l,j,i)
      ii = natm + 1
      jj = natm + 1 
      kk = natm + 1 
      m = 0
      if (a/=b) then 
          ii = a 
          jj = b
          m = m + 2 
      else 
          ii = a 
          m = m + 1 
      endif 

      if (l/=ii) then 
          if (l/=jj) then 
              kk = l 
              m = m + 1 
          endif
      endif
      tem = tem/m 

      energy_h1e(ii) = energy_h1e(ii) + tem
      energy_h1e(jj) = energy_h1e(jj) + tem
      energy_h1e(kk) = energy_h1e(kk) + tem
    enddo
  enddo
enddo
call system_clock(time2)
write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e

