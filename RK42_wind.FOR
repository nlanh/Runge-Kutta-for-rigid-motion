c	Runge Kutta cho he cac phuong trinh vi phan cap mot
      PROGRAM RK2
      IMPLICIT none
      
c     N: so luong phuong trinh, nsteps: so diem, tstep: buoc chia
c     y(1): vi tri dau, y(2): van toc dau
      REAL*8 t, tstep, y(5), yi, yf, ywi, yw, ya
      INTEGER N,j,nsteps
      N = 2
      nsteps = 2000
      tstep = 0.01
      y(1) = 0.
      PRINT*, 'Nhap omega_0 (rad/s)'
      READ*, y(2)
      
c     Mo file
      OPEN(6, FILE = 'RK4.dat')
      WRITE (6,'(4F20.5)') 0.d0, y(1), y(2), 0.

c     Vong lap n buoc cua thuat toan Runge-Kutta
      
      yi = y(1)
      ywi = y(2)
      ya = 0.
      DO 60 j = 1, nsteps
       t=j*tstep
       CALL rk4(t, y, tstep, N)
       yf = y(1)
       yw = (yf - yi)/tstep
       yi = yf
       ya = (yw - ywi)/tstep
       ywi = yw
       WRITE (6,'(4F20.10)') t, y(1), yw, ya
   60 CONTINUE
      CLOSE(6)
      STOP
      END
c--------------------Ket thuc chuong trinh chinh------------------------

c     Chuong trinh con Runge-Kutta 
      SUBROUTINE rk4(t, y, tstep, N)
      IMPLICIT none

      REAL*8 DERIV, h, t, tstep, y(5) 
      REAL*8 k1(5), k2(5),k3(5), k4(5), temp1(5), temp2(5), temp3(5)
      INTEGER i, N
      h = tstep/2.0

      DO 10 i = 1,N
       k1(i) = tstep * DERIV(t, y, i)
       temp1(i) = y(i) + 0.5*k1(i)
   10 CONTINUE

      DO 20 i = 1,N
       k2(i) = tstep * DERIV(t+h, temp1, i)
       temp2(i) = y(i) + 0.5*k2(i)
   20 CONTINUE

      DO 30 i = 1,N
       k3(i) = tstep * DERIV(t+h, temp2, i)
       temp3(i) = y(i) + k3(i)
   30 CONTINUE

      DO 40 i = 1,N
       k4(i) = tstep * DERIV(t+tstep, temp3, i)
       y(i) = y(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
   40 CONTINUE
 
      RETURN
      END

c     Ham so dao ham
      FUNCTION DERIV(t,temp,i)
      IMPLICIT none
      
      REAL*8 DERIV, t, temp(5)
      INTEGER i
      REAL*8 M, g, L, D, IT
      REAL*8 w0
      
      M = 4./3.*3.14159*0.02**3*7000.
      g = 9.8
      L = 1.0
      D = 12*1.293*0.02**2
      IT = 2./5.*M*0.02**2 + M*L**2
      w0 = 20.
       
      IF (i .EQ. 1) DERIV = temp(2)
      IF (i .EQ. 2) DERIV = -M*g*L/IT*sin(temp(1))
     &                      -L*D/IT
     &      *sqrt(abs((temp(2)*L)**2+w0**2-2*temp(2)*L*w0*cos(temp(1))))
     &       *(temp(2)*L-w0*cos(temp(1)))

      RETURN
      END
