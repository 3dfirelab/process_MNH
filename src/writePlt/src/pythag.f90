

      real FUNCTION PYTHAG (A, B)

!***PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.

      real,intent(in)  :: A,B
      real             :: P,Q,R,S,T

      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         if(p.gt.1.e10) then
         r=0.
         else
         R = (Q/P)**2
         end if
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
