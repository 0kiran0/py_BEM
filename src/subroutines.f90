subroutine cubicspline(yout,xin,na,xcp,ycp,ncp)
! Creates a cubic B-spline using the control points.
! input: control points (xcp,ycp) and ncp is the no. of control points
! output: Spline points (xs,ys) and nspline is the number of splined points

implicit none

integer i,j,k,ncp1,nbsp,nspline,np,nx,nax,ia,na,ncp
!f2py intent(in) na, ncp

parameter(np=50,nx=500,nax=100)

real*8 xc(ncp+2), yc(ncp+2)
real*8 xbs(nx), ybs(nx), xs(nx), ys(nx)
real*8 t(nx), T1(nx), T2(nx), T3(nx), T4(nx)
real*8 y_spl_end(ncp),xin(na),yout(na),xcp(ncp), ycp(ncp)
!f2py intent(in) xcp, ycp, xin
!f2py intent(out) yout
!f2py depend(na) xin, yout
!f2py depend(ncp) xcp, ycp
real*8 d1_B11,B11,d1_B22,B22,d1_B33,B33,d1_B44,B44
real*8 ys_0,d1_ys_0,tt_0,tt

!t(1) = 0.
!t(np) = 1.0
do i=1,np
  t(i)  = (1.0/(np-1))*(i-1)
  T1(i) = ((-t(i)**3) + (3*t(i)**2) - (3*t(i)) + 1)/6
  T2(i) = ((3*t(i)**3) - (6*t(i)**2) + 4)/6
  T3(i) = ((-3*t(i)**3) + (3*t(i)**2) + (3*t(i)) + 1)/6
  T4(i) = (t(i)**3)/6
enddo
! Making the start and end points as collocation points----
! Start point: Considering the 1st CP as the mid point of the 2nd CP and the start point----------
xc(1) = 2*xcp(1) - xcp(2)
yc(1) = 2*ycp(1) - ycp(2)
! End point: Considering the last CP as the mid point of the end point and the penultimate CP ----
xc(ncp+2) = 2*xcp(ncp) - xcp(ncp-1)
yc(ncp+2) = 2*ycp(ncp) - ycp(ncp-1)
do i=1,ncp
  xc(i+1) = xcp(i)
  yc(i+1) = ycp(i)
enddo
! write(*,*)
! print*,'New control points with start and end points as collocation points:'
! do i=1,ncp+2
! print*,xc(i),yc(i) ! span,property
! enddo
! write(*,*)
ncp1 = ncp + 2 ! includes the start and end points
!print*,'ncp-new:',ncp1
nbsp = (np)*(ncp1 - 3)
!print*,'nbsp:',nbsp
!--------------------------------------------------------------------------------------------------
!constructing the b-spline curve(cubic b-spline) using 4 points P0(x0,y0),P1(x1,y1),P2(x2,y2),P3(x3,y3)------
!-----B(t) =(((1-t)^3)*P0) + (3*((1-t)^2)*t*P1) + (3*(1-t)*(t^2)*P2) + ((t^3)*P3); 0.le.t.le.1 ---------
do j = 1, ncp1-3
 do i = 1, np
  xs(i) = (xc(j)*T1(i)) + (xc(j+1)*T2(i)) + (xc(j+2)*T3(i)) + (xc(j+3)*T4(i))
  ys(i) = (yc(j)*T1(i)) + (yc(j+1)*T2(i)) + (yc(j+2)*T3(i)) + (yc(j+3)*T4(i))
 ! print*,xs(i),ys(i)
  k = i + (np)*(j - 1)! Transforming from 2D array to 1D array
  xbs(k) = xs(i)
  ybs(k) = ys(i)
 ! print*,xbs(k),ybs(k)
 enddo
 !picking up the endpoints of each spline for Newton interpolation
 if (j==1) then
  y_spl_end(1) = ys(1)
 endif		   
 y_spl_end(j+1) = ys(np)
 !00000000000000000000000000000000000000000000000000000
 ! print*,y_spl_end(j)
enddo
! print*,xbs(nbsp+2),ybs(nbsp+2)
! write(*,*)
! do i = 1, nbsp
	! print*,i,ybs(i),xbs(i)
! enddo 
! write(*,*)
!
nspline = nbsp!+2
! Newton method to find yout corresponding to xin
! for 1st control point
yout(1) = xbs(1)
				   
! for last control point
yout(na) = xbs(nspline) !xbs(np*(ncp-3))
! print*,
! print*, y_spl_end
! print*,
do j = 1, ncp1-3 ! including the phantom CPs
  do i = 2, na-1
	if ((xin(i) .ge. y_spl_end(j)).and.(xin(i) .lt. y_spl_end(j+1)))then
       ! print*,"i,xin(i),y_spl_end(j):",i,xin(i),y_spl_end(j)
	   tt_0 = 0.3
	   do k = 1, 10
	   ! Basic functions:
		B11 = ((-tt_0**3)   + (3*tt_0**2) - (3*tt_0) + 1)/6
		B22 = ((3*tt_0**3)  - (6*tt_0**2) + 4)/6
		B33 = ((-3*tt_0**3) + (3*tt_0**2) + (3*tt_0) + 1)/6
		B44 = (tt_0**3)/6			   
	   ! first derivative:
		d1_B11 = ((-3*tt_0**2) + (6*tt_0) - (3))/6
		d1_B22 = ((9*tt_0**2)  - (12*tt_0))/6
		d1_B33 = ((-9*tt_0**2) + (6*tt_0) + (3))/6
		d1_B44 = (3*tt_0**2)/6
		! xs(t_0)
		  ys_0  = yc(j)*B11    + yc(j+1)*B22    + yc(j+2)*B33    + yc(j+3)*B44
		! d1_xs(t_0)
		d1_ys_0 = yc(j)*d1_B11 + yc(j+1)*d1_B22 + yc(j+2)*d1_B33 + yc(j+3)*d1_B44
		! Newton's interpolation:
		 tt = tt_0 + (xin(i)-ys_0)/d1_ys_0
		 ! print*, tt
		 if (abs(tt-tt_0)<1e-16) then
		   B11 = ((-tt**3)   + (3*tt**2) - (3*tt) + 1)/6
		   B22 = ((3*tt**3)  - (6*tt**2) + 4)/6
		   B33 = ((-3*tt**3) + (3*tt**2) + (3*tt) + 1)/6
		   B44 = (tt**3)/6			
			goto 20
		 endif
		 tt_0 = tt
		enddo
20 	    yout(i) = xc(j)*B11    + xc(j+1)*B22    + xc(j+2)*B33    + xc(j+3)*B44
     ! else
        ! print*,"xin not in the spline segment!"
	 endif			 
	 ! print*,'yout',yout(i)
 enddo
enddo
! do i = 1, na
    ! print*,'xin,yout:',xin(i),yout(i)
! enddo      
return
end subroutine cubicspline
!***************************************************************
!***************************************************************

!***************************************************************
subroutine cubicbspline_intersec(y_spl_end,xcp,ycp,ncp,xin,yout,na,xbs,ybs)
! Curve line intersection procedure using Newton-Raphson method....Kiran Siddappaji and Ahmed Nemnem 4/17/13
!Curve: F(y) (or F(x)). Line y = c or x = c.
! The root(s) of the F(y) - y_line (or F(x) - x_line) is the intersection point.
!inputs:
! (xcp,ycp): control points 
! ncp: no. of control points
! (xbs,ybs): spline points
! nbs: no. of spline points
!xin: x (or y) coordinate of the intersection point.
!outputs:
!yout: y (or x) coordinate of the intersection point.

implicit none

integer i,j,k,ia,na,ncp,nbs,nx,np

parameter(np=50,nx = 500)

real*8 xcp(ncp),ycp(ncp),xinarray(na),youtarray(na)
real*8 xin(na),yout(na), xp2,yp2,xp1,yp1, xvar
real*8 xbs(np*(ncp-3)),ybs(np*(ncp-3)) ! spline coordinates
real*8 y_spl_end(ncp-2)
real*8 d1_B11,B11,d1_B22,B22,d1_B33,B33,d1_B44,B44
real*8 ys_0,d1_ys_0,tt_0,tt
 
 !print*,'y_spl_end',y_spl_end

! Newton method to find yout corresponding to xin
    ! for 1st control point
		       yout(1) = xbs(1)
               			   
    ! for last control point
			   yout(na) = xbs(np*(ncp-3))

      do j=1,ncp-3
          do i=2,na-1
            if ((xin(i) > y_spl_end(j)).and.(xin(i) < y_spl_end(j+1)))then
			   tt_0 = 0.3
			   do k =1,10
			   ! Basic functions:
				B11 = ((-tt_0**3) + (3*tt_0**2) - (3*tt_0) + 1)/6
				B22 = ((3*tt_0**3) - (6*tt_0**2) + 4)/6
				B33 = ((-3*tt_0**3) + (3*tt_0**2) + (3*tt_0) + 1)/6
				B44 = (tt_0**3)/6			   
			   ! first derivative:
				d1_B11 = ((-3*tt_0**2) + (6*tt_0) - (3))/6
				d1_B22 = ((9*tt_0**2) - (12*tt_0))/6
				d1_B33 = ((-9*tt_0**2) + (6*tt_0) + (3))/6
				d1_B44 = (3*tt_0**2)/6
				! xs(t_0)
				  ys_0 = ycp(j)*B11+ycp(j+1)*B22+ycp(j+2)*B33+ycp(j+3)*B44
				! d1_xs(t_0)
				d1_ys_0= ycp(j)*d1_B11+ycp(j+1)*d1_B22+ycp(j+2)*d1_B33+ycp(j+3)*d1_B44
				! Newton's interpolation:
			     tt = tt_0 + (xin(i)-ys_0)/d1_ys_0
	             if (abs(tt-tt_0)<1e-16) then
				   B11 = ((-tt**3) + (3*tt**2) - (3*tt) + 1)/6
				   B22 = ((3*tt**3) - (6*tt**2) + 4)/6
				   B33 = ((-3*tt**3) + (3*tt**2) + (3*tt) + 1)/6
				   B44 = (tt**3)/6			
				    goto 20
	             endif
			     tt_0 = tt
                enddo
20     		   yout(i)=xcp(j)*B11+xcp(j+1)*B22+xcp(j+2)*B33+xcp(j+3)*B44
             endif			 
!			 print*,'yout',yout(i)
	     enddo
      enddo
 !     do i = 1, na
!	    print*,'xin',xin
!		print*,'yout',yout
  !	  enddo	
return
end subroutine cubicbspline_intersec
!***************************************************************
!***************************************************************
!**********************************************************************
subroutine curv_line_inters(yout,xin,nspan,xbs,ybs,nspline)
! Calculates the intersection point between the line and the curve
! input: curve points (xbs,ybs), no. of curve points (nspline),xin from the line
!output: yout for the line on the curve.

implicit none

integer ii,j,k, nspline, nx,nspan
parameter(nx=1000)
real*8 xbs(nspline),ybs(nspline), xin(nspan), yout(nspan)
real*8 min, max, xmax, xmin, xint
!print*,'xin:',xin
!f2py intent(in) xin,xbs,ybs
!f2py intent(out) yout
!f2py depend(nspan) yout, xin
!f2py depend(nspline) xbs,ybs
do j = 1, nspan
 do ii = 1, nspline-1
  xint = xbs(ii+1) + (ybs(ii+1) - xin(j))*(((xbs(ii) - xbs(ii+1))/(ybs(ii+1) - ybs(ii))))
!  print*,'xint',xint
  xmax = max(xbs(ii),xbs(ii+1))
  xmin = min(xbs(ii),xbs(ii+1))
!  print*,xmin,xmax
!  xmin1 = xmin(i)
!  xmax1 = xmax(i)
  if (xint.ge.xmin.and.xint.le.xmax)then
     yout(j) = xint
  elseif(xint.eq.0.)then
     print*," Curve-Line Intersection error."
     print*," Intersection point was not found"
     print*,"xint: ",xint
     STOP
  endif
 enddo
!write(*,*)
! Forcing Start and End points to be the same as control points at start and end.
! if(j.eq.1)then
 ! yout(j) = xbs(1)
! elseif(j.eq.nspan)then
 ! yout(j) = xbs(nspline)
! endif 
! Forcing only the end point to be the same as the point on the curve
!if(j.eq.nspan) yout(j) = xbs(nspline)
!print*,xmin,xmax
! print*,yout(j)
enddo
return
end subroutine curv_line_inters
!***************************************************************
!**********************************************************************