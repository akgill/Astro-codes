;4th order Runge-Kutta integration of a closed, precessing,
;elliptical orbit
;
;Amandeep Gill
;*************************************************************
function func,y,vc,r_peri,a

;
;Lorenz model rhs's
;
   r_apo=3.73d*r_peri
   f = dblarr(2)
   E = vc^2*(alog(2d)+alog(3.73d)*((3.73d)^2)/((3.73d)^2-1))
   L2 = 2d*vc^2*alog(r_apo/r_peri)*r_apo^2*r_peri^2/(r_apo^2-r_peri^2)
   r = y[0]
   phi = y[1]
   V = (vc^2)*alog(r)
;

   if (a) then f[0] = -sqrt(2*(E-V) - L2/r^2) else f[0] = sqrt(2*(E-V) - L2/r^2)
   f[1] = sqrt(L2)/r^2
;
   return,f
end
;*************************************************************
pro rkstep,yold,ynew,step,vc,r_peri,a
;
   step2 = step/2d
   step6 = step/6d
;
   k1 = func(yold,vc,r_peri,a)
   k2 = func(yold+step2*k1,vc,r_peri,a)
   k3 = func(yold+step2*k2,vc,r_peri,a)
   k4 = func(yold+step*k3,vc,r_peri,a)
;
   ynew = yold + step6*(k1+2*k2+2*k3+k4)
;
end
;*************************************************************
function drive
;
;drive for R-K soln of a rosette orbit
;
   print,' Enter vc,r_peri,dt,nsteps'
   read,vc,r_peri,dt,nsteps
;
   yold = dblarr(2)
   ynew = dblarr(2)
;
   print,' Enter initial conditions r(0),phi(0),approaching'
   read,a,b,approaching
   yold = [a,b]

   data=dblarr(2,nsteps)
   for i=0,nsteps-1 do begin
       rkstep,yold,ynew,dt,vc,r_peri,approaching
       data[*,i]=ynew[*]
       if (approaching) then begin
           if (ynew[0] le (0.01+r_peri)) then approaching=0
       endif else begin
           if (ynew[0] ge (3.73d*r_peri-0.01)) then approaching=1
       endelse
       yold=ynew
   endfor
   sigr2=(total(((data[0,1:nsteps-1]-data[0,0:nsteps-2])/dt)^2)/(nsteps-1))
   sigt2=(total(((data[1,1:nsteps-1]-data[1,0:nsteps-2])*data[0,1:nsteps-1]/dt)^2)/(nsteps-1))
   print,'sigma_r^2: ',sigr2
   print,'sigma_t^2: ',sigt2
   return,data
end
;**********************************************************
pro rkplot,data,sigr2,sigt2

   r=data[0,*]
   phi=data[1,*]

   set_plot,'ps'
   device,file='rosette.ps',ysize=25,xsize=17,xoff=1,yoff=1

   !p.multi=0
   str="Rosette Orbit,       sigma_r^2: "+strtrim(string(sigr2))+"       sigma_t^2/2 :"+strtrim(string(sigt2/2))
	print,str
   plot,/polar,r,phi,title=str
   device,/close

   set_plot,'x'

   plot,/polar,r,phi,title=str
   

end
