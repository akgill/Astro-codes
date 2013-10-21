

; main function for non-linear convection term project
pro main, numz, numt, n_rayleigh, n_prandtl, c_param, a_param, bc_type=bc_type

; define variables in common block
	common share, n_r, n_pr, a, c, nz, nt, dz, dt, z, time, $
		alpha, t_mean, phi, theta, w, bc

; define static variables
	if ~keyword_set(bc_type) then bc='rigid' else bc=bc_type
	n_r = n_rayleigh		; rayleigh number
	n_pr = n_prandtl		; prandtl number
	a = a_param			; horizontal wavenumber
	c = c_param			; planform configuration, c=0 or c=1./sqrt(6)
	nz = numz+1			; number of z grid points
	nt = numt+1			; number of timesteps
	dz = 1./numz			; size of z mesh step

;	dt = 1d-5
	dt = 0.25*dz^2/n_pr		; define time step size such that 
					; n_pr*dt/dz^2 = 0.25 for stability criterion
					; in FTCS numerical integration

	z = findgen(nz)*dz		; array of z grid points
	time = findgen(nt)*dt		; an array of timesteps

	k = 2.+a^2*dz^2			; a factor that appears in the approximation of w
					; from phi

	alpha = fltarr(nz-2)		; an array of constants that is needed to
					; calculate w's from phi's

	alpha[0] = 1 & alpha[1] = k	; alpha's that are not defined recursively
	for i=2,nz-3,1 do begin		; alpha's that are defined recursively
		alpha[i] = k*alpha[i-1] - 1
	endfor
	
; define evolving variables
	t_mean = fltarr(nz,3)		; an array of t_mean at all grid points at all times
	phi = fltarr(nz,3)		; an array of phi at all grid points at all times
	theta = fltarr(nz,3)		; an array of theta at all grid points at all times
	w = fltarr(nz,3)		; an array of w at all grid points at all times

; define initial conditions
;	t_mean[0,0] = 1.		; switch on heat plate at t=0, 
;	t_mean[1:*, 0] = 0.		; for z>0, t_mean=0 at t=0
	t_mean[*,0] = 1.-z			; t_mean starts as a linear gradient with z
;	w[*,0] = randomn(seed1, nz)*dz/a	; random, small amplitude velocity perturbations
;	phi[*,0] = randomn(seed2, nz)*dz/a	; random, small amplitude phi perturbations
;	theta[*,0] = randomn(seed3, nz)*dz	; initially quiescent, no temperature perturbations
	if (bc eq 'rigid') then w[*,0] = (dz/a)*(sin(!pi * z))^2 else w[*,0] = (dz/a)*(sin(!pi * z))
	phi[1:nz-2,0] = (w[0:nz-3,0] - 2*w[1:nz-2,0] + w[2:nz-1,0])/dz^2 - a^2*w[1:nz-2,0]
	theta[*,0] = (dz/a)*sin(!pi * z)

	impose_boundary,t=0

	phi[*,2] = phi[*,0]
	theta[*,2] = theta[*,0]
	t_mean[*,2] = t_mean[*,0]
	w[*,2]=w[*,0]

	for i=1,nt-1,1 do begin
		update_vars
		impose_boundary
		update_w

		phi[*,0] = phi[*,1]
		theta[*,0] = theta[*,1]
		t_mean[*,0] = t_mean[*,1]
		w[*,0]=w[*,1]

		plot,t_mean[*,1]/max(t_mean[*,1]),linestyle=0
		oplot,theta[*,1]/max(theta[*,1]),psym=4
		oplot,w[*,1]/max(w[*,1]),psym=3
	endfor	

	stop

;	!P.MULTI[0,2,2,0,1]

end



; procedure that enforces the boundary conditions at timestep t
pro impose_boundary,t=t
	common share

	if ~keyword_set(t) then t=1

	if (bc ne 'free') then begin
	
	; phi for rigid boundary condition
		phi[0,t] = (2./dz^2)*(w[1,t]) 		; phi simplified by requirement that
		phi[nz-1,t] = (2./dz^2)*(w[nz-2,t]) 	; DW = 0 at boundaries, implying that
							; w[0,*] = w[1,*] and 
							; w=[nz-1,*] = w[nz-2,*]
	endif else begin

	; phi for free boundary condition
		phi[0,t] = 0.		; D^2W = 0 at z=0
		phi[nz-1,t] = 0.	; D^2W = 0 at z=1
	endelse

; general boundary conditions
	theta[0,t] = 0.			; temperature fluctuations vanish
	theta[nz-1,t] = 0.		; at boundaries
	w[0,t] = 0.			; velocity fluctuations vanish
	w[nz-1,t] = 0.			; at boundaries
	t_mean[0,t] = 1.		; t_mean = 1 at z=0
	t_mean[nz-1,t] = 0.		; t_mean = 0 at z=1

end



; procedure that updates phi, theta, and t_mean at the timestep t
pro update_vars
	common share

	for i=1, nz-2, 1 do begin
		phi1 = (phi[i-1,0] - 2*phi[i,0] + phi[i+1,0])/dz^2 	; d^2/dz^2 phi
		phi2 = -a^2*(phi[i,0] - n_r*theta[i,0]) 			; phi a^2 terms

		
		theta1 =  (theta[i-1,0] - 2*theta[i,0] + theta[i+1,0])/dz^2; d^2/dz^2 theta
		theta2 = -a^2*theta[i,0] - $					; theta a^2 and d/dz T terms
			w[i,0]*(t_mean[i+1,0] - t_mean[i-1,0])/(2.*dz)

		temp1 = (t_mean[i-1,0] - 2*t_mean[i,0] + t_mean[i+1,0])/dz^2; d^2/dz^2 t_mean
		temp2 = -w[i,0]*(theta[i+1,0] - theta[i-1,0])/(2.*dz) 	; w*(d/dz theta)
		temp3 = -theta[i,0]*(w[i+1,0] - w[i-1,0])/(2.*dz) 		; theta*(d/dz w)


		if (c eq 0) then begin 				; planform dependent term
			phi3 = 0 
			theta3 = 0
		endif else begin
			phi3 = -(c/n_pr)*((2.*phi[i,0]*(w[i+1,0] - w[i-1,0])) + $
				w[i,0]*(phi[i+1,0] - phi[i-1,0]))/(2.*dz)

			theta3 = -(c)*((theta[i,0]*(w[i+1,0] - w[i-1,0])) + $
				2.*w[i,0]*(theta[i+1,0] - theta[i-1,0]))/(2.*dz)
		endelse

		phi[i,1] = (phi1+phi2+phi3)*n_pr*dt + phi[i,0]		; update phi
		theta[i,1] = (theta1+theta2+theta3)*dt + theta[i,0]; update theta
		t_mean[i,1] = (temp1+temp2+temp3)*dt + t_mean[i,0]	; update t_mean
	;	print,'i=',i
	;	print,'temp1: ',temp1
	;	print,'temp2: ',temp2
	;	print,'temp3: ',temp3
	;	print,'phi1: ',phi1
	;	print,'phi1: ',phi2
	;	print,'phi1: ',phi3
	;	print,'theta1: ',theta1
	;	print,'theta1: ',theta2
	;	print,'theta1: ',theta3
	endfor
end



; procedure that updates w's at the timestep t using the updated phi's
pro update_w
	common share
	
	w[nz-1,1] = 0. 		; boundary condition specifies
	w[0,1] = 0.		; w=0 at z=0,1

	sum_terms = alpha[0:nz-3]*phi[1:nz-2,1]*dz^2

	for i=nz-2,1,-1 do begin
		sum = total(sum_terms[0:i-1])
		w[i,1] = (w[0,1]-sum+w[i+1,1])/alpha[i-2]
		;print,w[i,1]

	endfor

end


