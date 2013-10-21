

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
	t_mean = fltarr(nz,nt)		; an array of t_mean at all grid points at all times
	phi = fltarr(nz,nt)		; an array of phi at all grid points at all times
	theta = fltarr(nz,nt)		; an array of theta at all grid points at all times
	w = fltarr(nz,nt)		; an array of w at all grid points at all times

; define initial conditions
;	t_mean[0,0] = 1.		; switch on heat plate at t=0, 
;	t_mean[1:*, 0] = 0.		; for z>0, t_mean=0 at t=0
	t_mean[*,0] = 1.-z			; t_mean starts as a linear gradient with z
;	w[*,0] = randomn(seed1, nz)*dz/a	; random, small amplitude velocity perturbations
;	phi[*,0] = randomn(seed2, nz)*dz/a	; random, small amplitude phi perturbations
;	theta[*,0] = randomn(seed3, nz)*dz	; initially quiescent, no temperature perturbations
	w[*,0] = (dz/a)*(sin(!pi * z))^2
	phi[1:nz-2,0] = (w[0:nz-3,0] - 2*w[1:nz-2,0] + w[2:nz-1,0])/dz^2 - a^2*w[1:nz-2]
	theta[*,0] = (dz/a)*sin(!pi * z)

	impose_boundary, 0

	for i=1,nt-1,1 do begin
		update_vars, i
		update_w, i
		impose_boundary, i
	endfor	

	stop

;	!P.MULTI[0,2,2,0,1]

end



; procedure that enforces the boundary conditions at timestep t
pro impose_boundary, t
	common share

	if (bc ne 'free') then begin
	
	; phi for rigid boundary condition
		phi[0,t] = (2./dz^2)*(w[1,t]-w[0,t]) 		; phi simplified by requirement that
		phi[nz-1,t] = (2./dz^2)*(w[nz-2,t]-w[nz-1,t]) 	; DW = 0 at boundaries, implying that
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
pro update_vars, t
	common share

	for i=1, nz-2, 1 do begin
		phi1 = (phi[i-1,t-1] - 2*phi[i,t-1] + phi[i+1,t-1])/dz^2 	; d^2/dz^2 phi
		phi2 = -a^2*(phi[i,t-1] - n_r*theta[i,t-1]) 			; phi a^2 terms

		
		theta1 =  (theta[i-1,t-1] - 2*theta[i,t-1] + theta[i+1,t-1])/dz^2; d^2/dz^2 theta
		theta2 = -a^2*theta[i,t-1] - $					; theta a^2 and d/dz T terms
			w[i,t-1]*(t_mean[i+1,t-1] - t_mean[i-1,t-1])/(2.*dz)

		temp1 = (t_mean[i-1,t-1] - 2*t_mean[i,t-1] + t_mean[i+1,t-1])/dz^2; d^2/dz^2 t_mean
		temp2 = -w[i,t-1]*(theta[i+1,t-1] - theta[i-1,t-1])/(2.*dz) 	; w*(d/dz theta)
		temp3 = -theta[i,t-1]*(w[i+1,t-1] - w[i-1,t-1])/(2.*dz) 		; theta*(d/dz w)


		if (c eq 0) then begin 				; planform dependent term
			phi3 = 0 
			theta3 = 0
		endif else begin
			phi3 = -(c/n_pr)*(2.*phi[i,t-1]*(w[i+1,t-1] - w[i-1,t-1])) + $
				w[i,t-1]*(phi[i+1,t-1] - phi[i-1,t-1])/(2.*dz)

			theta3 = -(c/n_pr)*(2.*phi[i,t-1]*(w[i+1,t-1] - w[i-1,t-1])) + $
				w[i,t-1]*(phi[i+1,t-1] - phi[i-1,t-1])/(2.*dz)
		endelse

		phi[i,t] = (phi1+phi2+phi3)*n_pr*dt + phi[i,t-1]		; update phi
		theta[i,t] = (theta1+theta2+theta3)*n_pr*dt + theta[i,t-1]	; update theta
		t_mean[i,t] = (temp1+temp2+temp3)*n_pr*dt + t_mean[i,t-1]	; update t_mean
	endfor
end



; procedure that updates w's at the timestep t using the updated phi's
pro update_w, t
	common share
	
	w[nz-1,t] = 0.
	w[0,t] = 0.

	sum_terms = alpha[0:nz-3]*phi[1:nz-2,t]*dz^2

	for i=nz-2,1,-1 do begin
		sum = total(sum_terms[0:i-1])
		w[i,t] = (w[0,t]-sum+w[i+1,t])/alpha[i-2]

	endfor
end


