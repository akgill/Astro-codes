FUNCTION convec, numz, numt, n_rayleigh, n_prandtl, c_param, a_param, bc_type=bc_type

	;assigning constants
	nz = numz+1      	; number of spatial steps
	sigma = n_prandtl  	; Prandtl number
	dz = (1d)/numz   	; spatial step size
	dt = 0.25*dz^2/sigma	; timestep size
	C = c_param 		; planform #
	a = a_param		; wavenumber
	R = n_rayleigh	; Rayleigh number
	t_final = numt	; final timestep
	if ~keyword_set(bc_type) then bc='free' else bc=bc_type


	k = 2.+a^2*dz^2	; a factor that appears in the approximation of w
			; from phi

	alpha = fltarr(nz-2)	; an array of constants that is needed to
				; calculate w's from phi's

	alpha[0] = 1 & alpha[1] = k	; alpha's that are not defined recursively
	FOR i=2,nz-3,1 DO BEGIN	; alpha's that are defined recursively
		alpha[i] = k*alpha[i-1] - 1
	ENDFOR

	;defining variables as 1-d grids
	phi = dblarr(nz,3)
	theta = dblarr(nz,3)
	temp = dblarr(nz,3)
	w = dblarr(nz,3)
	z = dindgen(nz)*dz	; this variable used for plotting purposes

	;setting the initial conditions
	temp[*,0]=1d0-z
	if (bc eq 'rigid') then w[*,0] = (dz/a)*(sin(!pi * z))^2 else w[*,0] = (dz/a)*(sin(!pi * z))
	w[0,0] = 0 & w[nz-1,0] = 0
	theta[*,0] = 0.01*sin(!pi*z)
	theta[0,0] = 0 & theta[nz-1,0] = 0


	; solving for phi from initial w values
	; using phi = d2wdz2 - a^2*w and finite difference approximations
	phi[1:nz-2,0] = (w[0:nz-3,0] - 2*w[1:nz-2,0] + w[2:nz-1,0])/dz^2 - a^2*w[1:nz-2,0]
	if (bc eq 'rigid') then begin
	
		; phi for rigid boundary condition
		phi[0,0] = (2./dz^2)*(w[1,0]) 		; phi simplified by requirement that
		phi[nz-1,0] = (2./dz^2)*(w[nz-2,0]) 	; DW = 0 at boundaries, implying that
							; w[0,*] = w[1,*] and 
							; w=[nz-1,*] = w[nz-2,*]
	endif else begin

		; phi for free boundary condition
		phi[0,0] = 0.		; D^2W = 0 at z=0
		phi[nz-1,0] = 0.	; D^2W = 0 at z=1
	endelse

	; initial conditions will be in 0th row, the most recent timestep in the 2nd
	; row and the timestep before it in the 1st row.  Putting initial conditions
	; in most recent timestep position. This will be copied into the 'old'
	; position at the beginning of time-iterative for loop

	phi[*,2] = phi[*,0]
	theta[*,2] = theta[*,0]
	temp[*,2] = temp[*,0]
	w[*,2] = w[*,0]

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;time evolution

	FOR t=0,t_final DO BEGIN           

		; move last timestep info to 'old' position
		phi[*,1] = phi[*,2]
		theta[*,1] = theta[*,2]
		temp[*,1] = temp[*,2]
		w[*,1] = w[*,2]
	
		; updating phi, theta, temp
		FOR i=1,nz-2 DO BEGIN
	
	
			; the different spatial derivative terms evaluated for each i
			dwdz = (w[i+1,1]-w[i-1,1])/(2d*dz)
			dphidz = (phi[i+1,1]-phi[i-1,1])/(2d*dz)
			d2phidz2 = (phi[i-1,1]-2d*phi[i,1]+phi[i+1,1])/(dz^2)
			dthetadz = (theta[i+1,1]-theta[i-1,1])/(2d*dz)
			d2thetadz2 = (theta[i-1,1]-2d*theta[i,1]+theta[i+1,1])/(dz^2)
			dtempdz = (temp[i+1,1]-temp[i-1,1])/(2d*dz)
			d2tempdz2 = (temp[i-1,1]-2d*temp[i,1]+temp[i+1,1])/(dz^2)
	
			phi[i,2] = phi[i,1] + sigma*dt*( -C/sigma*( 2d*dwdz*phi[i,1] + $
				w[i,1]*dphidz ) + d2phidz2 - a^2*phi[i,1] - R*a^2*theta[i,1] )
	
			theta[i,2] = theta[i,1] + dt*( -C*( 2d*w[i,1]*dthetadz + $
	                   theta[i,1]*dwdz ) + d2thetadz2 - a^2*theta[i,1] - dtempdz*w[i,1] )
	
			temp[i,2] = temp[i,1] + dt*(-theta[i,1]*dwdz - w[i,1]*dthetadz + $
	                  d2tempdz2)
	
		ENDFOR  ;end of time update of phi, theta, temp
	
		if (bc eq 'rigid') then begin
		
			; phi for rigid boundary condition
			phi[0,2] = (2./dz^2)*(w[1,1]) 		; phi simplified by requirement that
			phi[nz-1,2] = (2./dz^2)*(w[nz-2,1]) 	; DW = 0 at boundaries, implying that
								; w[0,*] = w[1,*] and 
								; w=[nz-1,*] = w[nz-2,*]
  		endif else begin

			; phi for free boundary condition
			phi[0,2] = 0.		; D^2W = 0 at z=0
			phi[nz-1,2] = 0.	; D^2W = 0 at z=1

		endelse

		temp[0,2] = 1d & temp[nz-1,2] = 0d
		theta[0, 2] = 0d & theta[nz-1, 2] = 0d

		w[nz-1,2] = 0d 	; boundary condition specifies
		w[0,2] = 0d		; w=0 at z=0,1

		sum_terms = alpha[0:nz-3]*phi[1:nz-2,2]*dz^2

		FOR i=nz-2,1,-1 DO BEGIN
			sum = total(sum_terms[0:i-1])
			w[i,2] = (w[0,2]-sum+w[i+1,2])/alpha[i-2]
		ENDFOR
	

	ENDFOR   ;end of time loop!
  
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   
	consts = create_struct('nz', nz, 'sigma', sigma, 'dz', dz, $
		'dt', dt, 'c', c, 'a', a, 'r', r, 't_final', $
		t_final, 'bc', bc, 'k', k, 'alpha', alpha, name='consts')
	result = create_struct('phi', phi, 'theta', theta, 'temp', $
		temp, 'w', w, 'consts', consts, name='r')

	stop

	return,result

END
