PRO convection1, numz, numt, n_rayleigh, n_prandtl, c_param, a_param, bc_type=bc_type

  ;assigning constants
  nz = numz+1      	; number of spatial steps
  sigma = n_prandtl  	; Prandtl number
  dz = (1d)/numz   	; spatial step size
  dt = 0.25*dz^2/sigma	; timestep size
  C = c_param 		; planform #
  a = a_param		; wavenumber
  R = n_rayleigh	; Rayleigh number
  t_final = numt	; final timestep



  k = 2.+a^2*dz^2	; a factor that appears in the approximation of w
					; from phi

  alpha = fltarr(nz-2)		; an array of constants that is needed to
				; calculate w's from phi's

  alpha[0] = 1 & alpha[1] = k	; alpha's that are not defined recursively
  FOR i=2,nz-3,1 DO BEGIN	; alpha's that are defined recursively
      alpha[i] = k*alpha[i-1] - 1
  ENDFOR

  ;defining variables as 1-d grids
  phi = dblarr(nz)
  theta = dblarr(nz)
  temp = dblarr(nz)
  W = dblarr(nz)
  z = findgen(nz)*dz	; this variable used for plotting purposes

  ;setting the initial conditions
  temp[0]=1
  W[*] = 0.01*sin(!pi * z)
  W[0] = 0                      
  W[nz-1] = 0
  theta[*] = 0.01*sin(!pi*z)
  theta[0] = 0
  theta[nz-1] = 0



  ;solving for phi from initial W values
  ; using phi = d2Wdz2 - a^2*W and finite difference approximations
  phi[1:nz-2] = (w[0:nz-3] - 2*w[1:nz-2] + w[2:nz-1])/dz^2 - a^2*w[1:nz-2]
  phi[0] = 2d*W[1]/dz^2               ;enforcing rigid boundary conditions
  phi[nz-1] = 2d*W[nz-2]/dz^2

  ; beginning the plotting routine
  ;plot, z, W, xtitle='height (z)', ytitle='W', yrange=[-0.4,1.0]
  ;oplot, z, theta, psym=2
  ;oplot, z, temp, psym=3

  ; defining 'old' and 'new' variable arrays.  Old are from the previous timestep
  ; New are for the next timestep and are calculated below.
  phi_old = phi
  phi_new = phi
  theta_old = theta
  theta_new = theta
  temp_old = temp
  temp_new = temp
  W_old = W
  W_new = W

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FOR LOOP time evolution

FOR t=0,t_final DO BEGIN           

; the 'new' variables from the previous timestep become the 'old' variables
phi_old = phi_new
theta_old = theta_new
temp_old = temp_new
W_old = W_new


  ; updating phi, theta, temp
  FOR i=1,nz-2 DO BEGIN


    ; the different spatial derivative terms evaluated for each i
    dWdz = (W_old[i+1]-W_old[i-1])/(2d*dz)
    dphidz = (phi_old[i+1]-phi_old[i-1])/(2d*dz)
    d2phidz2 = (phi_old[i-1]-2d*phi_old[i]+phi_old[i+1])/(dz^2)
    dthetadz = (theta_old[i+1]-theta_old[i-1])/(2d*dz)
    d2thetadz2 = (theta_old[i-1]-2d*theta_old[i]+theta_old[i+1])/(dz^2)
    dtempdz = (temp_old[i+1]-temp_old[i-1])/(2d*dz)
    d2tempdz2 = (temp_old[i-1]-2d*temp_old[i]+temp_old[i+1])/(dz^2)

    phi_new[i] = phi_old[i] + sigma*dt*( -C/sigma*( 2d*dWdz*phi_old[i] + $
                 W_old[i]*dphidz ) + d2phidz2 - a^2*phi_old[i] - $
                 R*a^2*theta_old[i] )

    theta_new[i] = theta_old[i] + dt*( -C*( 2d*W_old[i]*dthetadz + $
                   theta_old[i]*dWdz ) + d2thetadz2 - a^2*theta_old[i] - $
                   dtempdz*W_old[i] )

    temp_new[i] = temp_old[i] + dt*( -theta_old[i]*dWdz - W_old[i]*dthetadz + $
                  d2tempdz2 )

  ENDFOR  ;end of time update of phi, theta, temp

  ; forcing the rigid boundary conditions
  ;phi_new[0] = 0
  ;phi_new[nz-1] = 0
  temp_new[0] = 1
  temp_new[nz-1] = 0
  theta_new[0] = 0
  theta_new[nz-1] = 0

  W_new[nz-1] = 0. 		; boundary condition specifies
  W_new[0] = 0.		; w=0 at z=0,1

  sum_terms = alpha[0:nz-3]*phi_new[1:nz-2]*dz^2

  FOR i=nz-2,1,-1 DO BEGIN
    sum = total(sum_terms[0:i-1])
    W_new[i] = (W_new[0]-sum+W_new[i+1])/alpha[i-2]
    ;print,w[i,1]
  ENDFOR

  ;enforcing rigid boundaries for phi
  phi_new[0] = 2d*W_new[1]/dz^2
  phi_new[nz-1] = 2d*W_new[nz-2]/dz^2



  ;plot, z, W_new, xtitle='height (z)', ytitle='W', yrange=[0,1]
  ;legend,['W','Theta','Temp'],psym=[0,2,1]
  ;oplot, z, theta_new, psym=2
  ;oplot, z, temp_new, psym=1


ENDFOR   ;end of time loop!
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   

 

stop

END
