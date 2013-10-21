function convection, numz, numt, n_rayleigh, n_prandtl, c_param, a_param, t_init=t_init, w_init=w_init, theta_init=theta_init, plot=plot, save=save

  ;assigning constants
  N = numz+1      ;number of spatial steps
  sigma = n_prandtl  ;Prandtl number
  dz = (1d)/numz   ; spatial step size
  dt = 0.25*dz^2/sigma     ; timestep size
  C = c_param ; planform #
  a = a_param        ; wavenumber
  R = n_rayleigh       ; Rayleigh number
  t_final = numt  ; final timestep

  ;defining variables as 1-d grids
  phi = dblarr(N)
  theta = dblarr(N)
  temp = dblarr(N)
  W = dblarr(N)
  z = findgen(N)*dz          ; this variable used for plotting purposes
  nusselt = dblarr(N)

  ;setting the initial conditions
  if ~keyword_set(t_init) then temp=1.-z else temp = t_init
  if ~keyword_set(w_init) then begin
  	W[*] = (0.01)*(sin(!pi * z))^2 
 	W[0] = 0                      
  	W[N-1] = 0
  endif else W = W_init
  if ~keyword_set(theta_init) then begin
  	theta[*] = 0.01*sin(!pi*z)
  	theta[0] = 0
  	theta[N-1] = 0
  endif else theta=theta_init

  dT_bardz = [0,((temp[2:N-2]-temp[0:N-4])/(2d*dz)),0]

  nusselt = w*theta - dT_bardz

  ;solving for phi from initial W values
  ; using phi = d2Wdz2 - a^2*W and finite difference approximations
  FOR i=1,N-2 DO BEGIN
    phi[i] = (W[i-1] - (2d)*W[i] + W[i+1])/(dz^2) - (a^2)*W[i]
  ENDFOR
  phi[0] = 2d*W[1]/dz^2               ;enforcing rigid boundary conditions
  phi[N-1] = 2d*W[N-2]/dz^2

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
  FOR i=1,N-2 DO BEGIN


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
  ;phi_new[N-1] = 0
  temp_new[0] = 1
  temp_new[N-1] = 0
  theta_new[0] = 0
  theta_new[N-1] = 0
  
  k = 2.+a^2*dz^2	; a factor that appears in the approximation of w
					; from phi

  alpha = fltarr(N-2)		; an array of constants that is needed to
				; calculate w's from phi's

  alpha[0] = 1 & alpha[1] = k	; alpha's that are not defined recursively
  FOR i=2,N-3,1 DO BEGIN	; alpha's that are defined recursively
      alpha[i] = k*alpha[i-1] - 1
  ENDFOR

  W_new[N-1] = 0. 		; boundary condition specifies
  W_new[0] = 0.		; w=0 at z=0,1

  sum_terms = alpha[0:N-3]*phi_new[1:N-2]*dz^2

  FOR i=N-2,1,-1 DO BEGIN
    sum = total(sum_terms[0:i-1])
    W_new[i] = (W_new[0]-sum+W_new[i+1])/alpha[i-2]
    ;print,w[i,1]
  ENDFOR

  ;enforcing rigid boundaries for phi
  phi_new[0] = 2d*W_new[1]/dz^2
  phi_new[N-1] = 2d*W_new[N-2]/dz^2

  dT_bardz = [0,((temp_new[2:N-1]-temp_new[0:N-3])/(2d*dz)),0]

  nusselt = w_new*theta_new - dT_bardz
  nussvar = max(nusselt) - min(nusselt[where(nusselt ne 0)])

  if ((t mod 1000.) eq 0) then print, mean(nusselt), stdev(nusselt), nussvar

;	plot,z,nusselt

	if keyword_set(plot) then begin
  		plot, z, W_new/max(W_new), xtitle='height (z)', ytitle='W', yrange=[0,1]
 		;legend,['W','Theta','Temp'],psym=[0,2,1]
		oplot, z, theta_new/max(theta_new), psym=2
 		oplot, z, temp_new/max(temp_new), psym=1
	endif


ENDFOR   ;end of time loop!
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   


  ;stop
  result = create_struct(name='r','theta_new',theta_new,'w_new',w_new,'temp_new',temp_new, 'nusselt', nusselt, 'nussvar', nussvar)

  if keyword_set(save) then begin
	name = './savfiles/c'+strtrim(string(c))+'a'+strtrim(string(a))+'r'+strtrim(string(R))+npr+strtrim(string(sigma))+'.sav'
	print,name
	save,filename=name,result
  endif

  return, result

END
