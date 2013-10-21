FUNCTION helium_saha_func, x

  COMMON params, R_1, R_2

  z_1 = x[0]
  z_2 = x[1]

  z_0 = 1.D0 - z_1 - z_2
  z_e = z_1 + 2.D0*z_2

  f = z_e*[z_1,z_2] - [R_1,R_2]*[z_0,z_1]

  return, f

END

PRO helium_saha, n, T, z_1, z_2, z_0, z_e

  COMMON params, R_1, R_2

  ; Constants

  K_BOLTZ = 1.38d-16
  M_EL = 9.10938D-28
  H_PLANCK = 6.626D-27
  EV_TO_ERG = 1.602D-12

  g_1e0 = 4
  g_2e1 = 1

  chi_1 = 24.6*EV_TO_ERG
  chi_2 = 54.4*EV_TO_ERG

  ; Loop over T

  n_T = N_ELEMENTS(T)

  z_1 = FLTARR(n_T)
  z_2 = FLTARR(n_T)

  x = [0.333d0,0.333d0]

  for i = 0,n_T-1 do begin

     ; Set up the right-hand sides

     R_1 = g_1e0*(2*!PI*M_EL*K_BOLTZ*T[i])^(1.5)/H_PLANCK^3*EXP(-chi_1/(K_BOLTZ*T[i]))/n
     R_2 = g_2e1*(2*!PI*M_EL*K_BOLTZ*T[i])^(1.5)/H_PLANCK^3*EXP(-chi_2/(K_BOLTZ*T[i]))/n

     ; Solve for the ionization state

     x = NEWTON(x, 'helium_saha_func', ITMAX=100, /DOUBLE, CHECK=chk)

     z_1[i] = x[0]
     z_2[i] = x[1]

  endfor

  z_0 = 1 - z_1 - z_2
  z_e = z_1 + 2*z_2

END
