pro co_excitation,ps=ps

;;;;;;;;;;;;;;;;;;;;;;     constants     :::::::::::::::::::::::
	c=299792458*100. ;cm/s
	pi=!dpi
	m_e = 9.11E-28 ;g
	h = 6.6261E-27 ; cgs
	k = 1.381E-16 ;cgs
	sigma_nn = 1.5E-15 ;cm^-2
	T_cmb = 2.73 ;K
	mu = 3.1E-24 ; reduced mass of CO and H2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;; givens ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	A = [7.2, 69., 250., 610., 1200., 2100., 3400.]*10.^(-8) ; s^-1
	nu = [115.3, 230.8, 346.0, 461.5, 576.9, 691.2, 808.5]*10.^9. ;Hz
	E = h*nu ;cgs
	T_gas = [5., 15., 45.] ;K
	J = findgen(8)
	g = 2*J+1
	n_phot = 8*pi*h*nu^3/c^3/(exp(h*nu/k*T_cmb)-1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;; calculated coefficients ;;;;;;;;;;;;;;;;;;;;;;;
	B_ul = c^3/(8*pi*h*nu^3)*A
	B_lu = g[1:n_elements(g)-1]*B_ul/g[0:n_elements(g)-2]
	k_ul = dblarr(n_elements(T_gas),n_elements(nu))
	k_lu = dblarr(n_elements(T_gas),n_elements(nu))

	for i=0,n_elements(nu)-1,1 do begin
		k_ul[*,i] = sqrt(8*k*T_gas/(pi*mu))*sigma_nn
		k_lu[*,i] = (g[i+1]/g[i])*k_ul[*,i]*exp(-E[i]/(k*T_gas))
	endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;; T_exc as a function of n_H2;;;;;;;;;;;;;;;;;
	n_H2 = findgen(1000000)*10.^(0.)
	T_exc = dblarr(n_elements(T_gas),n_elements(nu),n_elements(n_H2))
	for i=0,n_elements(T_gas)-1,1 do begin
		for j=0,n_elements(nu)-1,1 do begin
			ratio = (A[j]+B_ul[j]*n_phot[j]+n_H2*k_ul[i,j])/(B_lu[j]*n_phot[j]+n_H2*k_lu[i,j])
			T_exc[i,j,*]=E[j]/(k*alog(g[j+1]*ratio/g[j]))
		endfor
	endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;; make plots ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	colors=[255.,200.*256.^2+256.*150+200.,200.*256.^2+256.*150.]

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename='ismhwk2.ps',decomposed=1,color=1,bits_per_pixel=8
	endif

	plot,n_H2[1:n_elements(n_H2)-1],t_exc[2,0,1:n_elements(n_H2)-1],/xlog,/ynozero, $
		xstyle=1,xrange=[1.,10.^6.],/ylog,ystyle=1,yrange=[0.4,50.], $
		xtitle="n_H2 [cm^(-3)]",ytitle="T [K]", $
		title="Excitation Temperature as a function of number density"

	for i=0,n_elements(T_gas)-1,1 do begin
		print,"T_gas = ",T_gas[i]," color = ",colors[i]
		for j=0,2,1 do begin
			oplot,n_H2[1:n_elements(n_H2)-1],t_exc[i,j,1:n_elements(n_H2)-1], $
				linestyle=j,color=colors[i]
			print,"/t level = ",j," linestyle = ",j
		endfor
	endfor

	if keyword_set(ps) then device,/close & set_plot,'x'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end