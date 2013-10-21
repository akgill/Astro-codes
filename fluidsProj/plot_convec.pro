; plotting program
pro plot_convec, ps=ps

    restore,'savfiles/c6a3r1e6npr68v2.sav'

    temp_new = r.temp_new
    w_new = r.w_new
    theta_new = r.theta_new
    nusselt = r.nusselt
    z=findgen(n_elements(temp_new))*1./(n_elements(temp_new)-1)

    ; set plot environment variables
    !P.Multi=[0,1,4,0,1]
    !P.COLOR=255+256*255+(256)^2*255
;    if keyword_set(ps) then !P.COLOR=0

    if keyword_set(ps) then begin
        set_plot,'ps'
        device,filename="thetac6varyamulti2.ps",decomposed=1,/color
    endif

    ymin=fltarr(3) & ymax = [max(theta_new)+0.1*max(theta_new),max(w_new)+0.5*max(w_new),max(temp_new)]

    titlestr = ['theta(z) for Stable Convection with c=6^(-1/2) with a from 0.1 to 9',$
		'w(z) for Stable Convection with c=6^(-1/2) with a from 0.1 to 9',$
		'T_bar(z) for Stable Convection with c=6^(-1/2) with a from 0.1 to 9']

    plot,z,theta_new,xrange=[0.,1.],yrange=[ymin[0],ymax[0]],xtitle="z",$
	ytitle="theta(z)",title=titlestr[0]

    pp = !P & xx = !X & yy = !Y

    plot,z,w_new,xrange=[0.,1.],yrange=[ymin[1],ymax[1]],xtitle="z",$
	ytitle="w(z)",title=titlestr[1]

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    plot,z,temp_new,xrange=[0.,1.],yrange=[ymin[2],ymax[2]],xtitle="z",$
	ytitle="T_bar(z)",title=titlestr[2]

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    plot,z,nusselt,xrange=[0.,1.],/ynozero,xtitle="z",$
	ytitle="Nusselt number",title="Nusselt number"

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    a_params = [0.1, 0.25, 0.5, 0.75, 1, 3, 5, 7, 9]
    inc=fix(256/n_elements(a_params))   
    count = -inc 

    for k=0,n_elements(a_params)-1,1 do begin
        count=count+inc
	print,count

	if k eq 0 then tstep=300000. else tstep=50000.
;	tstep=300000.
	r=convection(100., tstep, 1.e6, 6.8, 6.^(-1/2), double(a_params[k]), t_init=temp_new, w_init=w_new, theta_init=theta_new)
;	r=convection(100., tstep, 8000., 6.8, 6.^(-1/2), double(a_params[k]), /save)

	temp_new = r.temp_new
	w_new = r.w_new
	theta_new = r.theta_new
	nusselt = r.nusselt

	col = color24([(count mod 255),255-(count mod 255),150-(count mod 100)])

	!P = pp[0] & !X = xx[0] & !Y = yy[0]
	oplot,z,theta_new,color=col

	!P = pp[1] & !X = xx[1] & !Y = yy[1]
	oplot,z,w_new,color=col

	!P = pp[2] & !X = xx[2] & !Y = yy[2]
	oplot,z,temp_new,color=col

	!P = pp[3] & !X = xx[3] & !Y = yy[3]
	oplot,z,nusselt,color=col

	print,nusselt
	print,'end of k = ',k
    endfor


    if keyword_set(ps) then device,/close & set_plot,'x'
    !P.Multi=0


    restore,'savfiles/c6a3r1e6npr68v2.sav'

    temp_new = r.temp_new
    w_new = r.w_new
    theta_new = r.theta_new
    nusselt = r.nusselt
    z=findgen(n_elements(temp_new))*1./(n_elements(temp_new)-1)

    ; set plot environment variables
    !P.Multi=[0,1,4,0,1]
    !P.COLOR=255+256*255+(256)^2*255
    if keyword_set(ps) then !P.COLOR=0

    if keyword_set(ps) then begin
        set_plot,'ps'
        device,filename="thetac6varyRmulti2.ps",decomposed=1,/color
    endif

    ymin=fltarr(3) & ymax = [2*max(theta_new),2*max(w_new),max(temp_new)]

    titlestr = ['theta(z) for Stable Convection with c=6^(-1/2) with a from 1e5 to 3e6',$
		'w(z) for Stable Convection with c=6^(-1/2) with R from 1e5 to 3e6',$
		'T_bar(z) for Stable Convection with c=6^(-1/2) with R from 1e5 to 3e6']

    plot,z,theta_new,xrange=[0.,1.],yrange=[ymin[0],ymax[0]],xtitle="z",$
	ytitle="theta(z)",title=titlestr[0]

    pp = !P & xx = !X & yy = !Y

    plot,z,w_new,xrange=[0.,1.],yrange=[ymin[1],ymax[1]],xtitle="z",$
	ytitle="w(z)",title=titlestr[1]

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    plot,z,temp_new,xrange=[0.,1.],yrange=[ymin[2],ymax[2]],xtitle="z",$
	ytitle="T_bar(z)",title=titlestr[2]

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    plot,z,nusselt,xrange=[0.,1.],/ynozero,xtitle="z",$
	ytitle="Nusselt number",title="Nusselt number"

    pp = [pp,!P] & xx = [xx,!X] & yy = [yy,!Y]

    r_params = [1.e5,2.5e5,5.e5,7.5e5,1.e6,2.e6,3.e6]    
    inc=fix(256/n_elements(r_params))   
    count = -inc 

    for k=0,n_elements(r_params)-1,1 do begin
        count=count+inc
	print,count

	if k eq 0 then tstep=300000. else tstep=50000.
;	tstep=300000.
	r=convection(100., tstep, r_params[k], 6.8, 6.^(-1/2), 3., t_init=temp_new, w_init=w_new, theta_init=theta_new)
;	r=convection(100., tstep, 8000., 6.8, 6.^(-1/2), double(a_params[k]), /save)

	temp_new = r.temp_new
	w_new = r.w_new
	theta_new = r.theta_new
	nusselt = r.nusselt

	col = color24([(count mod 255),255-(count mod 255),150-(count mod 100)])

	!P = pp[0] & !X = xx[0] & !Y = yy[0]
	oplot,z,theta_new,color=col

	!P = pp[1] & !X = xx[1] & !Y = yy[1]
	oplot,z,w_new,color=col

	!P = pp[2] & !X = xx[2] & !Y = yy[2]
	oplot,z,temp_new,color=col

	!P = pp[3] & !X = xx[3] & !Y = yy[3]
	oplot,z,nusselt,color=col

	print,nusselt
	print,'end of k = ',k
    endfor


    if keyword_set(ps) then device,/close & set_plot,'x'
    !P.Multi=0

 end


