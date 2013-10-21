;+
; NAME:
;
; raytrace.pro
;
; PURPOSE: 
;
; To produce spot diagrams from ray tracing through an f/20 cassegrain
; telescope system with a 1m f/3 primary mirror and a confocal
; hyperbolic secondary mirror.
;
; I also just learned the IDL is an Object Oriented Programming language
; and now it seems a lot more cool and useful than it did before.  I've 
; also used this assignment as an exercise for myself in the OO
; related syntax of IDL, making it way more complicated than it actually
; needs to be in order to accomplish it's purpose.  But at least you can
; be sure I didn't copy it from someone else.  The code also should be more 
; adaptable to slightly different optical systems than if I had not written 
; it in an OO way.
;
; CATEGORY:
; Homework #2, Instrumentation, taught by Jason Glenn, Spring 2011
; Department of Astrophysical and Planetary Sciences
; University of Colorado, Boulder
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
; raytrace
;
; CLASSES:
;	LightRay		LightPacket		Mirror
;		define			define			define
;		init			init			init
;		cleanup			cleanup			cleanup
;		set_angle		plot			reflect
;		get_angle		reflectm1		get_diam
;		set_direction		intersectm2		plot
;		get_direction		reflectm2
;		set_position
;		get_position
;		focus
;		plot_ray
;		plot_spot
;	
;	Hyperboloid		Paraboloid
;		define			define
;		init			init
;		cleanup			cleanup
;		get_normal		get_normal
;		get_foc1		get_foc
;		get_foc2		get_z
;		get_a2			plot
;		get_c2
;		get_b2
;		get_dirz
;		plot
;		intersect
;		get_hyp_diam
;
; 
; OTHER METHODS:
;	
; get_mag
; raytrace
;
; MODIFICATION HISTORY:
; v 1.0 written by Amandeep Gill (01/27/2011 - 02/08/2011)
;-


;
; returns the magnitude of a vector, v (used for normalization)
;
FUNCTION get_mag,v
	return,sqrt(v[0]^2+v[1]^2+v[2]^2)
END

; define the LightRay class
;
; dir = a unit vector specifying the direction of the LightRay
; pos = a vector specifying the location of the LightRay with respect
; 	to the origin
; angle = the initial angle of incidence of a LightRay
;
PRO LightRay__DEFINE
	LR = { LightRay, dir:[0d,0d,0d], pos:[0d,0d,0d], angle:0d }
END
;
; LightRay lifecycle functions
;
;
; LightRays must be initialized with a vector specifying position and can optionally
; be initialized with an angle (radians) and with a vector specifying direction.  
; Default angle is 0 and default direction is [0, sin(self.angle), cos(self.angle)] 
; which equals [0, 0, 1] when angle equals 0. 
;
FUNCTION LightRay::init, position_vector, angle=angle, direction_vector=direction_vector
	self.angle = 0d
	if keyword_set(angle) then self.angle=angle
	
	self.dir = [0,sin(self.angle),(cos(self.angle))]
	if keyword_set(direction_vector) then self.dir = direction_vector/get_mag(direction_vector)

	self.pos = position_vector
	return,1
END
;
PRO LightRay::cleanup
END
;
; LightRay accessor methods
;
PRO LightRay::set_angle, t
	self.angle = double(t)
END
;
FUNCTION LightRay::get_angle
	return,self.angle
END
;
PRO LightRay::set_direction, dir_vect
	if (get_mag(dir_vect) eq 0) then self.dir = dir_vect else begin
		self.dir = dir_vect/get_mag(dir_vect)
	endelse
END
;
FUNCTION LightRay::get_direction
	return,self.dir
END
;
PRO LightRay::set_position, pos_vect
	self.pos = pos_vect
END
;
FUNCTION LightRay::get_position
	return,self.pos
END
;
; the focus function propagates a LightRay along its direction from its
; position until it reaches the "focal_z"
;
FUNCTION LightRay::focus, focal_z

; defining components of the direction vector, v, the position vector, and defining
; the propagation parameter, t, such that LightRay propagates to focal_z
	vx = self.dir[0] & vy = self.dir[1] & vz = self.dir[2]
	x0 = self.pos[0] & y0 = self.pos[1] & z0 = self.pos[2]
	t = (focal_z - z0)/vz

; defining the final position coordinate of the propagated LightRay
	x = x0 + vx*t
	y = y0 + vy*t
	z = focal_z

	return,[x, y, z]
END

; plot_ray plots a LightRay from it's current location to zf, optional parameter
; oplot allows over plotting an existing plot
;
PRO LightRay::plot_ray, zf, oplot=oplot
	xi = self.pos[0] & yi = self.pos[1] & zi = self.pos[2]
	vx = self.dir[0] & vy = self.dir[1] & vz = self.dir[2]
	t = (zf - zi)/vz
	yf = yi + vy*t & xf = xi + vx*t

	if ~keyword_set(oplot) then plot,[zi,zf],[yi,yf],linestyle=0 else begin
		oplot,[zi,zf],[yi,yf],linestyle=0
	endelse
END

; plot_spot plots the [x, y] coordinate of a LightRay at its current z position
;
PRO LightRay::plot_spot, oplot=oplot
	if keyword_set(oplot) then oplot,[self.pos[0]],[self.pos[1]],psym=2 $
	else plot,[self.pos[0]],[self.pos[1]],psym=2
END


; Defining the LightPacket class
;
;
; angle = angle of initialization of LightRays in LightPacket (in arcseconds)
; packet = an 11 x 11 array of LightRays
; mirror = points to the primary mirror object LightPacket is initialized on
;
;
PRO LightPacket__DEFINE
	LP = { LightPacket, angle:0d, packet:objarr(11,11), mirror1:obj_new(), mirror2:obj_new()}
END
;
; LightPacket lifecycle functions
;
; LightPacket must be initialized with an array of pointers to the Mirrors in the
; system
;
; LightPackets can optionally be initialized with an angle which must be
; specified in arcseconds.  The constructor converts it to radians.
;
FUNCTION LightPacket::init, mirror_pointer, angle=angle
	self.mirror1=mirror_pointer[0]
	self.mirror2=mirror_pointer[1]
	self.angle = 0d
	if keyword_set(angle) then self.angle=double(angle*!dpi/(3600d*180d))
	
	dx = self.mirror1->get_diam()/(sqrt(2d)*10d)
	for i=0d,10d,1d do for j=0d,10d,1d do begin
		x = -5d*dx + i*dx & y = -5d*dx + j*dx
		z = self.mirror1->get_z(x, y)
		pos = [x, y, z]
		self.packet[i,j] = obj_new("LightRay", pos, angle=self.angle)
	endfor
	return,1
END
;
PRO LightPacket::cleanup
END

FUNCTION LightPacket::get_angle
	return,(3600d*180d*self.angle/!dpi)
END

PRO LightPacket::plot, initial=initial
	if keyword_set(initial) then begin
		self.mirror1->plot,/faceon
		for i=0,10,1 do for j=0,10,1 do self.packet[i,j]->plot_spot, /oplot
	endif
END

PRO LightPacket::reflectm1
	for i=0,10,1 do for j=0,10,1 do begin
		oldpos=self.packet[i,j]->get_position()
		olddir=self.packet[i,j]->get_direction()
		self.packet[i,j]->set_direction,(self.mirror1->reflect(oldpos, olddir))
	endfor
END

FUNCTION LightPacket::intersectm2
	intersections=dblarr(11, 11, 3)
	for i=0,10,1 do for j=0,10,1 do begin
		oldpos=self.packet[i,j]->get_position()
		olddir=self.packet[i,j]->get_direction()
		intersections[i, j, *] = self.mirror2->intersect(oldpos, olddir)
	endfor
	return,intersections
END

PRO LightPacket::set_positions,pos_arr
	for i=0,10,1 do for j=0,10,1 do begin
		self.packet[i, j]->set_position,pos_arr[i,j,*]
	endfor
END

PRO LightPacket::reflectm2
	for i=0,10,1 do for j=0,10,1 do begin
		oldpos=self.packet[i,j]->get_position()
		olddir=self.packet[i,j]->get_direction()
		self.packet[i,j]->set_direction,(self.mirror2->reflect(oldpos, olddir))
	endfor
END

PRO LightPacket::plot_packet_rays, zarr
	for i=0,10,1 do for j=0,10,1 do begin
		pos = self.packet[i,j]->get_position()
		dir = self.packet[i,j]->get_direction()
		xi = pos[0] & yi = pos[1] & zi = pos[2]
		vx = dir[0] & vy = dir[1] & vz = dir[2]
		xf = zarr[i,j,0] & yf = zarr[i,j,1] & zf = zarr[i,j,2] 
		oplot,[zi,zf],[yi,yf],linestyle=2
;		self.packet[i,j]->plot_ray,zarr[i,j],/oplot
	endfor
END

PRO LightPacket::plot_packet_spots, spots_pos
	xmin = min(spots_pos[*,*,0]) & xmax = max(spots_pos[*,*,0]) & xoff = (xmax-xmin)/10d
	ymin = min(spots_pos[*,*,1]) & ymax = max(spots_pos[*,*,1]) & yoff = (ymax-ymin)/10d

	plot,[spots_pos[*,*,0]],[spots_pos[*,*,1]],psym=2, $
	xrange=[xmin-xoff,xmax+xoff], yrange=[ymin-yoff,ymax+yoff], xstyle=1, ystyle=1
END

FUNCTION LightPacket::focus, focal_z
	packet_pos = dblarr(11, 11, 3)
	for i=0,10,1 do for j=0,10,1 do begin 
		packet_pos[i,j,*] = self.packet[i,j]->focus(focal_z)
	endfor
	return,packet_pos
END

FUNCTION LightPacket::rms, z
	center = self.find_center(z)
	ms=0d
	for i=0,10,1 do for j=0,10,1 do begin 
		point = self.packet[i,j].focus(z)
		dist = sqrt((point[0]-center[0])^2 + (point[1]-center[1])^2)
		ms = ms + dist^2
;		print,"dist: ",dist
	endfor
	rms = sqrt(ms/121d)

	return,rms
END

FUNCTION LightPacket::find_center,z
	xsum = 0d & ysum = 0d

	for i=0,10,1 do for j=0,10,1 do begin 
		point = self.packet[i,j]->focus(z)
		xsum = xsum + point[0]
		ysum = ysum + point[1]
	endfor
	
	center = [xsum/121d, ysum/121d]
	return,center
END

FUNCTION LightPacket::optimize_focal_z, z_init
	threshhold = (10d)^(-15d)
	zlo = 0d
	zhi = double(z_init)
	repeat begin
		rmshi = self.rms(zhi) & rmslo = self.rms(zlo)
		if (rmshi lt rmslo) then zlo = (zhi+zlo)/2d else zhi = (zhi+zlo)/2d
	endrep until (rmslo le threshhold) or ((zhi - zlo) le threshhold)
	return,zlo
END

; Defining the Mirror class - Mirror is a parent class of Paraboloid and Hyperboloid
;
;
; diam = Mirror diameter
;
PRO Mirror__DEFINE
	M = {Mirror, diam:0d}
END
;
; Mirror lifecycle functions
;
FUNCTION Mirror::init, diameter
       self.diam=double(diameter)
       return, 1
END

PRO Mirror::cleanup
END

; reflect calculates the trajectory of an outgoing light ray given the
; direction vector of the incoming light ray and it's position of incidence
; 
; Even though Mirror classes do not have get_normal methods, the daughter
; classes do, and their methods are used to properly calculate reflection.
;
FUNCTION Mirror::reflect, light_pos, light_dir
	pos = light_pos
	dir = light_dir
	norm = self->get_normal(pos[0], pos[1], pos[2])
	
	dotprod = (norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2])
	new_dir = dir - 2d*dotprod*norm
	new_dir_unit = new_dir/get_mag(new_dir)
	return, new_dir_unit
END

; Mirror accessor function
;
FUNCTION Mirror::get_diam
	return,self.diam
END

PRO Mirror::plot, oplot=oplot
END


; Defining the Paraboloid class - daughter of the Mirror class
;
; F = focal length
; diam = diameter
; 
;
PRO Paraboloid__DEFINE
	P = { Paraboloid, F:0d, inherits Mirror}
END
;
; Paraboloid lifecycle functions
;
;
; Paraboloids must be initialized with specified focal length, and diameter
;
FUNCTION Paraboloid::init, focal_length, diameter
	self.F = double(focal_length)
	self.diam = double(diameter)
	return, 1
END
;
PRO Paraboloid::cleanup
END
;
; calculates the normal to a paraboloid oriented with rotational symmetry
; about the z axis, opening up in the -z direction, given a position for
; which to find the normal
;
; returns the gradient vector in the cartesian basis
;
FUNCTION Paraboloid::get_normal, x, y, z
	grad=[(-x/(2d*self.F)), (-y/(2d*self.F)), (-1d)]
	return, (grad/get_mag(grad))
END
;
; Paraboloid accessor methods
;
FUNCTION Paraboloid::get_foc
	return,self.F
END

FUNCTION Paraboloid::get_z, x, y
	return, (-(x^2+y^2)/(4d*self.F))
END

PRO Paraboloid::plot, oplot=oplot, faceon=faceon
	if (keyword_set(faceon)) then begin
		x = dindgen(1000)/100d - 1d
		y = sqrt(0.25d - (x^2))
		plot,x,y,xrange=[-1,1],yrange=[-1,1]
		y = -sqrt(0.25d - (x^2))
		oplot,x,y
	endif else begin
		;z = (x^2 + y^2)/(4d*self.F)
		z = -dindgen(26)/(25d*48d)
		y = -sqrt(-z*4d*self.F)
;		print,self.F
;		print,z,y
		plot,z,y,xrange=[-3.5,0.5],yrange=[-0.7,0.7],xstyle=1,ystyle=1
		y = sqrt(-z*4d*self.F)
		oplot,z,y
	endelse
END
;
; Defining the Hyperboloid class - daughter of the Mirror class
;
; for a Hyperboloid defined by the equation (z-z0)^2/a^2 - (x^2+y^2)/b^2 = 1
; where z0 is the location of the center of the hyperbola, a is the distance
; from the center to the hyperbola vertex, c is the distance of the center to a
; focus, and a^2 + b^2 = c^2
;
; F1 = absolute value of the z position of focus before the hyperboloid
; F2 = absolute value of the z position of the focus after the hyperboloid
; a2 = a^2 in the hyperbola definition
; c2 = c^2 in the hyperbola definition
; Z = z0 in the hyperbola definition
; diam = diameter
;
;
PRO Hyperboloid__DEFINE
	H = { Hyperboloid, F1:0d, F2:0d, a2:0d, c2:0d, Z:0d, inherits Mirror}
END
;
; Hyperboloid lifecycle functions
;
;
FUNCTION Hyperboloid::init, focal_length1, focal_length2, asq, fnumber
	self.F1 = abs(double(focal_length1))
	self.F2 = abs(double(focal_length2))
        self.a2 = double(asq)
        self.c2 = ((abs(self.F1)+abs(self.F2))/(2d))^2
	self.Z = ((abs(self.F1)+abs(self.F2))/(2d))-self.F1
	fnum=double(fnumber)
	self.diam = self.get_hyp_diam(self.F1, self.F2, fnum)
	return, 1
END
;
PRO Hyperboloid::cleanup
END
;
; Hyperboloid accessor methods
;
FUNCTION Hyperboloid::get_foc1
	return,self.F1
END

FUNCTION Hyperboloid::get_foc2
	return,self.F2
END

FUNCTION Hyperboloid::get_a2
	return,self.a2
END

FUNCTION Hyperboloid::get_c2
	return,self.c2
END

FUNCTION Hyperboloid::get_b2
	return,(self.c2-self.a2)
END

FUNCTION Hyperboloid::get_dirz
	return,self.Z
END

PRO Hyperboloid::plot, oplot=oplot
END
;
; calculates the normal to a paraboloid oriented with rotational symmetry
; about the z axis, opening up in the -z direction
;
; returns the gradient vector in the cartesian basis
;
FUNCTION Hyperboloid::get_normal, x, y, z
	a2 = self.a2 
	b2 = self.c2-self.a2
;	z = - sqrt((a2)*(1d + (x^2 + y^2)/(b2))) + self.Z
	grad=[(2d*x/(b2)), (2d*y/(b2)), (-2d*(z-self.Z)/(a2))]
	return, (grad/get_mag(grad))
END
;
; calculates the position at which a LightRay with given initial position and
; direction of travel will intersect the Hyperboloid
;
; for the relevant mathematics, see attached work
;
FUNCTION Hyperboloid::intersect, light_pos, light_dir
	a2 = self.a2 & b2 = self.c2-self.a2
	vx = light_dir[0] & vy = light_dir[1] & vz = light_dir[2]
	x0 = light_pos[0] & y0 = light_pos[1] & z0 = light_pos[2]
	z1 = self.Z

	alpha = vz^2/a2 - vx^2/b2 - vy^2/b2
	beta = 2d*(vz*z0/a2 - vz*z1/a2 - vx*x0/b2 - vy*y0/b2)
	gamma = (z0^2 - 2d*z0*z1 + z1^2)/a2 - (x0^2 + y0^2)/b2 - 1d

	tsol1 = (-beta + sqrt(beta^2 - 4d*alpha*gamma))/(2d*alpha)
	tsol2 = (-beta - sqrt(beta^2 - 4d*alpha*gamma))/(2d*alpha)

	if (abs(tsol1) gt abs(tsol2)) then t = tsol1 else t = tsol2
;	print,t

	intersection = [(x0+vx*t), (y0+vy*t), (z0+vz*t)]
	return, intersection
END
;
;
; calculates the diameter of the Hyperboloid based on the location of
; the foci and on the system f#
;
FUNCTION Hyperboloid::get_hyp_diam, f1, f2, fnum
	f1=f1 & f2=f2 & fnum=fnum
	c=(f1+f2)/(2d) & z0 = c-f1
	x1=0d & y1=0.5d & z1=-(x1^2d + y1^2d)/12d
	xf1 = 0d & yf1 = 0d & zf1 = -f1
	vx = x1 - xf1 & vy = y1 - yf1 & vz = z1 - zf1
; y = y1 + vy*t
; t = (y - y1)/vy
; z = z1 + vz*t
; z = z1 + (vz/vy)(y - y1)

	l1m = (vz/vy)
	l1b = 2d*(-1d/48d + 3d)*(-0.5d) - 1d/48d ;z1 - l1m*y1

; D = 2y
; F = z0 + c - z
; fnum = F/D
; 2*fnum*y = z0 + c - z
; z = -2*fnum*y +z0 + c

	l2m = -2d*fnum
	l2b = z0+c

	yint = (l2b-l1b)/(l1m-l2m)
	zint = l2m*yint + l2b

        return, 2d*yint
	
END






; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; 
; ; ;     raytrace is the driver program    ; ;
; ; ;    THIS IS WHERE THE MAGIC HAPPENS!   ; ;
; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ; ;
;
PRO raytrace

; here I create the optical system

; create a Paraboloid object with focal length = 3 and diameter = 1
	P = obj_new("Paraboloid",3d,1d)

; create a Hyperboloid object with first focal point at -3 (my methods 
; use the absolute value of focal positions), the second focal point at 
; 0.1, an a_squared value I calculated outside of IDL and the system 
; fnumber = 20
        a_squared_hyp = ((961d*(943201d - 6960d*sqrt(1601d)))/(486643600d))
	fnum=20d
	H = obj_new("Hyperboloid", 3d, 0.1d, a_squared_hyp, fnum)
	
; plotting stuff
	!P.Multi=0
	!P.Multi=[0,1,2,0,2]

	set_plot,'ps'
;	set_plot,'x'

;
; create an array of LightPackets, each with different incidence angles	
;
	angs = [0, 10, 30, 120, 600, 900, 1800]
;	angs = dindgen(180)*10d
	numpacks = n_elements(angs)
	packets = objarr(numpacks)
	for i=0, numpacks-1, 1 do packets[i] = obj_new("LightPacket", [P,H], angle=angs[i])
;	packets[0]->plot,/initial

	avg_plate_scale=0d

	foc_plane = dblarr(numpacks,2)

	for i=0, numpacks-1, 1 do begin

; reflect packets off of the primary mirror
		packets[i]->reflectm1

; find the intersections of the reflected packet with the secondary mirror
		intersections = packets[i]->intersectm2()
;		print,intersections

; file to plot to
	num=strtrim(angs[i],2)
	filename="spotdiag_z_opt"+num+".ps"
	device,file=filename

; plot the primary and the path of the primary-reflected packets for i=0
;		if (i eq 0) then begin
			P.plot
			packets[i]->plot_packet_rays,intersections
;		endif

; update packet positions to hyperboloid intersections
		packets[i]->set_positions,intersections

; reflect packet off of the secondary mirror
		packets[i]->reflectm2

; get packet positions at specified z and plot rays
		z_slice_pos = packets[i]->focus(H.get_foc2())
		z_opt = packets[i]->optimize_focal_z(z_slice_pos[0,0,2])
		zoptrms = packets[i]->rms(z_opt)
		rmsstring = "rms= "+strtrim(zoptrms)
		xyouts,-3,-1.3,rmsstring
		packets[i]->plot_packet_rays,packets[i]->focus(z_opt)
;		wait,1

; plot spot diagrams, print out rms, calculate plate scale, and optimal focal_z
; calculated plate scale by measuring the offset of the lightpacket centroid
; from 0 in mm and then dividing that number into the lightpacket angle (in arcseconds)

		packets[i]->plot_packet_spots,packets[i]->focus(z_opt)
		center = packets[i]->find_center(z_slice_pos[0,0,2])
		angle = packets[i]->get_angle()
		ps = angle/(sqrt(center[0]^2+center[1]^2)*1000d)
		avg_plate_scale = avg_plate_scale+ps
		print,"angle:",angle," rms: ",packets[i]->rms(z_slice_pos[0,0,2]),$
		" plate scale: ",ps," mm/arcsec  z_opt: ",z_opt," center: ",center
;		wait,3

		foc_plane[i,*] = [angle,z_opt]

; close plotted file
	device,/close	

	endfor
	avg_plate_scale = avg_plate_scale/(n_elements(angs)-1)
	print,"avg. plate scale: ",avg_plate_scale

	fit = poly_fit(foc_plane[*,0],foc_plane[*,1],2,/double)
	print,""
	print,"Parabolic Fit to Focal Plane:"
	print,"optimal z = ",fit[2]," * angle^2 ",fit[1]," * angle + 0.1"

END
