;;
;; jschmidt Tue Nov 10 23:05:25 EET 2015
;;
;; 
;;
;;

;------------------------------------------------
function rhs,tt,vec
;
; equations of motion
;
	common globals,mu,j2,Rp

;
; split up the phase space vector
;
	xx=vec(0)	
	yy=vec(1)	
	zz=vec(2)
	uu=vec(3)	
	vv=vec(4)	
	ww=vec(5)

;
; define aux quantities to save floating point operations
;
	; note: dont do the square as r2=xx^2+yy^2+zz^2
	; (power takes more time than multiplication)
	; (your own implementation of the J2 term
	; will not necessarily use these terms here)
	r2=xx*xx+yy*yy+zz*zz
	rr=sqrt(r2)
	r3=rr*r2
	r4=r2*r2
	r5=r3*r2
	z2=zz*zz
	z2byr2=z2/r2
	
;
; Newtons dynamic law as 1st order ODEs
;
; YOUR TASK 1(a): SPLIT UP THE SECOND DERIVATIVE OF THE RADIUS VECTOR
;		  INTO A SET OF FIRST ORDER DIFFERENTIAL EQUATIONS
;	          FOR THE CARTESIAN COMPONENTS

	; put here the central gravity of a spherical planet

	dxdt = uu
	dydt = vv
	dzdt = ww

	dudt = -1d0*(mu/r3)*xx
	dvdt = -1d0*(mu/r3)*yy
	dwdt = -1d0*(mu/r3)*zz

	if j2 gt 0 then begin
		;
		; include here the effect of non-axisymmetric gravity
		; in terms of the spherical harmonic coefficient j2
		; up to the second order in the Legendre polynomials
		; YOUR TASK 1(b)

		

	endif	

; RETURN HERE A VECTOR OF FIRST DERIVATIVES OF THE CARTESIAN COMPONENTS
; OF THE PHASE SPACE VECTOR. NOTE: IF YOU NAMED THEM DIFFERENTLY IN YOUR
; CODING ABOVE YOU MUST CHANGE THE RETURN STATEMENT ACCORDINGLY. 	
return,[dxdt,dydt,dzdt,dudt,dvdt,dwdt]
end

;------------------------------------------------
function PredictorCorrectorStep,dt,tt,vec,rhs_func

	dvecdt=call_function(rhs_func,tt,vec)
	vec_new_one=vec+dt*dvecdt
	dvecdt2=call_function(rhs_func,tt+dt,vec_new_one)

	vec_new=vec+0.5d0*dt*(dvecdt+dvecdt2)

return,vec_new
end

;------------------------------------------------
;------------------- MAIN -----------------------
;------------------------------------------------
pro PlanetAndSatellites,noj2=noj2

	device,decomposed=0,retain=2
	tek_color
	!p.charsize=2.3
	thick=2
	!x.thick=thick
	!y.thick=thick
	!p.thick=thick

;
; common variables
;
	common globals,mu,j2,Rp

;
; presets of planet and satellite parameters
; Saturn and its moon Mimas
;
	msaturn=5.683d26
	Rsaturn=60268d3
	j2saturn=0.016298d0
	mmimas=3.79d19
	amimas=3.0783d0*Rsaturn

	mjupiter=1.898d27
	Rjupiter=69911d3
	j2jupiter=14736d-6
	mio=893.2d20
	aio=5.91d0*RJupiter

;
; general settings and derived quantities
;	
	
	mode =1 
	if mode eq 1 then begin
	msatellite=mmimas
	asatellite=amimas
	mplanet=msaturn
	Rp=Rsaturn
	j2=j2saturn
	tag='SATURN AND MIMAS'
	endif

	if mode eq 2 then begin
	msatellite=mio
	asatellite=aio
	mplanet=mjupiter
	Rp=Rjupiter
	j2=j2jupiter
	tag='JUPITER AND IO'
	endif


	if keyword_set(noj2) then j2=0d0
	GG=6.673d-11
	mu=gg*(mplanet+msatellite)
	omega_k=sqrt(mu/asatellite^3)
	period=2d0*!dpi/omega_k
	
;
; energy of the orbit
;	
	energy0=-mu/2d0/asatellite

;
; cartesian starting velocities
;
	; this only works if ystart=0 and zstart=0
	delta=5d-3
	ustart=sqrt((j2*(rp/asatellite)^2+delta)*mu/asatellite)
	vstart=omega_k*asatellite*sqrt(1d0-delta)
	wstart=0d0
	v2=ustart*ustart+vstart*vstart+wstart*wstart

;
; cartesiant starting locations
;
	; this only works if ystart=0 and zstart=0
	xstart=asatellite
	ystart=0d0
	zstart=0d0
	rstart2=xstart*xstart+ystart*ystart+zstart*zstart
	rstart=sqrt(rstart2)

;
;  angular momentum and eccentricity
;  (note: this value for angular momentum is the keplerian one, i.e. without j2)
;
	h2=v2*rstart2-(xstart*ustart+ystart*vstart+zstart*wstart)^2
	estart=sqrt((1d0-h2/mu/asatellite)>0d0)
	print,'starting eccentricity:'+string(estart)

;
; time for integration
;
	n_perorbit=100000L
	n_orbits=50L
	tmax=double(n_orbits*period)
	tstart=0d0
	nt=n_perorbit*n_orbits
	dt=tmax/double(nt)

;
; ploting
;
	; define a circle symbol
	phi=!dpi*2d0*dindgen(40)/39d0
	USERSYM,cos(phi),sin(phi),/fill
	symsize=2

	; a window
	window,1,xsize=800,ysize=800
	range=asatellite/rp*(1d0+estart)+1d0
	; plot frame
	plot,[0],[0],xrange=[-1,1]*range,yrange=[-1,1]*range,/nod, $
		xtitle='X [RP]',ytitle='Y [RP]',xst=1,yst=1,title=tag
	; circle denoting the semi-major axis
	oplot,asatellite/rp*cos(phi),asatellite/rp*sin(phi),lin=2
	; circle denoting the apocenter distance
	oplot,asatellite/rp*(1d0+estart)*cos(phi),asatellite/rp*(1d0+estart)*sin(phi),lin=2,col=3
	; circle denoting the pericenter distance
	oplot,asatellite/rp*(1d0-estart)*cos(phi),asatellite/rp*(1d0-estart)*sin(phi),lin=2,col=3
	; circle denoting the planet
	oplot,[0],[0],psym=8,symsize=16,col=7

;
; set starting phase space vector
;
	vec=[xstart,ystart,zstart,ustart,vstart,wstart]	
	rr=sqrt(xstart*xstart+ystart*ystart+zstart*zstart)

	tt=tstart
	xkeep=[xstart]/rp
	ykeep=[ystart]/rp
	
	negative=-1
	positive=1
	rdot=0d0

	print,'PERIOD','|DELTA E /E|','pomega [DEG]' ,f='(a8,a16,a16)'

	for it=1L,nt-1 do begin
;
; main loop
;
		tt=tt+dt
		
		; get updated phase space vector
		vec=PredictorCorrectorStep(dt,tt,vec,'rhs')

		; store radius and sign of drdt from previous step
		rkeep=rr
		rdotkeep=rdot

		update= ((it mod 5) eq 0)
		xx=vec(0)			
		yy=vec(1)			
		zz=vec(2)			
		if update then begin
		; this is for updating the symbol for the satellite
			; first delete the previous symbol (color = 0)
			oplot,xkeep,ykeep,psym=8,col=0,symsize=symsize
			xkeep=[xx/rp]
			ykeep=[yy/rp]
			; plot the symbol at the new location
			oplot,xkeep,ykeep,psym=8,col=2,symsize=symsize
		endif
		r2=xx*xx+yy*yy+zz*zz
		rr=sqrt(r2)
		; check, if we are moving inbound or outbound
		if rr gt rkeep then rdot=positive
		if rr lt rkeep then rdot=negative
		if ((rdot*rdotkeep lt 0) and (rr lt asatellite)) then begin
		; if we are changing the motion from inbound to outbound
		; we are at the pericenter
			; pericenter angle pomega (this is only valid for z=0)
			pomega=atan(yy,xx)
			while pomega lt 0 do pomega=pomega+!dpi*2d0
			; plot a line that denotes the direction to the pericenter
			rtab=dindgen(100)/99d0*range
			xxx=rtab*cos(pomega)
			yyy=rtab*sin(pomega)
			oplot,xxx,yyy,col=9
		endif
		if ((it mod n_perorbit) eq 0) then begin
			v2=vec(3)*vec(3)+vec(4)*vec(4)+vec(5)*vec(5)
			energy=v2/2d0-mu/rr*(1d0-0.5d0*j2*(Rp/rr)^2*(3d0*zz*zz/r2-1d0))
			if it eq n_perorbit then begin
				pomtab=pomega
				ttab=tt
			endif else begin
				pomtab=[pomtab,pomega]
				ttab=[ttab,tt]
			endelse
			; output to screen the number of the orbit, the relative energy change
			; and the pericenter angle in degrees
			print,tt/period,abs((energy-energy0)/energy0),pomega*180d0/!dpi,f='(i8,e16.2,f16.2)'
		endif

	endfor

	pomegadot=1.5d0*omega_k*j2*(rp/asatellite)^2
	print,mean(deriv(ttab,pomtab))/pomegadot

end
