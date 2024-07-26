;;
;; jschmidt Tue Nov 17 18:26:29 EET 2015
;;
;;

;------------------------------------------------
function rhs,tt,vec
;
; equations of motion
;
	common globals,mu1,mu2


;
; split up the phase space vector
;
	xx=vec(0)	
	yy=vec(1)	
	uu=vec(2)	
	vv=vec(3)	

	x1=-mu2
	x2=mu1

        r13=((xx-x1)*(xx-x1)+yy*yy)^(3.d0/2)
        r23=((xx-x2)*(xx-x2)+yy*yy)^(3.d0/2)
        n=1.d0
;
; define aux quantities to save floating point operations
;


; TASK: PUT HERE THE EQUATIONS OF MOTION
;
; RTB as 1st order ODEs
;
        dxdt=uu
        dydt=vv

        dudt=2.d0*n*dydt+n*n*xx+mu1*(x1-xx)/r13+mu2*(x2-xx)/r23
        dvdt=-2.d0*n*dxdt+n*n*yy-mu1*yy/r13-mu2*yy/r23
		
return,[dxdt,dydt,dudt,dvdt]
end

;------------------------------------------------
function PredictorCorrectorStep,dt,tt,vec,rhs_func

	dvecdt=call_function(rhs_func,tt,vec)
	vec_new_one=vec+dt*dvecdt
	dvecdt2=call_function(rhs_func,tt+dt,vec_new_one)

	vec_new=vec+0.5d0*dt*(dvecdt+dvecdt2)

return,vec_new
end

;---------------------------------------------------------------------------
function PseudoPotential,xx,yy,mu1,mu2

	x1=-mu2
	x2=mu1

	y2=yy*yy

	dx1=xx-x1
	r12=dx1*dx1+y2
	r1=sqrt(r12)	

	dx2=xx-x2
	r22=dx2*dx2+y2
	r2=sqrt(r22)	

return,0.5d0*(xx*xx+y2)+mu1/r1+mu2/r2
end

;---------------------------------------------------------------------------
function CJacobi,xx,yy,mu1,mu2
;
; jacobi integral of the circular three body problem,
; evaluated in the whole plane in [xmin,xmax] and [ymin,ymax]
;

	cj=dblarr(n_elements(xx),n_elements(yy))

	for ix=0,n_elements(xx)-1 do begin

		; here for each xx(ix) the whole yy vector
		; is given as input
		cj(ix,*)=2d0*PseudoPotential(xx(ix),yy,mu1,mu2)

	endfor

return,cj	
end

;------------------------------------------------
;------------------- MAIN -----------------------
;------------------------------------------------
pro RTBDynamics,mubar0=mubar0,cj00=cj00,horseshoe0=horseshoe0,tadpole0=tadpole0,inertial=inertial

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
	common globals,mu1,mu2

;
; some parameter presets
;

;
; default starting parameters
;
	mubar=0.02d0
	xstart=-0.6d0
	ystart=0d0
	ustart=0d0
	cj0=3.04d0
	vstart=-sqrt(2d0*PseudoPotential(xstart,ystart,1d0-mubar,mubar)-ustart*ustart-cj0)
	n_perorbit=9800L
	xmin=-4
	xmax=4
	ymin=-4
	ymax=4

;
; overruled by keywords
;
	horseshoe=0
	tadpole=0

	if keyword_set(horseshoe0) and keyword_set(tadpole0) then message,'choose either horseshoe or tadpole'
	if keyword_set(horseshoe0) then begin 
		horseshoe=horseshoe0
		tadpole=0
	endif
	if keyword_set(tadpole0) then begin
		tadpole=tadpole0
		horseshoe=0
	endif

	if horseshoe eq 1 then begin
	mubar=0.000953875d0
	xstart=-0.97668d0
	ystart=0d0
	ustart=0d0
	vstart=-0.06118d0
	cj0=[3.005d0]
	n_perorbit=800L
	if keyword_set(inertial) then n_perorbit=n_perorbit*5L
	xmin=-2
	xmax=2
	ymin=-2
	ymax=2
	endif

	if horseshoe eq 2 then begin
	mubar=0.000953875d0
	xstart=-1.02745d0
	ystart=0d0
	ustart=0d0
	vstart=0.04032d0
	cj0=[3.005d0]
	n_perorbit=800L
	if keyword_set(inertial) then n_perorbit=n_perorbit*5L
	xmin=-2
	xmax=2
	ymin=-2
	ymax=2
	endif

	if tadpole eq 1 then begin
	mubar=1d-3
	x0=0.5d0-mubar
	y0=sqrt(3d0)/2d0
	xstart=x0+0.0065d0
	ystart=y0+0.0065d0
	ustart=0d0
	vstart=0d0
	cj0=[3.0000d0]
	n_perorbit=800L
	if keyword_set(inertial) then n_perorbit=n_perorbit*5L
	xmin=-2
	xmax=2
	ymin=-2
	ymax=2
	endif

	if tadpole eq 2 then begin
	mubar=1d-3
	x0=0.5d0-mubar
	y0=sqrt(3d0)/2d0
	xstart=x0+0.008d0
	ystart=y0+0.008d0
	ustart=0d0
	vstart=0d0
	cj0=[3.0000d0]
	n_perorbit=800L
	if keyword_set(inertial) then n_perorbit=n_perorbit*5L
	xmin=-2
	xmax=2
	ymin=-2
	ymax=2
	endif

;
; assignment of primary masses
;
	mu1=1.d0-mubar
	mu2=mubar

;
; determine for given mu1, mu2 the value of the 
; jacobi integral of the circular three body problem,
; evaluated in the whole plane in [xmin,xmax] and [ymin,ymax]
;
	nx=400
	ny=400
	xtab=gentab(nx,xmin,xmax,/double)
	ytab=gentab(ny,ymin,ymax,/double)
	cjtab=CJacobi(xtab,ytab,mu1,mu2)

;
; the level of the Jacobi integral selected for
; dislpay as a contour
;
	if keyword_set(cj00) then cj0=cj00
	anno=string(cj0,f='(f7.3)')

;
; time for integration
;
	n_orbits=350L
	period=1d0
	tmax=double(n_orbits*period)
	tstart=0d0
	nt=n_perorbit*n_orbits
	dt=tmax/double(nt)

	; define a circle symbol
	phi=!dpi*2d0*dindgen(40)/39d0
	USERSYM,cos(phi),sin(phi),/fill
	symsize=4

;
; first plot
;
	nwin,xs=800,ys=800
	nodata=0
	frame='ROTATING FRAME'
	if keyword_set(inertial) then begin
		 nodata=1
		frame='INERTIAL FRAME'
	endif
	contour,cjtab,xtab,ytab,levels=cj0,c_anno=anno,c_charsize=!p.charsize, $
		c_charthick=thick,xra=[xmin,xmax],yra=[ymin,ymax],xst=1,yst=1,nodata=nodata,$
		xtitle='X',ytitle='Y',title='mu2='+string(mu2,f='(e16.8)')+', '+frame
	plots,-mu2,0,psym=8,col=4,symsize=symsize+4
	plots,mu1,0,psym=8,col=2,symsize=symsize+2

	vec=[xstart,ystart,ustart,vstart]	
	cjstart=2d0*PseudoPotential(xstart,ystart,mu1,mu2)-vec(2)*vec(2)-vec(3)*vec(3)
	xyouts,0.19,0.90,/norm,'CJ='+string(cjstart,f='(f7.3)')

	tt=tstart

;
; keep the starting positions
; to be able to erase them later,
; when the plot is updated
;
	if keyword_set(inertial) then begin
; TASK: PUT HERE THE TRANSFORMATIONS TO THE INERTIAL FRAME 
	; transform back to inertial frame
		xkeep=[xstart]*cos(tt)-[ystart]*sin(tt)
		ykeep=[xstart]*sin(tt)+[ystart]*cos(tt)
	; positions of the primary masses
		x1=-mu2*cos(tt)
		y1=-mu2*sin(tt)
		x2=mu1*cos(tt)
		y2=mu1*sin(tt)
	endif else begin 	
	; keep the plain starting positions
		xkeep=[xstart]
		ykeep=[ystart]
		x1=-mu2
		y1=0d0
		x2=mu1
		y2=0d0

	endelse


	for it=1L,nt-1 do begin
;
; main loop
;
		tt=tt+dt
		
		; get updated phase space vector
		vec=PredictorCorrectorStep(dt,tt,vec,'rhs')

		update= ((it mod 150) eq 0)
		xx=vec(0)			
		yy=vec(1)			
		if update then begin
		; this is for updating the symbol for the test mass
		; first delete the previous symbol (color = 0)
			oplot,xkeep,ykeep,psym=8,col=0,symsize=symsize
		; overplot a new frame with the primary masses
			contour,cjtab,xtab,ytab,levels=cj0,c_anno=anno, $
				c_charsize=!p.charsize,c_charthick=thick,xra=[-1,1]*2, $
				yra=[-1,1]*2,xst=1,yst=1,/noerase,/overplot,nodata=nodata
			plots,x1,y1,psym=8,col=0,symsize=symsize+4
			plots,x2,y2,psym=8,col=0,symsize=symsize+2

			if keyword_set(inertial) then begin
; TASK: PUT HERE THE TRANSFORMATIONS TO THE INERTIAL FRAME 
			; transform back to the inertial frame
				xkeep=[xx]*cos(tt)-[yy]*sin(tt)
				ykeep=[xx]*sin(tt)+[yy]*cos(tt)

				x1=-mu2*cos(tt)
				y1=-mu2*sin(tt)
				x2=mu1*cos(tt)
				y2=mu1*sin(tt)

			endif else begin 	
				xkeep=[xx]
				ykeep=[yy]
				x1=-mu2
				y1=0d0
				x2=mu1
				y2=0d0
			endelse
			; plot the symbol of the test mass at the new location
			oplot,xkeep,ykeep,psym=8,col=7,symsize=symsize
			plots,x1,y1,psym=8,col=4,symsize=symsize+4
			plots,x2,y2,psym=8,col=2,symsize=symsize+2
			empty
		endif
		if ((it mod 5600) eq 0) then begin
		; print progress, accuracy on the screen
		cj=2d0*PseudoPotential(vec(0),vec(1),mu1,mu2)-vec(2)*vec(2)-vec(3)*vec(3)
		print,it,nt,abs((cj-cjstart)/cjstart),f='(2i9,e9.1)'
		endif
	; bail out, if we are far out the plot region
	if abs(xx) gt 5 then break
	if abs(yy) gt 5 then break
	endfor

end
