;------------------------------------------------------------------;
; Celestial mechanics II - Exercise 3
;------------------------------------------------------------------;
; 3 - 2a) Harmonic oscillator
; What to do: Make a copy of HarmOscillator.pro and modify it to
; produce a numerical solution of the simplified equation for the
; orbital decay of the (normalized) semimajor axis

; dalpha/dt=-sqrt(alpha)

; Also overplot the analytical solution we had obtained
; alpha=(1-tau/2)^2

; Hint: define scaled time tau in the range [0,2]

;------------------------------------------------------------------;
;;
;; jschmidt Wed Oct  7 09:40:06 EEST 2015
;;
;; test routine for simple integration of 
;; ordinary differentialequations (ODE)
;;
;; solve the equation for an undamped harmonic 
;; oscillator numerically and compare to 
;; the analytic solution
;;
;;
;; first use an Euler Step for the time advancing
;; then try a precictor-corrector step 
;;

;------------------------------------------------
function rhs,tt,vec
;
; This is the function that povides the Right Hand Sides 
; of the differential equation at hand, i.e.:
; df/dt=RHS
; note that f may be a vector, which is the case here
;
; in fact for the harmonic oscillator we have 
; d(dx/dt)/dt =- x (assuming unit frequency)
; this is a SECOND ORDER ODE, i.e. their occurs
; a second derivative with respect to time
;

;-------------------------------------------------------;
; MODIFY: dalpha/dt=-sqrt(alpha)
; Left-hand side is ok.
; RHS: -sqrt(alpha) in parametric form, how??
;-------------------------------------------------------;
  
; this can be written as two FIRST ORDER ODE by
; denoting with v=dx/dt the velocity 
; Thus we can write the harmonic oscillator as
; dx/dt=v
; dv/dt=-x
; 
; thus we can write f=[x,v] and RHS=[v,-x]
;
; because our inpout vector vec is representing f
; we have
;
	xx=vec(0)
	vv=vec(1)
;        vv=-sqrt(vec(1))
        
;        dxdt=vv
        dxdt=-sqrt(vec(0))
;	dvdt=-xx
        dvdt=-sqrt(vec(1))
;        dvdt=0.d0*vec(1)
        
; this routine should return the RHS so
return,[dxdt,dvdt]
end

;------------------------------------------------
function EulerStep,deltat,tt,vec,rhs_func
;
; this is the simplest possible numerical time
; advancing method for ODE
;
; unfortunately it is unstable.  but see yourself.
;
; It uses the first order taylor expansion of a 
; function to approximate the time derivative
;
; with Taylors theorem we write
; f(t+delta t)=f(t)+df/dt delta t
; the df/dt we get from the function rhs
; that provides the right hand sides
	dvecdt=call_function(rhs_func,tt,vec)
; so that the function at the later time
; is approximated by:
	vec_new=vec+deltat*dvecdt

return,vec_new
end

;------------------------------------------------
function PredictorCorrectorStep,deltat,tt,vec,rhs_func
;
; this is a stable method
;
; it first computes, from the rhs, the derivative
; at the given time tt (just as for the Euler step):
	dvecdt=call_function(rhs_func,tt,vec)
	
; then it computes an updated solution for the function
; from an Euler step
	vec_new_one=vec+deltat*dvecdt
; but this is not stable, as we have heard
; to do better, we now compute with rhs
; the derivative at tt+deltat
	dvecdt2=call_function(rhs_func,tt+deltat,vec_new_one)

; now we go back, and advance the function at time tt
; by using the AVERAGE of the two derivatives:
	vec_new=vec+0.5d0*deltat*(dvecdt+dvecdt2)

; the drawback is that this needs the double amount of computations.
; the advantage is that we have a stable, thus useful, scheme for
; time advancing instead of a faster but useless one.
return,vec_new
end


;------------------------------------------------
;pro HarmOscillator
pro exercise3
;
; some preliminary settings. don't worry about this.
;
	device,decomposed=0,retain=2
	tek_color
	thick=2
	!x.thick=thick
	!y.thick=thick
	!p.thick=thick

;
; this sets the period of the oscillation to 2 Pi. 
; thus the frequency is unity
;
	period=!dpi*2d0
;
; choose a time step. we like to resolve the oscillations
; over one period, so having a time step of one percent
; of the period should be enough. but change this number
; and look what happens.
;
	dt=period/100d0

;
; initial conditions, starting velocity is zero, starting at x=1 at zero time
;
	vv=0d0
	xx=1d0
	tt=0d0

;
; this is the analytic solution
;
	; a time array of 10 periods length
	; having 1000 entrys altogether
;	ttt=dindgen(1000)/999d0*10d0

;MODIFIED: scaled time [0,2]
        ttt=dindgen(2001)/1000.d0

; the known analytic solution for these
	; starting conditions
;	xxx=cos(ttt*2d0*!dpi)
;MODIFIED: for alpha
        xxx=(1.d0-ttt/2.d0)^2.d0

;
; plot the analytic solution
;
	window,1
	plot,ttt,xxx,/nod,xtitle='TIME',ytitle='XX';,xrange=[0,2]

;
; overplot, with symbols, the instantanous state of our numerical
; solution. this is for the moment tt=0 and xx=1:
	plots,tt/period,xx,psym=4

	while tt lt period*10d0 do begin
;
; this is the main loop. now we combine xx and vv into one vector:
;		
		vec=[xx,vv]
;
; and pass it to the routine that does the time advancing:
;
; TRY HERE ALTERNATIVELY THE TWO METHODS FOR TIME ADVANCING
; COMMENTING OUT ONE OF THEM AND RE-RUNNING THE ROUTINE
		res=PredictorCorrectorStep(dt,tt,vec,'rhs')
		;res=EulerStep(dt,tt,vec,'rhs')
;
; splitting up the result again in xx and vv, we have the 
; updated, time advanced, solutions:
;
		xx=res(0)
		vv=res(1)

;
; incrementing the time with the timestep
;
		tt=tt+dt	
;
; we can now plot our next symbol for the numerical sulution:
;
		plots,tt/period,xx,psym=4
		
	endwhile
;
; we overplot the analytic solution once more
; (because this looks better if the curve is on top
; of the symbols, not vice versa
;
	oplot,ttt,xxx,col=2

stop
end
