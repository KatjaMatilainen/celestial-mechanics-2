;------------------------------------------------------------------------;
; Use the subroutine PsPlot to save results in a postscript plot 
; (written by Heikki Salo)
;------------------------------------------------------------------------;
pro PsPlot,routine,filename
	thisdir=getenv('PWD')+'/'
	psopen,/color,dir=thisdir,filename
	call_procedure,routine
	psclose		
end

function GetFFTAmplitude,nt,dt,signal

;	table of frequencies
	frequency=lindgen(nt)

;
; this is from the IDL manual, on the ordering of frequencies:
;
	n21=Nt/2L+1L
	frequency[n21]=N21-nt+lindgen(n21-2L)		
	frequency=frequency/(nt*dt)

;
; get the (complex) fourier amplitude
;
	famp=fft(signal,-1)

;
; this is from the IDL manual, so that the most negative frequency comes first
;
	amplitude=shift(abs(famp),-n21)
;
; apply the same shift to frequencies and converting to periods
;
	period=1d0/shift(frequency,-N21)

;
; restrict to positive frequencies (signal is real)
;
	ind=where(period gt 0)
	; factor of 2 in amplitudes: half of amplitude is in
	; negative frequencies for this real signal
	amplitude=2d0*amplitude(ind)
	period=period(ind)

return,{period:period,amplitude:amplitude}
end

;-----------------------------------------
pro FFTExample


	tek_color
	!p.charsize=2

; 	number of time steps
	nt=100000L
;	time increment
	dt=.05d0

	ttab=dt*dindgen(nt)

;
; construct a (real) test signal
;
	period1=10d0
	period2=20d0
	
	omega1=!dpi*2d0/period1
	omega2=!dpi*2d0/period2

	phi1= 1d0*!dpi/3d0
	phi2=0d0*!dpi/17d0

	a1=5d0
	a2=30d0
	signal=a1*cos(omega1*ttab-phi1)+a2*cos(omega2*ttab-phi2)

	res=GetFFTAmplitude(nt,dt,signal)

	window,22
	plot,ttab,signal,xra=[0,150]
	
	window,11
	plot,res.period,res.amplitude,xlog=0,ylog=1,psym=0, $
		xra=[0,50],yra=[1d-20,100],yst=1, $
		xtitle='period',ytitle='FFT amplitude'
;	res= LNP_Test(ttab,signal,wk1=freq,wk2=lomb)
;	nwin 
;	plot,1d0/freq,lomb,xlog=1,ylog=1

end

;--------------------------------------------------------------------;
; Save the results to a PostScript file using PsPlot
;--------------------------------------------------------------------;
pro Plot_everything
PsPlot, 'FFTExample', 'FFTExample.ps'
end
