function gentab,n,min_in,max_in,double=double,openmax=openmax,log=log,power=power
;
; jschmidt 22.05.05
; build an array of given size, minimum, and maximum
;
if n_params() le 0 then begin
print,'function gentab'
print,'tab=gentab(n,min,max,double=double)'
print,' build aan array of given size, minimum, and maximum'
print,''
print,'input:'
print,'n                length of the array'
print,'min              minimum value of the array'
print,'max              maximum value of the array'
print,''
print,'keywords:'
print,'double		return double array'
print,'openmax		return intervall [nmin,nmax) instead of [nmin,nmax]'
print,''
print,'output:          tab'
print,''

return,0
endif

	on_error,2

	min=min_in
	max=max_in

	n=long(n(0))
	if n eq 1 then begin

		tab=min

	endif else begin

		nnorm=n-1
		if keyword_set(openmax) then nnorm=n

		if keyword_set(log) then begin


			  if min le 0 then message, $
				'LOG SCALE ONLY FOR POSITIVE DATA POINTS'
			  if max le 0 then message, $
				'LOG SCALE ONLY FOR POSITIVE DATA POINTS'

			min=alog(min)
			max=alog(max)

		endif

			min=double(min)
			max=double(max)
			tab=dindgen(n)/double(nnorm)*(max-min)+min
	
		if keyword_set(log) then tab=exp(tab)
		if keyword_set(power) then begin 
			apu=(tab-min)/(max-min)
			apu=apu^(1d0/double(power))
			tab=apu*(max-min)+min
		endif

	endelse

return,tab
end
