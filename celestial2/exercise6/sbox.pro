;***************************
pro sbox,size,xshift=xshift,yshift=yshift,aspect0=aspect0
;modified 9.11 to take into account the shape of the window

xfac=1.
yfac=1.

if(!d.name ne 'ps') then begin
xs=float(!d.x_size)
ys=float(!d.y_size)

aspect=1.
if(keyword_set(aspect0)) then aspect=aspect0

;ikkuna tavallista leveampi?
if(xs gt 1.25*aspect*ys) then xfac=ys/xs*1.25
;ikkuna tavallista kapeampi?
if(xs le 1.25*aspect*ys) then yfac=xs/ys/1.25
;print,'xs,ys,xfac,yfac'
;print,xs,ys,xfac,yfac
endif

xwidth=0.64*size*1.08939*xfac
ywidth=0.8*size*1.08939*yfac*aspect

x1=.5-xwidth/2.
x2=.5+xwidth/2.
y1=.5-ywidth/2.
y2=.5+ywidth/2.
if(keyword_set(xshift)) then begin
x1=x1+xshift
x2=x2+xshift
endif
if(keyword_set(yshift)) then begin
y1=y1+yshift
y2=y2+yshift
endif

!p.position=[x1,y1,x2,y2]
if(size eq 0) then !p.position=0

return
end
