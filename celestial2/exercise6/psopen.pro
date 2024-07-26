;***************************************
;  PSOPEN.PRO
;  directs output to ps-file
;  no color
   pro psopen,file,eps=eps,ctek=ctek,$
   color=color,vfont=vfont,sbox=sbox,copy=copy,up=up,right=right,$
   help=help,info=info,portrait=portrait,dir=dir

;**************************************************************

if(n_params() le 0 or keyword_set(info) or keyword_set(help)) then begin
print,' pro psopen,file,eps=eps,ctek=ctek,color=color,vfont=vfont'
print,' sbox=sbox,copy=copy,up=up,right=right,help=help,info=info'
print,' '
print,' direct output to file in landscape mode'
print,' /eps -> encapsulated ps :asks for more info'
print,' /color -> color ps'
print,' /ctek -> good at ctek-printer ?'
print,' /vfont -> use vector-fonts'
print,' /sbox -> force square shaped plot'
print,' /copy -> copy current colortable stretch (def=no)'
print,' up-> move upward by UP cm '
print,' right-> move right by RIGHT cm'
print,' /info or /help -> this message'
print,' HS 7.8.96  /appl1/heikki/UNAMIDL/psopen.pro'
print,' '
return
endif 
;**************************************************************

icopy=0
device,decomposed=1
if(keyword_set(copy)) then icopy=1
set_plot,'ps',copy=icopy


filepath='~/tmp/'
if keyword_set(dir) then filepath = dir


if not keyword_set(portrait) then begin
    device,FILENAME=filepath+file,/LANDSCAPE,encapsulated=0,bits_per_pixel=8,color=0,scale=1.
endif
if keyword_set(portrait) then begin
    device,FILENAME=filepath+file,/portrait,bits_per_pixel=8,encapsulated=0,color=0,scale=1.
endif

; good values for sparc
xoff=1
yoff=24
;good values for ury
yoff=28
if(keyword_set(up)) then xoff=xoff+up
if(keyword_set(right)) then yoff=yoff-right


if not keyword_set(portrait) then device,xsize=27.,ysize=18.,yoff=yoff,xoff=xoff
if keyword_set(portrait) then device,ysize=25.5,xsize=17.,yoff=2.,xoff=2.

if(keyword_set(ctek)) then device,yoff=27.
if(keyword_set(color)) then device,/color
sbox,0
if(keyword_set(sbox)) then sbox,1


if(keyword_set(eps)) then begin
device,/ENCAPSULATED
print,'ENCAPSULATED POSTRIPT CHOSEN'
if n_elements(eps) eq 1 then begin
print,'PORTRAIT-MODE (1-yes)'& read,vast & if(vast eq 1) then device,/PORTRAIT
print,'XSIZE (cm) ' & read,xsize & device,XSIZE=xsize
print,'YSIZE (cm) ' & read,ysize & device,YSIZE=ysize
print,'XOFF (cm) ' & read,xoff & device,XOFFSET=xoff
print,'YOFF (cm) ' & read,yoff & device,YOFFSET=yoff
endif else begin
if eps(0) eq 1 then device,/portrait
XSIZE=eps(1) & device,XSIZE=xsize
ysize=eps(2) & device,YSIZE=ysize
xoffset=eps(3) & device,XOFFSET=xoff
yoffset=eps(4) & device,YOFFSET=yoff
endelse
endif

!P.FONT=0
;device,/TIMES, /BOLD
;device,/COURIER, /BOLD
;device,/HELVETICA, /BOLD
device,/HELVETICA	
;device,/SCHOOLBOOK, /BOLD	

if(keyword_set(vfont)) then !p.font=-1
end







