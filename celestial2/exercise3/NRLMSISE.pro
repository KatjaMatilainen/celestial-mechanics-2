;;
;; jschmidt Mon Oct  5 15:30:35 EEST 2015 
;;
;; atmoshperic density from
;; the NRLMSIS-00 Model
;; http://ccmc.gsfc.nasa.gov/modelweb/
;;
;;
;---------------------------------------------------
function ReadData,file
;
; read the data file
;

	; open the file
	openr,lun,/get_lun,file

	; skip preliminary lines
	lines_2_skip=16
	line=''
	for ii=0,lines_2_skip-1 do begin
		readf,lun,line
	endfor

	; dummy double numbers to build the arrays
	; of data
	hkm=-1d0
	rhocgs=-1d0

	; a dummy input variable, reading each line
	; that consists of the two doubles
	; (see the data file)
	aa=dblarr(2)
	ff='(f8.1,e9.1)'

	while not eof(lun) do begin
	; loop over entrys of the data file
		readf,lun,aa,f=ff
		; append the arrays successively:
		hkm=[hkm,aa(0)]
		rhocgs=[rhocgs,aa(1)]
	endwhile

	; ignore the first (dummy -1) entry and convert 
	; to meters and kg/m**3
	hh=hkm(1:*)*1d3
	rho=rhocgs(1:*)*1000d0

	; close files
	close,lun
	free_lun,lun

; return the result in a structure
return,{h:hh,rho:rho}
end

;---------------------------------------------------
pro NRLMSISE

	device,decomposed=0,retain=2
	tek_color

	thick=2
	!p.charsize=1.4
	!p.thick=thick
	!x.thick=thick
	!y.thick=thick

;
; read the data
;
	data_file='NRLMSISE-00-Model_Density_vs_Height.txt'
	data=ReadData(data_file)
	title='NRLMSIS-00 Model, from http://ccmc.gsfc.nasa.gov/modelweb/'

;
; plot position
;

	pos=my_tile_to_position(1,1,bottom=3.5,left=1.2,top=1.1)

;
; plot the curve of density vs height
;
	window,1
	plot,data.h/1d3,data.rho,xlog=0,xra=[0,1000],ylog=1, $
		xtitle='ALTITUDE [km]',ytitle='DENSITY [kg/m**3]', $
		title=title,yra=[1d-20,1d0],xst=1,yst=1,pos=pos
;
; overplot in color the ranges of various parts of the atmosphere
;
	xtab=[700,1000,1000,700]
	ytab=[0,0,1,1]
	polyfill,xtab,ytab,/fill,col=2
	xyouts,/norm,ori=90,col=2,0.85,0.05,'EXOSPHERE';,charthick=thick

	xtab=[80,700,700,80]
	ytab=[0,0,1,1]
	polyfill,xtab,ytab,/fill,col=4
	xyouts,/norm,ori=90,col=4,0.4,0.05,'THERMOSPHERE';,charthick=thick

	xtab=[50,80,80,50]
	ytab=[0,0,1,1]
	polyfill,xtab,ytab,/fill,col=3
	xyouts,/norm,ori=90,col=3,0.18,0.05,'MESOSPHERE';,charthick=thick

	xtab=[12,50,50,12]
	ytab=[0,0,1,1]
	polyfill,xtab,ytab,/fill,col=7
	xyouts,/norm,ori=90,col=7,0.15,0.05,'STRATOSPHERE';,charthick=thick

	xtab=[0,12,12,0]
	ytab=[0,0,1,1]
	polyfill,xtab,ytab,/fill,col=6
	xyouts,/norm,ori=90,col=6,0.13,0.05,'TROPOSPHERE';,charthick=thick

;
; also give a typical height for the ISS
;
	hISS=350d0
	oplot,[1,1]*hISS,[1d-50,1d50],lin=2
	xyouts,hISS*1.01,0.96,ori=-90,'ISS',/data

;
; overplot once more, using noerase keyword, to get 
; the tickmaks visible
;
	plot,data.h/1d3,data.rho,xlog=0,xra=[0,1000],ylog=1, $
		yra=[1d-20,1d0],xst=1,yst=1,/noerase,pos=pos

end
