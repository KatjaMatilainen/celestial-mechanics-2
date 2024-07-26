;***************************************
;  PSCLOSE.PRO
;  closes ps-file
;  directs output to tek
pro psclose,junk=junk,turbops=turbops,lpr=lpr,zille=zille
;***************************************

device,/CLOSE_FILE
;device,/close
;!P.FONT=-1


if(keyword_set(junk)) then begin
file =filepath+'junk'
    if(not keyword_set(zille)) then begin
                                ;if(not keyword_set(lpr)) then 
        spawn,'lpr -h -Pury '+file
                                ;if(    keyword_set(lpr)) then spawn,'lpr -h -Pury junk'
        print,'junk sent to ury'
    endif else begin
        spawn,'lpr -Pzille '+file
        print,'junk sent to zille'
    endelse
endif else begin
    ;print,'ps-plot into file'
endelse

;if(keyword_set(file)) then begin
;      
;   if(not keyword_set(zille)) then begin
;       ;if(not keyword_set(lpr)) then 
;       spawn,'lpr -Pury '+file
;       ;if(    keyword_set(lpr)) then spawn,'lpr -Pury '+file
;       print,file,' sent to ury'
;   endif else begin
;       spawn,'lpr -Pzille '+file
;       print,file+' sent to zille'
;   endelse
;endif
!p.font=-1
set_plot,'X'
device,decomposed=0
end



