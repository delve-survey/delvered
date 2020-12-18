pro delvered_refcat_prep,exposure,refcatfile,offset=offset

if n_elements(offset) eq 0 then offset=0.2
print,'Exposure = ',exposure
print,'Reference catalog file = ',refcatfile
print,'Offset = ',strtrim(offset,2),' deg'

refcat = MRDFITS(refcatfile,1,/silent)
;; coordinates relative to the center of the field
cendec = mean(minmax(refcat.dec))
if range(refcat.ra) gt 100 then begin
  ra = refcat.ra
  bd = where(ra gt 180,nbd)
  if nbd gt 0 then ra[bd]-=360
  cenra = mean(minmax(ra))
  if cenra lt 0 then cenra+=360
endif else cenra=mean(minmax(refcat.ra))
ROTSPHCEN,refcat.ra,refcat.dec,cenra,cendec,reflon,reflat,/gnomic

field = first_el(strsplit(exposure,'-',/extract))

;; Loop over chips
For c=1,62 do begin
  schip = string(c,format='(i02)')
  chipfile = field+'/chip'+schip+'/'+exposure+'_'+schip+'.fits'
  print,chipfile
  if file_test(chipfile) gt 0 then begin
    ;; Save the reference catalog for this chip
    hd = PHOTRED_READFILE(chipfile,/header)
    nx = sxpar(hd,'naxis1')
    ny = sxpar(hd,'naxis2')
    HEAD_XYAD,hd,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
    ROTSPHCEN,vra,vdec,cenra,cendec,vlon,vlat,/gnomic
    gdrefcat = where(reflon ge min(vlon)-offset and reflon le max(vlon)+offset and $
                     reflat ge min(vlat)-offset and reflat le max(vlat)+offset,ngdrefcat)
    refcat1 = refcat[gdrefcat]
    chiprefcatfile = field+'/chip'+schip+'/'+exposure+'_'+schip+'_refcat.fits'
    MWRFITS,refcat1,chiprefcatfile,/create
    SPAWN,['gzip','-f',chiprefcatfile],/noshell
  endif else print,chipfile,' NOT FOUND'
Endfor  ; chip loop

;stop

end
