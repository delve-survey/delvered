function load_chipcat,chstr,usealf=usealf

;; Load the single ALS file and calibrate it the information in CHSTR

;; Not enough inputs
if n_elements(chstr) eq 0 then begin
  print,'Syntax - cat=load_chipcat(chstr)'
  return,-1
endif

;; Check that we have all of the necessary columns
tags = tag_names(chstr)
cols = ['FILE','EXPTIME','APCOR','CALIB_ZPTERM','CALIB_ZPTERMSIG',$
        'BASE','FILTER','EXPNUM','CHIP','UTDATE','UTTIME']
for i=0,n_elements(cols)-1 do begin
  tind = where(tags eq cols[i],ntind)
  if ntind eq 0 then begin
    print,'CHSTR must have '+cols[i]+' column'
    return,-1
  endif
endfor

fitsfile = strtrim(chstr.file,2)

;; Load ALS file
if not keyword_set(usealf) then begin
  alsfile = repstr(fitsfile,'.fits','.als')
endif else begin
  ;; Using ALF file
  alsfile = chstr.alffile
endelse
if file_test(alsfile) eq 0 then begin
  print,alsfile+' NOT FOUND'
  return,-1
endif
LOADALS,alsfile,als,count=nals
if nals eq 0 then begin
  print,alsfile+' IS EMPTY'
  return,-1
endif
;; Load FITS header
head = PHOTRED_READFILE(fitsfile,/header)

;; Trim out bad half of DECam ccdnum 31
dateobs = chstr.utdate+'T'+chstr.uttime
mjd = date2jd(dateobs,/mjd)
if chstr.chip eq 31 and mjd gt 56660 then begin
  print,'Masking bad half of DECam chip 31'
  ;; X: 1-1000 okay
  ;; X: 1000-2049 bad
  bdind = where(als.x ge 1000,nbdind,comp=gdind,ncomp=ngdind)
  if nbdind gt 0 then begin   ; some bad ones found
    if ngdind eq 0 then begin   ; all bad
      print,'NO useful measurements in ',fitsfile
      return,-1
    endif else begin
      print,'Removing '+strtrim(nbdind,2)+' bad measurements, '+strtrim(ngdind,2)+' left.'
      REMOVE,bdind,als
      nals = n_elements(als)
    endelse
  endif  ; some bad ones to remove
endif

;; Calibrate the photometry
;; exptime, aperture correction, zero-point
;; aperture correction is SUBTRACTIVE, makes it brighter
;; ZPTERM is a SUBTRACTIVE constant offset
cmag = als.mag + 2.5*alog10(chstr.exptime) - chstr.apcor - chstr.calib_zpterm
;; Add zero-point error in quadrature
cerr = sqrt(als.err^2+chstr.calib_zptermsig^2)

;; Coordinates
HEAD_XYAD,head,als.x-1,als.y-1,ra,dec,/deg

;; Create the new catalog
schema = {id:'',objid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,forced:0B,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,$
          imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0}
cat = replicate(schema,nals)
cat.id = strtrim(chstr.expnum,2)+'_'+strtrim(chstr.chip,2)+'.'+strtrim(als.id,2)
cat.exposure = strtrim(chstr.base,2)
cat.ccdnum = chstr.chip
cat.filter = strtrim(chstr.filter,2)
cat.mjd = date2jd(strtrim(chstr.utdate,2)+'T'+strtrim(chstr.uttime,2),/mjd)
cat.x = als.x
cat.y = als.y
cat.ra = ra
cat.dec = dec
cat.imag = als.mag
cat.ierr = als.err
cat.mag = cmag
cat.err = cerr
cat.sky = als.sky
cat.chi = als.chi
cat.sharp = als.sharp

return,cat

end
