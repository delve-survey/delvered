pro delvered_initpsfstars,file

;; This creates an initial list of PSF stars using Gaia DR2

if n_elements(file) eq 0 then begin
  print,'Syntax - delvered_initpsfstars,file'
  return
endif

print,'Creating initial PSF stars using Gaia DR2 for ',file 

dir = file_dirname(file)
if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then $
  base=FILE_BASENAME(file,'.fits.fz') else $
  base=FILE_BASENAME(file,'.fits')

;; Check that the FITS file exists
if file_test(file) eq 0 then begin
  print,file,' NOT FOUND'
  return
endif
;; Load the header
head = PHOTRED_READFILE(file,/header)
nx = sxpar(head,'naxis1')
ny = sxpar(head,'naxis2')

; Check if the reference file exists
reffile = dir+'/'+base+'_refcat.fits'
if file_test(reffile) eq 0 then reffile+='.gz'
if file_test(reffile) eq 0 then begin
  print,reffile,' NOT FOUND'
  return
endif
refcat = mrdfits(reffile,1,/silent)

;; Get the coordinates in X/Y
HEAD_ADXY,head,refcat.ra,refcat.dec,xref,yref,/deg
gd = where(xref ge 0 and xref le (nx-1) and yref gt 0 and yref le (ny-1) and refcat.gmag lt 50,ngd)
if ngd eq 0 then begin
  print,'No good Gaia DR2 PSF stars'
  return
endif
print,strtrim(ngd,2),' good stars found'

;; Write the output file
print,'Writing initial list of PSF stars to ',dir+'/'+base+'.cmn.lst'
WRITECOL,dir+'/'+base+'.cmn.lst',indgen(ngd)+1,xref[gd],yref[gd],refcat[gd].gmag,fltarr(ngd),fltarr(ngd),fmt='(I7,5f9.3)'
head =  ['NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD',$
         ' 3  '+strtrim(nx,2)+'  '+strtrim(ny,2)+'   361.3 180395.   92.43    3.00    1.03    5.96    3.89','']
WRITELINE,dir+'/'+base+'.cmn.lst',head,/prepend

;stop

end
