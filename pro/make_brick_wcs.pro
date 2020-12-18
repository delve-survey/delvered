function make_brick_wcs,brickstr

;; This creates a brick WCS given the brick structure

;; Each brick covers 0.25 x 0.25 deg at a pixel scale of 0.262"/pix
;; with 3600x3600 pixels that gives us an overlap of ~82 pixels on
;; each side.  The *UNIQUE* brick area (0.25x0.25 deg) is defined by
;; the BRAMIN/MAX and BDECMIN/MAX keywords.

if n_elements(brickstr) eq 0 then begin
  print,'tilestr = make_brick_wcs(brickstr)'
  return,-1
endif

; Make the tiling file
;---------------------
undefine,lines
; Lines with the tiling scheme first
nx = 3600
ny = 3600
step = 0.262 / 3600    ; 0.262" per pixel, DECam pixel scale
xref = nx/2
yref = ny/2

;;  Make the header as well
MKHDR,tilehead,fltarr(5,5)
SXADDPAR,tilehead,'NAXIS1',nx
SXADDPAR,tilehead,'CDELT1',step
SXADDPAR,tilehead,'CRPIX1',xref+1L
SXADDPAR,tilehead,'CRVAL1',brickstr.ra
SXADDPAR,tilehead,'CTYPE1','RA---TAN'
SXADDPAR,tilehead,'NAXIS2',ny
SXADDPAR,tilehead,'CDELT2',step
SXADDPAR,tilehead,'CRPIX2',yref+1L
SXADDPAR,tilehead,'CRVAL2',brickstr.dec
SXADDPAR,tilehead,'CTYPE2','DEC--TAN'
SXADDPAR,tilehead,'BRAMIN',brickstr.ra1,'RA min of unique brick area'
SXADDPAR,tilehead,'BRAMAX',brickstr.ra2,'RA max of unique brick area'
SXADDPAR,tilehead,'BDECMIN',brickstr.dec1,'DEC min of unique brick area'
SXADDPAR,tilehead,'BDECMAX',brickstr.dec2,'DEC max of unique brick area'
EXTAST,tilehead,tileast
tileast.equinox = 2000

; Create the TILE structure
tilestr = {type:'WCS',naxis:long([nx,ny]),cdelt:double([step,step]),crpix:double([xref+1L,yref+1L]),$
           crval:double([brickstr.ra,brickstr.dec]),ctype:['RA--TAN','DEC--TAN'],$
           head:tilehead,ast:tileast,xrange:[0,nx-1],yrange:[0,ny-1],nx:nx,ny:ny}

return, tilestr

end
