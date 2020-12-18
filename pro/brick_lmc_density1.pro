pro brick_lmc_density1,brick,redo=redo

;; Calculate the density of LMC stars in this brick

delvedir = '/net/dl2/dnidever/delve/'
brickdir = delvedir+'bricks/'

nbrick = n_elements(brick)
;; Multiple bricks input
if nbrick gt 1 then begin
  for i=0,nbrick-1 do brick_lmc_density1,brick[i],redo=redo
  return
endif

subdir = brickdir+strmid(brick,0,4)+'/'
bdir = subdir+brick+'/'
objfile = bdir+brick+'_object.fits.gz'
if file_test(objfile) eq 0 then begin
  print,objfile+' NOT FOUND'
  return
endif

outfile = brickdir+'summary/lmc1/'+brick+'.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile+' EXISTS already and /redo NOT set'
  return
endif

;; Load the object file
obj = mrdfits(objfile,1,/silent)
nobj = n_elements(obj)
print,brick,' ',strtrim(nobj,2)
tags = tag_names(obj)

;; Load the meta file
meta = mrdfits(bdir+brick+'_meta.fits',1,/silent)

brkstr = {brick:brick,ra:0.0d0,dec:0.0d0,area:0.0,gdepth95:99.99,rdepth95:99.99,idepth95:99.99,$
          gmin:99.99,rmin:99.99,imin:99.99,$
          gdepth5:99.99,rdepth5:99.99,idepth5:99.99,$
          nobj:0L,nblue:-1L,nlmc1:-1L,nlmc2:-1L,nlmc220:-1L,nlmc225:-1L,nlmc230:-1L,nlmc235:-1L,nback1:-1L,nback2:-1L}
brktags = tag_names(brkstr)
cenra = mean(minmax(obj.ra))
if range(obj.ra) gt 100 then begin
  b = where(obj.ra gt 180,nb)
  if nb gt 0 then obj[b].ra-=360
  cenra = mean(minmax(obj.ra))
  if cenra lt 0 then cenra+=360
endif
cendec = mean(minmax(obj.dec))
brkstr.ra = cenra
brkstr.dec = cendec
brkstr.nobj = n_elements(obj)

;; Check that we have the bands
gmagind = where(tags eq 'GMAG',ngmagind)
imagind = where(tags eq 'IMAG',nimagind)
if ngmagind eq 0 then begin
  print,'NO g-band information'
  return
endif
if nimagind eq 0 then begin
  print,'NO i-band information'
  return
endif

  ;; Extinction coefficients
  ;;u  4.239
  ;;g  3.303
  ;;r  2.285
  ;;i  1.698
  ;;z  1.263

gd = where(obj.gmag lt 50 and obj.imag lt 50 and abs(obj.sharp) lt 1 and obj.chi lt 3 and obj.prob gt 0.5,ngd)
if ngd gt 0 then begin
  gmag = obj[gd].gmag-obj[gd].ebv*3.303
  gi = (obj[gd].gmag-obj[gd].ebv*3.303)-(obj[gd].imag-obj[gd].ebv*1.698)

  ;; lmc density, maybe with various mag ranges
  gdlmc1 = where(gi ge -0.10 and gi le 0.56 and gmag gt 21.8 and gmag lt 22.8,ngdlmc1)
  brkstr.nlmc1 = ngdlmc1
  ;; 22<g<24
  gdlmc2 = where(gi ge -0.10 and gi le 0.56 and gmag gt 22.0 and gmag lt 24.0,ngdlmc2)
  brkstr.nlmc2 = ngdlmc2

  ;; 22.0<g<22.5
  gdlmc220 = where(gi ge -0.10 and gi le 0.50 and gmag gt 22.0 and gmag lt 22.5,ngdlmc220)
  brkstr.nlmc220 = ngdlmc220
  ;; 22.5<g<23.0
  gdlmc225 = where(gi ge 0.0 and gi le 0.50 and gmag gt 22.5 and gmag lt 23.0,ngdlmc225)
  brkstr.nlmc225 = ngdlmc225
  ;; 23.0<g<23.5
  gdlmc230 = where(gi ge 0.0 and gi le 0.60 and gmag gt 23.0 and gmag lt 23.5,ngdlmc230)
  brkstr.nlmc230 = ngdlmc230
  ;; 23.5<g<24.0
  gdlmc235 = where(gi ge 0.0 and gi le 0.70 and gmag gt 23.05and gmag lt 24.0,ngdlmc235)
  brkstr.nlmc235 = ngdlmc235

  ;; blue stars
  gdblue = where(gi ge -1 and gi lt 0.0 and gmag gt 18 and gmag lt 21.8,ngdblue)
  brkstr.nblue = ngdblue

  ;; density in various background regions
  ;; MW halo
  gdback1 = where(gi ge 0.60 and gi le 0.80 and gmag ge 22.5 and gmag le 23.0,ngdback1)
  brkstr.nback1 = ngdback1

  ;; MW disk
  gdback2 = where(gi ge 0.40 and gi le 0.60 and gmag ge 20.0 and gmag le 20.5,ngdback2)
  brkstr.nback2 = ngdback2
endif

;; depths in each band
bands = ['g','r','i']
for i=0,n_elements(bands)-1 do begin
  magind = where(tags eq strupcase(bands[i])+'MAG',nmagind)
  minind = where(brktags eq strupcase(bands[i])+'MIN',nminind)
  depth5ind = where(brktags eq strupcase(bands[i])+'DEPTH5',ndepth5ind)
  depth95ind = where(brktags eq strupcase(bands[i])+'DEPTH95',ndepth95ind)
  gdmag = where(obj.(magind[0]) lt 50,ngdmag)
  if ngdmag gt 0 then brkstr.(minind[0]) = min(obj[gdmag].(magind[0]))
  if ngdmag gt 5 then begin
    ;; Get 95% percentile depth
    mag = obj[gdmag].(magind[0])
    si = sort(mag)
    mag = mag[si]
    depth95 = mag[round(0.95*ngdmag)-1]
    brkstr.(depth95ind[0]) = depth95
    ;; Get 5% percentile depth
    depth5 = mag[round(0.05*ngdmag)-1]
    brkstr.(depth5ind[0]) = depth5
  endif
endfor

;; unique area
rotsphcen,obj.ra,obj.dec,cenra,cendec,lon,lat,/gnomic
undefine,dum,im
dx = 0.005
hess,lon,lat,dum,im,dx=dx,dy=dx,/noplot
npix = total(im gt 0)
area = npix * dx^2
brkstr.area = area

;; Save
print,'Writing to ',outfile
MWRFITS,brkstr,outfile,/create

end
