pro brick_lmc_density2,brick,redo=redo

;; Calculate the density of LMC stars in this brick

delvedir = '/net/dl2/dnidever/delve/'
brickdir = delvedir+'bricks/'
subdir = brickdir+strmid(brick,0,4)+'/'
bdir = subdir+brick+'/'
;objfile = bdir+brick+'_object.fits.gz'
objfile = bdir+brick+'_joint_object.fits.gz'
if file_test(objfile) eq 0 then begin
  print,objfile+' NOT FOUND'
  return
endif

outfile = brickdir+'summary/lmc2/'+brick+'.fits'
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
          gdepth10:99.99,rdepth10:99.99,idepth10:99.99,nobj:0L,nlmc:-1L,nback:-1L}
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

;; Sometimes these coordinates can be WAY off
;; Use the ones from the master brick list


  ;; Extinction coefficients
  ;;u  4.239
  ;;g  3.303
  ;;r  2.285
  ;;i  1.698
  ;;;z  1.263

gd = where(obj.gmag lt 50 and obj.rmag lt 50 and abs(obj.sharp) lt 1 and obj.chi lt 3 and obj.prob gt 0.5,ngd)
if ngd gt 0 then begin
  rmag = obj[gd].rmag-obj[gd].ebv*2.285
  gr = (obj[gd].gmag-obj[gd].ebv*3.303)-(obj[gd].rmag-obj[gd].ebv*2.285)

  ;; lmc density, maybe with various mag ranges
  gdlmc = where(gr ge 0.10 and gr le 0.40 and rmag gt 21.8 and rmag lt 22.7,ngdlmc)
  brkstr.nlmc = ngdlmc

  ;; MW disk
  gdback = where(gr ge 0.20 and gr le 0.45 and rmag ge 19.0 and rmag le 21.2,ngdback)
  brkstr.nback = ngdback
endif

;; depths in each band
bands = ['g','r','i']
for i=0,n_elements(bands)-1 do begin
  magind = where(tags eq strupcase(bands[i])+'MAG',nmagind)
  errind = where(tags eq strupcase(bands[i])+'ERR',nerrind)
  depthind = where(brktags eq strupcase(bands[i])+'DEPTH95',ndepthind)
  depth10ind = where(brktags eq strupcase(bands[i])+'DEPTH10',ndepth10ind)
  gdmag = where(obj.(magind[0]) lt 50,ngdmag)
  if ngdmag gt 5 then begin
    ;; Get 95% percentile depth
    mag = obj[gdmag].(magind[0])
    si = sort(mag)
    mag = mag[si]
    depth95 = mag[round(0.95*ngdmag)-1]
    brkstr.(depthind[0]) = depth95
    ;; Get 10 sigma depth
    ;;  S/N = 1.087/err
    ;;  so S/N=5 is for err=1.087/5=0.2174
    ;;  S/N=10 is for err=1.087/10=0.1087
    ;depth10sig = 99.99
    depind = where(obj.(magind[0]) lt 50 and obj.(magind[0]) gt depth95-3.0 and obj.(errind[0]) ge 0.0987 and obj.(errind[0]) le 0.1187,ndepind)
    if ndepind lt 5 then depind = where(obj.(magind[0]) lt 50 and obj.(magind[0]) gt depth95-3.0 and obj.(errind[0]) ge 0.0787 and obj.(errind[0]) le 0.1387,ndepind)
    if ndepind gt 5 then begin
      depth10sig = median([obj[depind].(magind[0])])
    endif else begin
      depind = where(obj.(magind[0]) lt 50,ndepind)
      if ndepind gt 0 then depth10sig=max([obj[depind].(magind[0])])
    endelse    
    brkstr.(depth10ind[0]) = depth10sig
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
