pro fix_badobjid,brick

;; Some bricks have a problem where the objids are corrupted, e.g.
;; 0666m472.72.000 and not uniue.
;; It looks like the .makemag file is okay but from .mag file onwards
;; there are problems.

print,'Fixing objids for ',brick

logfile = -1

bdir = '/net/dl2/dnidever/delve/bricks/'+strmid(brick,0,4)+'/'+brick+'/'
cd,current=curdir
cd,bdir

mchfiles = file_search(bdir+'F*_comb.mch',count=nmchfiles)
if nmchfiles eq 0 then begin
  print,'No mch files'
  return
endif
mchbase = file_basename(mchfiles[0],'_comb.mch')
combbase = mchbase+'_comb'


;; Code from end of ALLFRAME.PRO

;; Load the SExtractor file
print,'Loading sex file'
sexfile = combbase+'_allf.sex'
sex = MRDFITS(sexfile,1,/silent)

;; Load the MAKEMAG file
print,'Loading makemag file'
LOADMAKEMAG,mchbase+'.makemag',mag,alfhead
nmag = n_elements(mag)

;; Match them with IDs
MATCH,mag.id,sex.number,ind1,ind2,count=nind

;; Add SExtractor information to mag file
sextags = tag_names(sex)

;; New columns
newcols = ['FLAGS','CLASS_STAR','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA_WORLD',$
           'A_WORLD','B_WORLD','THETA_WORLD','ELLIPTICITY','FWHM_WORLD']
newname = ['FLAG','PROB','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA',$
           'ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM']
for k=0,n_elements(newcols)-1 do begin
  colind = where(sextags eq newcols[k],ncolind)
  if ncolind gt 0 then begin
    add_tag,mag,newname[k],fix('',type=size(sex[0].(colind),/type)),mag
    mag[ind1].(n_tags(mag)-1) = sex[ind2].(colind)
    ;; convert to arcsec
    if newcols[k] eq 'A_WORLD' or newcols[k] eq 'B_WORLD' or newcols[k] eq 'FWHM_WORLD' then mag[ind1].(n_tags(mag)-1) *= 3600
  endif
endfor

;; Write the final output file
finalfile = mchbase+'.mag'
print,'Writing ',finalfile
;; FITS has a limit 999 columns/fields for binary tables, use ASCII                                                                                                                 
;; if over that limit                                                                                                                                                               
catformat = 'FITS'
if catformat eq 'FITS' and n_tags(mag) gt 999 then printlog,logf,'Cannot use FITS output because number of columns>999.  Using ASCII instead'
if (catformat eq 'FITS') and (n_tags(mag) lt 1000) then begin
  MWRFITS,mag,finalfile,/create,/silent
endif else begin  ; ASCII                                                                                                                                                           
  PRINTSTR,mag,finalfile,/silent
endelse



;; Code from end of delvered_forcebrick.pro
magfile = finalfile
mchfile = mchbase+'.mch'
chstr = mrdfits(brick+'_meta.fits',1)
chstr.base = strtrim(chstr.base,2)
nchstr = n_elements(chstr)

delvereddir = '/home/dnidever/projects/delvered/'
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1,/silent)
bind = where(brickstr.brickname eq brick,nbind)
brickstr1 = brickstr[bind[0]]
tilestr = MAKE_BRICK_WCS(brickstr1)

;; Step 4: Calculate coordinates
;;-------------------------------
printlog,logfile,'Step 4: Adding coordinates'
; Load the MCH file
LOADMCH,mchfile,alsfiles
nalsfiles = n_elements(alsfiles)
; Load the photometry file
instphot = PHOTRED_READFILE(magfile)
ninstphot = n_elements(instphot)
printlog,logfile,'Nstars = '+strtrim(ninstphot,2)
;; Converting to IDL X/Y convention, starting at (0,0)
;; DAOPHOT has X/Y start at (1,1)
HEAD_XYAD,tilestr.head,instphot.x-1.0,instphot.y-1.0,ra,dec,/degree


;; Step 5: Calibrating photometry with zero-points
;;-------------------------------------------------
printlog,logfile,'Step 5: Calibrating photometry with zero-points'
cmag = fltarr(ninstphot,nchstr)+99.99
cerr = fltarr(ninstphot,nchstr)+9.99
;; Chip loop
for i=0,nchstr-1 do begin
  ; id, x, y, ra, dec, unsolved magnitudes, chi, sharp
  imag = instphot.(2*i+3)
  ierr = instphot.(2*i+4)
  gdmag = where(imag lt 50,ngdmag)
  if ngdmag gt 0 then begin
    ;; exptime, aperture correction, zero-point
     ;; aperture correction is SUBTRACTIVE, makes it brighter
    ;; ZPTERM is a SUBTRACTIVE constant offset
    cmag[gdmag,i] = imag[gdmag] + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].calib_zpterm
    ;; Add zero-point error in quadrature
    cerr[gdmag,i] = sqrt(ierr[gdmag]^2+chstr[i].calib_zptermsig^2)
  endif
endfor
;; Calculate average photometry per filter
ufilt = chstr[uniq(chstr.filter,sort(chstr.filter))].filter
nufilt = n_elements(ufilt)
avgmag = fltarr(ninstphot,nufilt)
avgerr = fltarr(ninstphot,nufilt)
ndet = lonarr(ninstphot,nufilt)
for i=0,nufilt-1 do begin
  gdf = where(chstr.filter eq ufilt[i],ngdf)
  ;; Single exposure
  if ngdf eq 1 then begin
    avgmag[*,i] = cmag[*,gdf[0]]
    avgerr[*,i] = cerr[*,gdf[0]]
    ndet[*,i] = (cmag[*,gdf[0]] lt 50)
  ;; Multiple exposures
  endif else begin
    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(ninstphot)
    totalfluxwt = dblarr(ninstphot)
    for k=0,ngdf-1 do begin
      gdmag = where(cmag[*,gdf[k]] lt 50,ngdmag)
      if ngdmag gt 0 then begin
        totalwt[gdmag] += 1.0d0/cerr[gdmag,gdf[k]]^2
        totalfluxwt[gdmag] += 2.5118864d^cmag[gdmag,gdf[k]] * (1.0d0/cerr[gdmag,gdf[k]]^2)
        ndet[gdmag,i]++
      endif
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif
    avgmag[*,i] = newmag
    avgerr[*,i] = newerr
  endelse
endfor
;; Measure scatter
scatter = fltarr(ninstphot,nufilt)+99.99
for i=0,nufilt-1 do begin
  gdf = where(chstr.filter eq ufilt[i],ngdf)
  if ngdf gt 1 then begin
    totaldiff = fltarr(ninstphot)
    for k=0,ngdf-1 do begin
      gdmag = where(cmag[*,gdf[k]] lt 50,ngdmag)
      if ngdmag gt 0 then totaldiff[gdmag] += (avgmag[gdmag,i]-cmag[gdmag,gdf[k]])^2
    endfor
    scatter[*,i] = sqrt( totaldiff/(ndet[*,i]>1) )
    bd = where(ndet[*,i] le 1,nbd)
    if nbd gt 0 then scatter[bd,i]=99.99
  endif
endfor
;; Create final catalog schema
newschema = {objid:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0}
;; Add columns for calibrated single-epoch photometry columns
cmagnames = strarr(nchstr)
cerrnames = strarr(nchstr)
for i=0,nufilt-1 do begin
  ind = where(chstr.filter eq ufilt[i],nind)
  cmagnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(nind)+1,2)+'MAG'
  cerrnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(nind)+1,2)+'ERR'
endfor
for i=0,nchstr-1 do newschema = create_struct(newschema,cmagnames[i],0.0,cerrnames[i],0.0)
;; Add columns for average photometry per filter
for i=0,nufilt-1 do newschema = create_struct(newschema,ufilt[i]+'MAG',0.0,ufilt[i]+'ERR',0.0,ufilt[i]+'SCATTER',0.0,'NDET'+ufilt[i],0L)
;; Extra columns
newschema = create_struct(newschema,'chi',0.0,'sharp',0.0,'prob',0.0,'ebv',0.0)
;; other SE columns
newschema = create_struct(newschema,'mag_auto',0.0,'magerr_auto',0.0,'asemi',0.0,'bsemi',0.0,'theta',0.0,'ellipticity',0.0,'fwhm',0.0)
;; in unique brick area
newschema = create_struct(newschema,'brickuniq',0B)
;; Create final catalog
phot = replicate(newschema,ninstphot)
struct_assign,instphot,phot,/nozero
phtags = tag_names(phot)
;; object IDs
phot.objid = brick+'.'+strtrim(instphot.id,2)
;; Stuff in the coordinates calculated above
phot.ra = ra
phot.dec = dec
;; Stuff in the calibrated single-epoch photometry columns
for i=0,nchstr-1 do begin
  magind = where(strupcase(phtags) eq cmagnames[i],nmagind)
  errind = where(strupcase(phtags) eq cerrnames[i],nerrind)
  phot.(magind) = cmag[*,i]
  phot.(errind) = cerr[*,i]  
endfor
;; Stuff in the average photometry per filter
for i=0,nufilt-1 do begin
  magind = where(strupcase(phtags) eq strupcase(ufilt[i])+'MAG',nmagind)
  errind = where(strupcase(phtags) eq strupcase(ufilt[i])+'ERR',nerrind)
  scatind = where(strupcase(phtags) eq strupcase(ufilt[i])+'SCATTER',nscatind)
  detind = where(strupcase(phtags) eq 'NDET'+strupcase(ufilt[i]),ndetind)
  phot.(magind) = avgmag[*,i]
  phot.(errind) = avgerr[*,i]
  phot.(scatind) = scatter[*,i]
  phot.(detind) = ndet[*,i]
endfor

;; Calculate SFD E(B-V)
GLACTC,phot.ra,phot.dec,2000.0,glon,glat,1,/deg
phot.ebv = dust_getval(glon,glat,/noloop,/interp)

;; THIS IS NOW BEING DONE IN DELVERED_FINALCAT.PRO THAT COMBINES ALL CATALOGS
;; Only include objects that are INSIDE the UNIQUE brick area
;;   INCLUSIVE at the lower RA and DEC limit
;; Getting objects that are in the UNIQUE brick area
if brickstr1.dec eq -90 then begin
  ;; the brick right at the pole does not have any RA limits
  ginside = where(phot.dec lt brickstr1.dec2,ninside)
endif else begin
  ginside = where(phot.ra ge brickstr1.ra1 and phot.ra lt brickstr1.ra2 and $
                  phot.dec ge brickstr1.dec1 and phot.dec lt brickstr1.dec2,ninside)
endelse
if ninside gt 0 then phot[ginside].brickuniq=1B

;; Get some meta-data
for i=0,nchstr-1 do begin
  alffile = chstr[i].base+'.alf'
  if file_test(alffile) eq 1 then chstr[i].alf_nsources=file_lines(alffile)-3
endfor


;; Make the exposure-level forced photometry catalog
;;--------------------------------------------------
;; Load the individual ALF files to get chi, sharp
;; Load TFR file to conver ALF IDs to final object ID
print,'Making the exposure-level forced photometry catalog'
tfrfile = mchbase+'_comb.tfr'
LOADTFR,tfrfile,alffiles,tfrstr
schema = {id:'',objid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,$
          imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0}
expcat = replicate(schema,long(total(cmag lt 50))+10000L)
cnt = 0LL
for i=0,nchstr-1 do begin
  base1 = chstr[i].base
  fitsfile = base1+'.fits'
  if file_test(alffiles[i]) eq 1 and file_test(fitsfile) eq 1 then begin
    LOADALS,alffiles[i],alf,count=nalf
    if nalf eq 0 then goto,BOMB2
    ;; Sometimes the rows are duplicated in the ALF file
    ui = uniq(alf.id,sort(alf.id))
    if n_elements(ui) lt nalf then begin
      alf = alf[ui]
      nalf = n_elements(alf)
    endif
    head = photred_readfile(fitsfile,/header)
    
    ;; Calibrate the photometry
    ;; exptime, aperture correction, zero-point
    ;; aperture correction is SUBTRACTIVE, makes it brighter
    ;; ZPTERM is a SUBTRACTIVE constant offset
    cmag1 = alf.mag + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].calib_zpterm
    ;; Add zero-point error in quadrature  
    cerr1 = sqrt(alf.err^2+chstr[i].calib_zptermsig^2)
    
    ;; Coordinates
    HEAD_XYAD,head,alf.x-1,alf.y-1,ra1,dec1,/deg
    
    ;; MATCH up ALF IDs to TFR INDEX
    ;; One row per unique object
    ;; the INDEX values are 1-based indices into the ALF files
    MATCH,tfrstr.index[i],lindgen(nalf)+1,ind1,ind2,/sort,count=nmatch
    objid = brick+'.'+strtrim(long(tfrstr[ind1].id),2)
    
    ;; Create the new catalog
    newcat = replicate(schema,nalf)
    newcat.objid = objid
    newcat.id = chstr[i].expnum+'_'+strtrim(chstr[i].chip,2)+'.'+strtrim(alf.id,2)
    newcat.exposure = chstr[i].base
    newcat.ccdnum = chstr[i].chip
    newcat.filter = chstr[i].filter
    newcat.mjd = date2jd(chstr[i].utdate+'T'+chstr[i].uttime,/mjd)
    newcat.x = alf.x
    newcat.y = alf.y
    newcat.ra = ra1
    newcat.dec = dec1
    newcat.imag = alf.mag
    newcat.ierr = alf.err
    newcat.mag = cmag1
    newcat.err = cerr1
    newcat.sky = alf.sky
    newcat.chi = alf.chi
    newcat.sharp = alf.sharp
    
    ;; Add more elements
    if cnt+nalf gt n_elements(expcat) then expcat=add_elements(expcat,100000L>nalf)
    
    ;; Add to global catalog
    expcat[cnt:cnt+nalf-1] = newcat
    cnt += nalf

    BOMB2:
endif else print,alffile+' NOT FOUND'
endfor
expcat = expcat[0:cnt-1]  ; this should not be needed                                                                                                                                 

;; Object catalog
;;---------------
lo = first_el(where(phtags eq 'DEC',nlo))
hi = first_el(where(strupcase(phtags) eq strupcase(ufilt[0])+'MAG',nhi))
obj_schema = create_struct(phtags[0],fix('',type=size(phot.(0),/type)))
for i=1,lo do obj_schema = create_struct(obj_schema,phtags[i],fix('',type=size(phot.(i),/type)))
for i=hi,n_elements(phtags)-2 do obj_schema = create_struct(obj_schema,phtags[i],fix('',type=size(phot.(i),/type)))
obj_schema = create_struct(obj_schema,'brickuniq',0B)
obj = replicate(obj_schema,ninstphot)
STRUCT_ASSIGN,phot,obj,/nozero

;; Saving final catalog
photfile = bdir+brick+'.fits'
if n_tags(phot) le 999 then begin
  printlog,logfile,'Writing photometry to '+photfile+'.gz'
  MWRFITS,phot,photfile,/create
  spawn,['gzip','-f',photfile],/noshell
endif else begin
  printlog,logfile,'Too many columns for FITS.  Saving as IDL SAVE file instead. '+bdir+brick+'.dat'
  SAVE,phot,file=bdir+brick+'.dat'
endelse

;; Saving object photometry catalog
objfile = bdir+brick+'_object.fits'
MWRFITS,obj,objfile,/create
spawn,['gzip','-f',objfile],/noshell

;; Saving exposure-level forced photometry catalog
expfile = bdir+brick+'_expforced.fits'
MWRFITS,expcat,expfile,/create
spawn,['gzip','-f',expfile],/noshell

;; Save metadata
metafile = bdir+brick+'_meta.fits'
printlog,logfile,'Writing meta-data to '+metafile
MWRFITS,chstr,metafile,/create

cd,curdir

end
