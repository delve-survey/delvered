pro remake_delvered_expforced,brick

;; The loadtfr.pro program was messed up.  Remake the expforced.fits catalogs for delvered_forcebrick.pro

print,'Remaking expcat.fits catalog for ',brick

bdir = '/net/dl2/dnidever/delve/bricks/'+strmid(brick,0,4)+'/'+brick+'/'

;; The fast way to correct the OBJIDs
expfile = bdir+brick+'_expforced.fits'
expcat = MRDFITS(expfile+'.gz',1)
;; the number part of ID and OBJID should be the same
;;   ID              STRING    '00619012_19.21677'
;;   OBJID           STRING    '0988m505.21677'
;; But the number in the old OBJID was wrong, e.g.
;;   ID              STRING    '00619012_19.21677'
;;   OBJID           STRING    '0988m505.37821'  37821 is WRONG!!
;; So just take the 21677 from ID and put it in OBJID with BRICK
id = strtrim(expcat.id,2)
num = reform((strsplitter(id,'.',/extract))[1,*])
objid = brick+'.'+num
expcat.objid = objid

;; Saving exposure-level forced photometry catalog
print,'Saving corrected file to ',expfile
MWRFITS,expcat,expfile,/create
spawn,['gzip','-f',expfile],/noshell

return


;; This code below is correct but overkill

cd,current=curdir
cd,bdir

chstr = mrdfits(brick+'_meta.fits',1)
chstr.base = strtrim(chstr.base,2)
chstr.expnum = strtrim(chstr.expnum,2)
chstr.filter = strtrim(chstr.filter,2)
chstr.utdate = strtrim(chstr.utdate,2)
chstr.uttime = strtrim(chstr.uttime,2)
nchstr = n_elements(chstr)
mchfile = file_search('F*_comb.mch',count=nmchfile)
mchbase = file_basename(mchfile,'_comb.mch')

;; Make the exposure-level forced photometry catalog
;;--------------------------------------------------
;; Load the individual ALF files to get chi, sharp
;; Load TFR file to conver ALF IDs to final object ID
tfrfile = mchbase+'_comb.tfr'
LOADTFR,tfrfile,alffiles,tfrstr
schema = {id:'',objid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,$
          imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0}
expcat = replicate(schema,long(total(chstr.alf_nsources))+10000L)
cnt = 0LL
for i=0,nchstr-1 do begin
  base1 = chstr[i].base
  fitsfile = base1+'.fits'
  if file_test(alffiles[i]) eq 1 and file_test(fitsfile) eq 1 then begin
    LOADALS,alffiles[i],alf,count=nalf
    if nalf eq 0 then goto,BOMB2
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
    objid = brick+'.'+strtrim(tfrstr[ind1].id,2)

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

stop

;; Saving exposure-level forced photometry catalog
expfile = bdir+brick+'_expforced.fits'
MWRFITS,expcat,expfile,/create
spawn,['gzip','-f',expfile],/noshell

stop

cd,curdir

end
