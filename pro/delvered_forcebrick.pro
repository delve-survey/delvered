;+
;
; DELVERED_FORCEBRICK
;
; Process a single DELVE brick and perform ALLFRAME FORCED photometry
;
; INPUTS:
;  brick   The DELVE brick name, e.g. 1234m045
;
; By D. Nidever  August 2019
;-

pro delvered_forcebrick,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,update=update,$
                        logfile=logfile,delvedir=delvedir

;; This bricks pre-processing script gets DELVE and community MC data ready
;; to run PHOTRED ALLFRAME on it.

t0 = systime(1)
CD,current=curdir

;; Not enough inputs
if n_elements(brick) eq 0 then begin
  print,'Syntax - delvered_forcebrick,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,logfile=logfile'
  return
endif

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
if FILE_TEST(delvedir,/directory) eq 0 then FILE_MKDIR,delvedir
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
;expfile = '/home/dnidever/projects/delvered/data/decam_mcs_20181009.fits.gz'
if n_elements(scriptsdir) gt 0 then scriptsdir=trailingslash(scriptsdir) else scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
if n_elements(irafdir) gt 0 then irafdir=trailingslash(irafdir) else irafdir='/home/dnidever/iraf/'
tempdir = '/tmp/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir = '/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
if file_test(brickdir,/directory) eq 0 then file_mkdir,brickdir
logsdir = brickdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
if n_elements(logfile) eq 0 then logfile=-1

; Start the logfile
;------------------
; format is delvered_forcebrick.brick.host.DATETIME.log
host = GETENV('HOST')
hostname = first_el(strsplit(host,'.',/extract))
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logtime = smonth+sday+syear+shour+sminute+ssecond
journalfile = logsdir+'delvered_brick.'+brick+'.'+hostname+'.'+logtime+'.log'
JOURNAL,journalfile

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1,/silent)

;; Get the brick information
bind = where(brickstr.brickname eq brick,nbind)
if nbind eq 0 then begin
  printlog,logfile,ibrick+' not in DELVE-MC brick list'
  return
endif
brickstr1 = brickstr[bind[0]]

;; Subdirectory is the first 4 digits of the brickname, e.g., 0952 of 0952m462, the RA portion of the name
subdir = brickdir+strmid(brick,0,4)+'/'
if file_test(subdir,/directory) eq 0 then file_mkdir,subdir
;; Create brick directory if it doesn't exist yet
bdir = subdir+brick+'/'
if file_test(bdir,/directory) eq 0 then file_mkdir,bdir
logfile = bdir+brick+'.'+logtime+'.log'

;; Check the output file
photfile = bdir+brick+'.fits'
if file_test(photfile+'.gz') and not keyword_set(redo) and not keyword_set(update) then begin
  printlog,logfile,photfile+'.gz EXISTS and /redo NOT set'
  return
endif

;; DECam imager
thisimager = {telescope:'BLANCO',instrument:'DECam',namps:62,separator:'_'}

;; Print information
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,' Run ALLFRAME forced photometry on DELVE Brick = ',brick
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,'DELVEDIR = ',delvedir
printlog,logfile,'SCRIPTSDIR = ',scriptsdir
printlog,logfile,'IRAFDIR = ',irafdir
printlog,logfile,'EXPOSUREDIR = ',expdir
printlog,logfile,'WORKINGDIR = ',workdir
printlog,logfile,'HOST = ',host



;; The photed.setup file
setup = ['##### REQUIRED #####',$
         'scriptsdir  '+scriptsdir,$
         'irafdir     '+irafdir,$
         'telescope   Blanco',$
         'instrument  DECAM',$
         'observatory CTIO',$
         'nmulti      10',$
         'nmulti_wcs       40',$
         'nmulti_daophot   30',$
         'nmulti_allframe  10',$
         'filtref     g,i,r,z,u',$
         'trans       delve.trans',$
         '##### OPTIONAL #####',$
         'sepfielddir  1',$
         'sepchipdir   1',$
         'keepmef      0',$
         'catformat    FITS',$
         'workdir      '+workdir,$
         'clean        1',$
         'skipcheck    1',$
         'redo         0',$
         'wcsrefname   GAIADR2',$
         'searchdist   20',$
         '#wcsrmslim   1.0',$
         'hyperthread   1',$
         'daopsfva      1',$
         'daofitradfwhm 1.0',$
         'psfcomsrc     0',$
         'psfcomglobal  0',$
         'psfcomgauss   0',$
         '#mchmaxshift  50.0',$
         'finditer      1',$      ; 2->1 on 7/22/20
         'alfdetprog  sextractor',$
         '#alfnocmbimscale 0',$
         'alftrimcomb   0',$
         '#ddo51radoffset  1',$
         'cmbforce      1',$
         'keepinstr     1',$
         'avgmag        1',$
         'avgonlymag    0',$
         'todered       u,g,r,i,z,g-i',$
         '##### STAGES #####',$
         '#rename',$
         '#split',$
         '#wcs',$
         '#daophot',$
         '#zeropoint',$
         '#match',$
         ' allframe',$
         '#apcor',$
         '#astrom',$
         '#calib',$
         '#combine',$
         '#deredden',$
         '#save',$
         '#html']
WRITELINE,bdir+'photred.setup',setup
PHOTRED_LOADSETUP,setup,setupdir=bdir

;; Print out some information about the brick
printlog,logfile,'RA = ',stringize(brickstr1.ra,ndec=5)
printlog,logfile,'DEC = ',stringize(brickstr1.dec,ndec=5)
printlog,logfile,'RA range  = [ ',stringize(brickstr1.ra1,ndec=5),',',stringize(brickstr1.ra2,ndec=5),' ]'
printlog,logfile,'DEC range = [ ',stringize(brickstr1.dec1,ndec=5),',',stringize(brickstr1.dec2,ndec=5),' ]'

;; Get the Brick WCS information
tilestr = MAKE_BRICK_WCS(brickstr1)


;; Step 1: Get the list of exposures/chips that overlap this brick
;;----------------------------------------------------------------
printlog,logfile,'Step 1: Get list of chip files that overlap this brick'
;printlog,logfile,'Getting chips that overlap this brick'
cenra = brickstr1.ra
cendec = brickstr1.dec
tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
;; /noshell causes problems on gp09 because it gives python2 instead of python3
spawn,delvereddir+'bin/query_delvered_summary_table '+strtrim(cenra,2)+' '+strtrim(cendec,2)+' '+tmpfile+' --lim 0.5',out,errout
info = file_info(tmpfile)
if info.size eq 0 then begin
  printlog,logfile,'No overlapping chips found'
  return
endif
chstr = MRDFITS(tmpfile,1,/silent)
file_delete,tmpfile,/allow
nchstr = n_elements(chstr)
;chstr.file = repstr(chstr.file,'/net/dl1/','/dl1/')   ;; fix /net/dl1 to /dl1
printlog,logfile,'Found ',strtrim(nchstr,2),' overlapping chips within 0.5 deg of brick center'

;; Make sure the chips are unique,  some were duplicated on SMASH nights
chid = chstr.expnum+'-'+strtrim(chstr.chip,2)
ui = uniq(chid,sort(chid))
chstr = chstr[ui]
nchstr = n_elements(chstr)

;; Do more rigorous overlap checking
;;  the brick region with overlap
print,'Performing more rigorous overlap checking'
HEAD_XYAD,tilestr.head,[0,tilestr.nx-1,tilestr.nx-1,0],[0,0,tilestr.ny-1,tilestr.ny-1],bvra,bvdec,/deg
olap = intarr(nchstr)
vxarr = fltarr(nchstr,4)
vyarr = fltarr(nchstr,4)
for i=0,nchstr-1 do begin
  if (i mod 100 eq 0) and (i gt 0) then print,i
  hd1 = PHOTRED_READFILE(chstr[i].file,/header)
  nx = sxpar(hd1,'naxis1')
  ny = sxpar(hd1,'naxis2')
  head_xyad,hd1,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
  olap[i] = dopolygonsoverlap(bvra,bvdec,vra,vdec)
  head_adxy,tilestr.head,vra,vdec,vx,vy,/deg
  vxarr[i,*] = vx
  vyarr[i,*] = vy
endfor
;; Require at least a 2 pixel overlap in X and Y
g = where(olap eq 1 and max(vxarr,dim=2) ge 2 and max(vyarr,dim=2) ge 2 and min(vxarr,dim=2) le tilestr.nx-3 and min(vyarr,dim=2) le tilestr.ny-3,ng)
if ng eq 0 then begin
  printlog,logfile,'No chips overlap this brick'
  return
endif
printlog,logfile,strtrim(ng,2),' chips overlap this brick'
chstr = chstr[g]
nchstr = ng

;; APPLY cuts on the exposures
;;-----------------------------
printlog,logfile,'Applying quality, zero-point, filter and exptime cuts'

; Zero-point structure, from NSC
zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},7)
zpstr.instrument = 'c4d'
zpstr.filter = ['u','g','r','i','z','Y','VR']
zpstr[0].amcoef = [-1.60273, -0.375253]   ; c4d-u
zpstr[1].amcoef = [0.277124, -0.198037]   ; c4d-g
zpstr[2].amcoef = [0.516382, -0.115443]   ; c4d-r
zpstr[3].amcoef = [0.380338, -0.067439]   ; c4d-i
zpstr[4].amcoef = [0.074517, -0.067031]   ; c4d-z
zpstr[5].amcoef = [-1.07800, -0.060014]   ; c4d-Y  
zpstr[6].amcoef = [1.111859, -0.083630]   ; c4d-VR

;; Convert to additive zero-point as used in NSC
zpterm = -chstr.calib_zpterm
;; Fix early DES exposures that used different units/gain
gdes = where(chstr.gain lt 2,ngdes)
if ngdes gt 0 then zpterm[gdes] -= 1.55
;; global zpterm and airmass correction
for i=0,n_elements(zpstr)-1 do begin
  ind = where(chstr.filter eq zpstr[i].filter,nind)
  if nind gt 0 then zpterm[ind] -= poly(chstr[ind].airmass,zpstr[i].amcoef)
endfor

fwhmthresh = 2.0                ; seeing 2.0" threshold
filt = strmid(chstr.filter,0,1)
gdch = where(chstr.fwhm*chstr.pixscale le fwhmthresh and chstr.exptime ge 90. and zpterm ge -0.5 and $
             finite(zpterm) eq 1 and finite(chstr.apcor) eq 1 and $
             (filt eq 'u' or filt eq 'g' or filt eq 'r' or filt eq 'i' or filt eq 'z' or filt eq 'Y'),ngdch)

if ngdch eq 0 then begin
  printlog,logfile,'No chips passed the cuts'
  return
endif
printlog,logfile,strtrim(ngdch,2),' chips passed the cuts'
chstr = chstr[gdch]
nchstr = ngdch
chstr.file = strtrim(chstr.file,2)
chstr.base = strtrim(chstr.base,2)
;add_tag,chstr,'newbase','',chstr


;; Check if we need to update
if keyword_set(update) and not keyword_set(redo) then begin
  ;; Load previous meta
  metafile = bdir+brick+'_meta.fits'
  if file_test(metafile) then begin
    meta0 = mrdfits(metafile,1,/silent)
    meta0.base = strtrim(meta0.base,2)
    MATCH,meta0.base,chstr.base,ind1,ind2,/sort,count=nmatch
    if nchstr eq n_elements(meta0) and nchstr eq nmatch then begin
      printlog,logfile,'Nothing to UPDATE'
      return
    endif
    ;; Moving previous catalogs to a backup
    oldfiles = file_search(bdir+[brick+'*fits*','allframe.opt','default.*'],count=noldfiles)
    if noldfiles gt 0 then begin
      bakdir = bdir+'bak'+smonth+sday+syear+shour+sminute+ssecond
      printlog,logfile,'Backing up old files to ',bakdir
      FILE_MKDIR,bakdir
      FILE_MOVE,oldfiles,bakdir
   endif else printlog,logfile,'No old files to backup'
  endif
  printlog,logfile,'There are new exposures to include.  UPDATING'
endif



;; Create temporary local directory to perform the work/processing
;;  copy everything to it
if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
procdir = first_el(MKTEMP('brk',outdir=workdir,/directory))+'/'
procdir = repstr(procdir,'//','/')
FILE_CHMOD,procdir,/a_execute
printlog,logfile,'Working in temporary directory '+procdir


;; Copy the files
;;----------------
printlog,logfile,'Copying over the necessary files to ',procdir
;; Chip loop
for i=0,nchstr-1 do begin
  printlog,logfile,'Copying over files for ',chstr[i].file
  ;; New filename
  ;chstr[i].newbase = chstr[i].base
  odir = file_dirname(chstr[i].file)+'/'
  obase = PHOTRED_GETFITSEXT(chstr[i].file,/basename)
  obase = first_el(obase)
  ;; Need symlinks to .psf, .als
  if file_test(odir+obase+'.psf') eq 0 then begin
    printlog,logfile,odir+obase+'.psf NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  if file_test(odir+obase+'.als') eq 0 then begin
    printlog,logfile,odir+obase+'.als NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_DELETE,procdir+chstr[i].base+['.psf','.als','.ap','.opt','.als.opt','.log'],/allow
  FILE_COPY,odir+obase+['.psf','.als','.ap','.opt','.als.opt','.log'],procdir,/allow_same,/overwrite
  FILE_CHMOD,procdir+chstr[i].base+['.psf','.als','.ap','.opt','.als.opt','.log'],'755'o   ; make sure they are writable
  ;; Copy the fits, fits resource file and header files locally
  if file_test(odir+obase+'.fits') eq 0 then begin
    printlog,logfile,odir+obase+'.fits NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_COPY,odir+obase+'.fits',procdir,/over
  if file_test(odir+'.'+obase+'.fits') eq 1 then FILE_COPY,odir+'.'+obase+'.fits',procdir,/over
  if file_test(odir+obase+'.fits.head') eq 1 then FILE_COPY,odir+obase+'.fits.head',procdir,/over
  FILE_CHMOD,procdir+chstr[i].base+'.fits','755'o
  BOMB1:
endfor


;; Step 2: Run DAOMATCH_TILE.PRO on the files
;;--------------------------------------------
printlog,logfile,'Step 2: Matching up objects with DAOMATCH_TILE'
cd,procdir
groupstr = {x0:0,y0:0}
DAOMATCH_TILE,chstr.base+'.als',tilestr,groupstr


;; Step 3: Run ALLFRAME
;;----------------------
printlog,logfile,'Step 3: Run ALLFRAME'
;; DO I NEED TO HAVE IT TRIM THE COMBINED IMAGE???
mchbase = procdir+chstr[0].base
mchfile = mchbase+'.mch'  ;; allframe needs absolute path
ALLFRAME,mchfile,tile=tilestr,setupdir=bdir,scriptsdir=scriptsdir,irafdir=irafdir,$
         logfile=logfile,catformat='FITS',imager=thisimager,geocoef=0
magfile = chstr[0].base+'.mag'
if file_test(magfile) eq 0 then begin
  printlog,logfile,magfile+' NOT FOUND'
  return
endif

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


;; Delete files we don't want to keep
;;-----------------------------------
;; Individual fits files
alsbase = file_basename(alsfiles,'.als')
;; if fits files have resources files then replace the fits file by a
;;  dummy file
for i=0,n_elements(alsbase)-1 do begin
  info1 = file_info(alsbase[i]+'.fits')
  if info1.exists and info1.size gt 1 then begin
    file_delete,alsbase[i]+'.fits'
    if file_test('.'+alsbase[i]+'.fits') eq 1 then WRITELINE,alsbase[i]+'.fits',''
  endif
endfor
;; Combined files
;;   _comb  lst, lst1, lst2, lst1.chi, grp, nst, lst2.chi, plst.chi, psfini.ap
;;   nei, als.inp, a.fits, cmn.log, cmn.coo, cmn.ap, cmn.lst,
;;   _sub.fits, _sub.cat, _sub.als, _all.coo, makemag
base = file_basename(mchfile,'.mch')
FILE_DELETE,base+'_comb'+['.lst','.lst1','.lst2','.lst1.chi','.lst2.chi','.grp','.nst','.plst.chi',$
                          '.nei','.als.inp','.cmn.log','.cmn.coo','.cmn.ap','.cmn.lst','a.fits',$
                          'a.fits.fz','_sub.fits','_sub.cat','_sub.als','_all.coo','.makemag'],/allow
FILE_DELETE,'check.fits',/allow

;; fpack _comb.fits and _combs.fits
FILE_DELETE,base+['_comb.fits.fz','_combs.fits.fz'],/allow
spawn,['fpack','-D','-Y',base+'_comb.fits'],/noshell
spawn,['fpack','-D','-Y',base+'_combs.fits'],/noshell
;; gzip _comb.mask.fits and _comb.bpm.fits
spawn,['gzip','-f',base+'_comb.mask.fits'],/noshell
spawn,['gzip','-f',base+'_comb.bpm.fits'],/noshell

;; Copy everything left back to the original directory
printlog,logfile,'Copying files back to '+bdir
allfiles = file_search(procdir+['*','.*.fits'],count=nallfiles)
FILE_MOVE,allfiles,bdir,/allow,/overwrite
FILE_DELETE,procdir  ; delete temporary processing directory

printlog,logfile,'DELVERED_FORCEBRICK done after '+strtrim(systime(1)-t0,2)+' sec.'

CD,curdir  ;; back to original directory

;; Create JOINT catalogs
printlog,logfile,''
printlog,logfile,'CREATE JOINT CATALOGS'
printlog,logfile,''

if keyword_set(redo) or keyword_set(update) then jntredo=1 else jntredo=0
DELVERED_JOINTBRICKCATS,brick,logfile=logfile,redo=jntredo

JOURNAL

end
