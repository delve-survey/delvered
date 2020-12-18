;+
;
; DELVERED_JOINTBRICKCATS
;
; Process a single DELVE brick and perform ALLFRAME FORCED photometry
;
; INPUTS:
;  brick   The DELVE brick name, e.g. 1234m045
;
; By D. Nidever  August 2019
;-

pro delvered_jointbrickcats,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,logfile=logfile

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
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
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
logsdir = brickdir+'logs/'
if n_elements(logfile) eq 0 then logfile=-1

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
;; Create brick directory if it doesn't exist yet
bdir = subdir+brick+'/'
;logfile = bdir+brick+'.'+logtime+'.log'
logfile = -1

;; Check output file
if file_test(bdir+brick+'_joint_object.fits.gz') eq 1 and not keyword_set(redo) then begin
  print,bdir+brick+'_joint_object.fits.gz EXISTS and /redo NOT set'
  return
endif


;; DECam imager
thisimager = {telescope:'BLANCO',instrument:'DECam',namps:62,separator:'_'}

;; Print information
printlog,logfile,'BRICK = ',brick
printlog,logfile,'DELVEDIR = ',delvedir
printlog,logfile,'SCRIPTSDIR = ',scriptsdir
printlog,logfile,'IRAFDIR = ',irafdir
printlog,logfile,'EXPOSUREDIR = ',expdir
printlog,logfile,'WORKINGDIR = ',workdir
printlog,logfile,'HOST = ',host

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
;gdch = where(chstr.fwhm*chstr.pixscale le fwhmthresh and chstr.exptime ge 90. and zpterm ge -0.5 and $
;             finite(zpterm) eq 1 and finite(chstr.apcor) eq 1 and $
;             (filt eq 'u' or filt eq 'g' or filt eq 'r' or filt eq 'i' or filt eq 'z' or filt eq 'Y'),ngdch)
gdch = where(chstr.fwhm*chstr.pixscale le fwhmthresh and zpterm ge -0.5 and $
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
chstr.fieldname = strtrim(chstr.fieldname,2)

;; Get the "reference" name
mchfile = file_search(bdir+'*_comb.mch',count=nmchfile)
if nmchfile eq 0 then begin
  printlog,logfile,'No comb.mch file found for '+brick
  return
endif
mchfile = mchfile[0]
mchbase = file_basename(mchfile,'_comb.mch')
combbase = mchbase+'_comb'


;; Load the forced photometry object catalog
objfile = bdir+brick+'_object.fits.gz'
printlog,logfile,'Loading forced photometry object catalog '+objfile
fobj = MRDFITS(objfile,1,/silent)
fobj.objid = strtrim(fobj.objid,2)
nfobj = n_elements(fobj)

;; Load the forced photometry measurement catalog
expfile = bdir+brick+'_expforced.fits.gz'
printlog,logfile,'Loading forced photometry measurement catalog '+expfile
fmeas = MRDFITS(expfile,1,/silent)
fmeas.id = strtrim(fmeas.id,2)
fmeas.objid = strtrim(fmeas.objid,2)
fmeas.exposure = strtrim(fmeas.exposure,2)
fmeas.filter = strtrim(fmeas.filter,2)
;; Create exposure index
fmeasexpindex = create_index(fmeas.exposure)

;; Load the forced meta-file
metafile = bdir+brick+'_meta.fits'  
fmeta = MRDFITS(metafile,1,/silent)
fmeta.file = strtrim(fmeta.file,2)
fmeta.base = strtrim(fmeta.base,2)
fmeta.expnum = strtrim(fmeta.expnum,2)


;; Add ALLFRAME detection iteration number
;;----------------------------------------
;; restore the SExtractor file
sexfile = file_search(bdir+'*_comb_allf.sex',count=nsexfile)
;; reverse sort by date if more than 1
if nsexfile gt 0 then begin
  sexinfo = file_info(sexfile)
  si = reverse(sort(sexinfo.ctime))
  sexfile = sexfile[si]
endif
sex = MRDFITS(sexfile[0],1,/silent)
nsex = n_elements(sex)
if tag_exist(sex,'NDETITER') eq 0 then begin
  printlog,logfile,'Adding ALLFRAME detection iteration number'
  ; get it from the log file
  logfiles = file_search(bdir+brick+'.????????????.log',count=nlogfiles)
  if nlogfiles eq 0 then begin
    print,'No logfile found'
    return
  endif
  lognlines = file_lines(logfiles)
  gd = where(lognlines gt 1,ngd)
  logfiles = logfiles[gd]
  nlogfiles = ngd
  info = file_info(logfiles)
  si = reverse(sort(info.mtime))
  logfiles = logfiles[si]
  info = info[si]
  logfile = logfiles[0]
  readline,logfile,loglines
  lo = where(stregex(loglines,'STEP 3: Running allframe prep',/boolean) eq 1,nlo)
  hi = where(stregex(loglines,'STEP 4: Running ALLFRAME',/boolean) eq 1,nhi)
  if nlo gt 0 and nhi gt 0 then begin
    loglines1 = loglines[lo+4:hi-5]
    ;; --Iteration 1--
    ;; Running SExtractor
    ;; SExtractor found 21671 sources
    ;; Running ALLSTAR
    ;; ALLSTAR found 18141 sources
    ;; 18141 new stars found
    ;; --Iteration 2--
    ;; Running SExtractor
    ;; SExtractor found 18307 sources
    ;; Running ALLSTAR
    ;; ALLSTAR found 26554 sources
    ;; 8413 new stars found
    sexind = where(stregex(loglines1,'^SExtractor found',/boolean) eq 1,nsexind)
    dum = strsplitter(loglines1[sexind],' ',/extract)
    sexnstars = long(reform(dum[2,*]))
    add_tag,sex,'NDETITER',0L,sex
    sex[0:sexnstars[0]-1].ndetiter = 1
    if nsex gt sexnstars[0] then sex[sexnstars[0]:*].ndetiter = 2
    ;; The code has problems sometimes
    ;scount = 0LL
    ;for i=0,nsexind-1 do begin
    ;  sex[scount:scount+sexnstars[i]-1].ndetiter = i+1
    ;  scount += sexnstars[i]
    ;endfor
  endif else begin
    ;; Something's wrong with the logfile
    printlog,logfile,'Cannont find ALLFPREP information in logfile.  Trying to get detection iteration information from SExtractor positions'

    sexfile = file_search(bdir+'*_comb_allf.sex',count=nsexfile)
    sex = MRDFITS(sexfile[0],1,/silent)
    add_tag,sex,'NDETITER',0L,sex
    ;; Between iterations the Y value should change by a lot
    bd = where(slope(sex.y_image) lt -max(sex.y_image)*0.5,nbd)
    sex[0:bd[0]].ndetiter = 1
    sex[bd[0]+1:*].ndetiter = 2
  endelse
endif

;; Now match to the object catalog
;; the obj.objid star number part should match with the SE number
;; 0988m505.404   0988m505.899   0988m505.1423
add_tag,fobj,'nalfdetiter',0,fobj
objnum = long(reform((strsplitter(fobj.objid,'.',/extract))[1,*]))
MATCH,sex.number,objnum,ind1,ind2,/sort,count=nmatch
fobj[ind2].nalfdetiter = sex[ind1].ndetiter


;; Remove measurements from duplicate chips
;;  this happened on some of the SMASH nights
chid = fmeta.expnum+'-'+strtrim(fmeta.chip,2)
ui = uniq(chid,sort(chid))
ui = ui[sort(ui)]  ; try to keep in same order
if n_elements(ui) lt n_elements(fmeta) then begin
  bdchip = lindgen(n_elements(fmeta))
  REMOVE,ui,bdchip
  nbdchip = n_elements(bdchip)
  printlog,logfile,'--- Removing measurements from '+strtrim(nbdchip,2)+' duplicate chips ---'
  ;; Remove duplicate measurements from fmeas
  MATCH,fmeasexpindex.value,fmeta[bdchip].base,ind1,ind2,/sort,count=nmatch
  for i=0,nmatch-1 do begin
    ind = fmeasexpindex.index[fmeasexpindex.lo[ind1[i]]:fmeasexpindex.hi[ind1[i]]]
    push,badmeasind,ind
  endfor
  printlog,logfile,strtrim(n_elements(badmeasind),2)+' duplcate measurements to remove'
  REMOVE,badmeasind,fmeas
  ;; Remove chips with no measurements from fmeta
  REMOVE,bdchip,fmeta
  ;; Remove fobj objects with no measurements
  ui = uniq(fmeas.objid,sort(fmeas.objid))
  fobjid = fmeas[ui].objid
  MATCH,fobj.objid,fobjid,ind1,ind2,/sort,count=nmatch
  if nmatch lt n_elements(fobj) then begin
    left = lindgen(n_elements(fobj))
    REMOVE,ind1,left
    printlog,logfile,'Removing '+strtrim(n_elements(left),2)+' forced objects with no measurements'
    REMOVE,left,fobj
  endif
  nfobj = n_elements(fobj)
  ;; Remake fmeasexpindex
  fmeasexpindex = create_index(fmeas.exposure)
endif


;; Remove measurements from bad half of chip 31 from forced measurements
;;----------------------------------------------------------------------
mjd = dblarr(n_elements(fmeta))
for i=0,n_elements(fmeta)-1 do mjd[i]=date2jd(fmeta[i].utdate+'T'+fmeta[i].uttime,/mjd)
bdchip = where(fmeta.chip eq 31 and mjd gt 56660,nbdchip)
undefine,badmeasind
if nbdchip gt 0 then begin
  MATCH,fmeasexpindex.value,fmeta[bdchip].base,ind1,ind2,/sort,count=nmatch
  for i=0,nmatch-1 do begin
    ind = fmeasexpindex.index[fmeasexpindex.lo[ind1[i]]:fmeasexpindex.hi[ind1[i]]]
    bdind = where(fmeas[ind].x gt 1000,nbdind,comp=gdind,ncomp=ngdind)
    if nbdind gt 0 then begin   ; some bad ones found
      push,badmeasind,ind[bdind]
      if ngdind eq 0 then begin   ; all bad
        printlog,logfile,'NO useful measurements in '+fmeasexpindex.value[ind1[i]]
        fmeta[bdchip[ind2[i]]].alf_nsources = 0
      endif else begin
        printlog,logfile,fmeasexpindex.value[ind1[i]]+' removing '+strtrim(nbdind,2)+' bad measurements, '+strtrim(ngdind,2)+' left'
        fmeta[bdchip[ind2[i]]].alf_nsources = ngdind
      endelse
    endif  ; some bad ones to remove
  endfor
  printlog,logfile,strtrim(n_elements(badmeasind),2)+' bad chip 31 measurements to remove'
  ;; Remove bad chip 31 measurements from fmeas
  if n_elements(badmeasind) gt 0 then begin
    REMOVE,badmeasind,fmeas
  endif
  ;; Remove chips with no measurements from fmeta
  bdfmeta = where(fmeta.alf_nsources eq 0,nbdfmeta)
  if nbdfmeta gt 0 then begin
    printlog,logfile,'Removing '+strtrim(nbdfmeta,2)+' chips with no measurements'
    REMOVE,bdfmeta,fmeta
  endif
  ;; Remove fobj objects with no measurements
  ui = uniq(fmeas.objid,sort(fmeas.objid))
  fobjid = fmeas[ui].objid
  MATCH,fobj.objid,fobjid,ind1,ind2,/sort,count=nmatch
  if nmatch lt n_elements(fobj) then begin
    left = lindgen(n_elements(fobj))
    REMOVE,ind1,left
    printlog,logfile,'Removing '+strtrim(n_elements(left),2)+' forced objects with no measurements'
    REMOVE,left,fobj
  endif
  nfobj = n_elements(fobj)
  ;; Remake fmeasexpindex
  fmeasexpindex = create_index(fmeas.exposure)
endif



;; Check the ALLFRAME astrometric solutions and removing data from
;;   bad chips from the forced measurements catalog
;;----------------------------------------------------------------
printlog,logfile,'--- Removing chips with bad ALLFRAME astrometric solutions ---'
mchfile = file_search(bdir+'*_comb.mch',count=nmchfile)
if nmchfile eq 0 then begin
  printlog,logfile,'No comb.mch file found for '+brick
  return
endif
mchfile = mchfile[0]
mchbase = file_basename(mchfile,'_comb.mch')
astcheck = CHECK_ALLFRAME_COORDTRANS(mchfile,/silent)
bdchip = where(astcheck.decstd gt 1.0 or finite(astcheck.decstd) eq 0,nbdchip)
if nbdchip eq 0 then begin
  printlog,logfile,'NO ALLFRAME astrometric solutions problems'  
endif else begin
  printlog,logfile,strtrim(nbdchip,2)+' chips found with bad ALLFRAME astrometric solutions'
  printlog,logfile,strjoin(file_basename(astcheck[bdchip].file),', ')
  printlog,logfile,'Removing from forced META structure'
  badchips = file_basename(astcheck[bdchip].file,'.alf')
  ;; Remove bad chips from fmeta
  MATCH,fmeta.base,badchips,ind1,ind2,/sort,count=nmatch
  ;; there might not be a match if it's just a ccdnum=31 chip
  ;; that was removed above.
  if nmatch gt 0 then begin
    fmeta0 = fmeta
    REMOVE,ind1,fmeta
    ;; Remove bad chip measurements from fmeas
    fmeaschipindex = create_index(fmeas.exposure)
    MATCH,fmeaschipindex.value,badchips,ind1,ind2,/sort,count=nmatch
    badind = lon64arr(total(fmeaschipindex.num[ind1],/int))
    cnt = 0LL
    for i=0,nmatch-1 do begin
      ind = fmeaschipindex.index[fmeaschipindex.lo[ind1[i]]:fmeaschipindex.hi[ind1[i]]]
      nind = n_elements(ind)
      badind[cnt:cnt+nind-1] = ind
      cnt += nind
    endfor
    printlog,logfile,'Removing '+strtrim(n_elements(badind),2)+' forced measurements from bad chips'
    REMOVE,badind,fmeas
    ;; Remove fobj objects with no measurements
    ui = uniq(fmeas.objid,sort(fmeas.objid))
    fobjid = fmeas[ui].objid
    MATCH,fobj.objid,fobjid,ind1,ind2,/sort,count=nmatch
    if nmatch lt n_elements(fobj) then begin
      left = lindgen(n_elements(fobj))
      REMOVE,ind1,left
      printlog,logfile,'Removing '+strtrim(n_elements(left),2)+' forced objects with no measurements'
      REMOVE,left,fobj
    endif
    nfobj = n_elements(fobj)
    ;; Remake fmeasexpindex
    fmeasexpindex = create_index(fmeas.exposure)
  endif
endelse


;; Merge close neighbors
;;----------------------
printlog,logfile,'--- Merging close neighbors ---'
combfits = bdir+combbase+'.fits.fz'
if max(fobj.nalfdetiter) gt 1 then begin
  DELVERED_MERGEALF_CLOSENEIGHBORS,combfits,sex,fobj,fmeas,newobj,newmeas,logfile=logfile
  fobj = newobj & undefine,newobj
  fmeas = newmeas & undefine,newmeas
  nfobj = n_elements(fobj)
  ;; Remake the exposure index
  fmeasexpindex = create_index(fmeas.exposure)
endif


;; Initialize the final measurement table
meas_schema = {id:'',objid:'',brick:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,forced:0B,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,$
               imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0}
meas = replicate(meas_schema,n_elements(fmeas))
struct_assign,fmeas,meas
meas.brick = brick
meas.forced = 1B
mcount = long64(n_elements(meas))

;; Initialize the final object table, with ALL BANDS
obj_schema = {objid:'',brick:'',forced:0B,x:999999.0,y:999999.0,ra:0.0d0,dec:0.0d0,nalfdetiter:0L,neimerged:0,$
              umag:99.99,uerr:9.99,uscatter:99.99,ndetu:0L,$
              gmag:99.99,gerr:9.99,gscatter:99.99,ndetg:0L,$
              rmag:99.99,rerr:9.99,rscatter:99.99,ndetr:0L,$
              imag:99.99,ierr:9.99,iscatter:99.99,ndeti:0L,$
              zmag:99.99,zerr:9.99,zscatter:99.99,ndetz:0L,$
              ymag:99.99,yerr:9.99,yscatter:99.99,ndety:0L,$
              chi:99.99,sharp:99.99,prob:99.99,ebv:99.99,mag_auto:99.99,magerr_auto:9.99,$
              asemi:999999.0,bsemi:999999.0,theta:999999.0,ellipticity:999999.0,fwhm:999999.9,$
              depthflag:0,brickuniq:0B}
; depthflag: 1-allstar, single processing; 2-forced photometry; 3-both
obj = replicate(obj_schema,nfobj)
struct_assign,fobj,obj,/nozero
obj.brick = brick
obj.depthflag = 2
obj.nalfdetiter = fobj.nalfdetiter
if tag_exist(fobj,'NEIMERGED') then obj.neimerged = fobj.neimerged
obj = add_elements(obj,100000L)
ocount = nfobj


;; Initialize the final meta structure
meta = fmeta
add_tag,meta,'nmeas',0L,meta
;; add NMEAS for the chips already used in the forced photometry
MATCH,meta.base,fmeasexpindex.value,ind1,ind2,/sort,count=nmatch
if nmatch gt 0 then meta[ind1].nmeas = fmeasexpindex.num[ind2]


;; Step 1: Loop over the chips that we used for the forced photometry
;;-------------------------------------------------------------------
;;  check for any measurements that were missed, e.g. for bright stars
printlog,logfile,' '
printlog,logfile,'Step 1: Adding ALLSTAR measurements for exposures already used in the forced photometry'
printlog,logfile,' '

ui = uniq(meta.expnum,sort(meta.expnum))
uexpnum = meta[ui].expnum
nuexpnum = n_elements(uexpnum)
printlog,logfile,strtrim(nuexpnum,2)+' exposures'
printlog,logfile,' '

;; Current objid count
;;  fobj.objid sometimes goes higher than the number of objects
dum = strsplitter(fobj.objid,'.',/extract)
objectidcount = max(long(reform(dum[1,*])))+1

;; Loop over forced photometry exposures
For e=0,nuexpnum-1 do begin
  eind = where(meta.expnum eq uexpnum[e],neind)
  printlog,logfile,strtrim(e+1,2)+' '+uexpnum[e]+'  '+strtrim(neind,2)+' chip files'

  measexpnew = replicate(meas_schema,100000L)
  mexpcount = 0LL

  ;; Loop over the chips in this exposure
  for i=0,neind-1 do begin
    meta1 = meta[eind[i]]
    chcat = LOAD_CHIPCAT(meta1)  ; forced=0
    if size(chcat,/type) ne 8 then goto,BOMB1
    nchcat = n_elements(chcat)
    printlog,logfile,'  chip '+strtrim(i+1,2)+' '+repstr(meta1.file,'.fits','.phot')+' '+strtrim(nchcat,2)
    ;; Only keep measurements INSIDE the brick/tile area
    HEAD_ADXY,tilestr.head,chcat.ra,chcat.dec,bx,by,/deg
    gdcat = where(bx ge 0 and bx le (tilestr.nx-1) and by ge 0 and by le (tilestr.ny-1),ngdcat)
    if ngdcat gt 0 then begin
      printlog,logfile,'  '+strtrim(ngdcat,2)+' ALS measurements overlap the brick'
      chcat = chcat[gdcat]
      nchcat = n_elements(chcat)
      ;; Check if they are already in the forced measurement catalog
      ;;  remove any overlaps
      fmeasexpind = where(fmeasexpindex.value eq meta1.base,nfmeasexpind)
      if nfmeasexpind gt 0 then begin
        ;; Find the measurements already in the forced measurement catalog
        allfmeasexpind = fmeasexpindex.index[fmeasexpindex.lo[fmeasexpind]:fmeasexpindex.hi[fmeasexpind]]
        fmeas1 = fmeas[allfmeasexpind]
        printlog,logfile,'  '+strtrim(n_elements(fmeas1),2)+' forced measurements from this chip'
        SRCMATCH,chcat.ra,chcat.dec,fmeas1.ra,fmeas1.dec,1.0,ind1,ind2,/sph,count=nmatch
        printlog,logfile,'  '+strtrim(nmatch,2)+' matches'
        ;; Remove any matches
        if nmatch gt 0 then begin
          ;; None left
          if nmatch eq nchcat then begin
            undefine,chcat
            nchcat = 0
          endif else begin
            REMOVE,ind1,chcat
            nchcat = n_elements(chcat)
          endelse
        endif  ; found some matches
      endif  ; found in forced measurement catalog

      ;; Add remaining measurements to MEASEXPNEW
      if nchcat gt 0 then begin
        printlog,logfile,'  Adding '+strtrim(nchcat,2)+' remaining measurements to MEASEXPNEW'
        if nchcat+mexpcount gt n_elements(measexpnew) then measexpnew=add_elements(measexpnew,100000L>nchcat)
        ;; Add all measurements
        temp = measexpnew[mexpcount:mexpcount+nchcat-1]
        struct_assign,chcat,temp,/nozero
        temp.brick = brick
        measexpnew[mexpcount:mexpcount+nchcat-1] = temp
        mexpcount += nchcat
        ;; Update meta
        meta[eind[i]].nmeas += nchcat
      endif else begin
        printlog,logfile,'  No measurements left to add'
      endelse

    ;; No chip measurements overlap the bricks
    endif else begin
      printlog,logfile,'  NO measurements overlap the brick'
    endelse
    printlog,logfile,' '
    BOMB1:
  endfor
  ;; No new measurements to add
  if mexpcount eq 0 then begin
    printlog,logfile,'NO new measurements from this exposure to add'
    goto,ENDBOMB1
  endif
  ;; trim MEASEXPNEW
  if mexpcount lt n_elements(measexpnew) then measexpnew=measexpnew[0:mexpcount-1]

  printlog,logfile,'  '+strtrim(mexpcount,2)+' measurements to add for this exposure'

  ;; ADD EXPOSURE MEASUREMENTS TO OBJECT CATALOG!!!
  SRCMATCH,obj[0:ocount-1].ra,obj[0:ocount-1].dec,measexpnew.ra,measexpnew.dec,1.0,ind1,ind2,/sph,count=nmatch

  ;; Some matches, add the photometry
  left = lindgen(mexpcount)
  nleft = mexpcount
  if nmatch gt 0 then begin
    printlog,logfile,'  '+strtrim(nmatch,2),' measurements matched to previous objects'
    tempobj = obj[ind1]
    tempmeas = measexpnew[ind2]
    ADD_MEAS2OBJ,tempobj,tempmeas
    tempobj.depthflag OR= 1
    obj[ind1] = tempobj
    ;; update meas OBJID
    measexpnew[ind2].objid = obj[ind1].objid
    ;; Remove these from the list
    if nmatch eq mexpcount then begin
      undefine,left
      nleft = 0
    endif else begin
      REMOVE,ind2,left
      nleft = n_elements(left)
    endelse
  endif

  ;; Some leftover, add new objects
  if nleft gt 0 then begin
    printlog,logfile,'  '+strtrim(nleft,2),' measurements added as new objects'
    newobj = replicate(obj_schema,nleft)
    tempmeas = measexpnew[left]
    struct_assign,tempmeas,newobj,/nozero
    ADD_MEAS2OBJ,newobj,tempmeas
    newobj.depthflag OR= 1
    ;; Assign new OBJID
    newobjid = brick+'.'+strtrim(lindgen(nleft)+1+objectidcount,2)
    objectidcount += nleft
    newobj.objid = newobjid
    measexpnew[left].objid = newobjid
    ;; Add to OBJ
    if ocount+nleft gt n_elements(obj) then obj=add_elements(obj,100000L>nleft)
    obj[ocount:ocount+nleft-1] = newobj
    ocount += nleft
  endif

  ;; Add to MEAS
  ;;  add elemeents
  if mexpcount+mcount gt n_elements(meas) then meas=add_elements(meas,100000L>mexpcount)
  meas[mcount:mcount+mexpcount-1] = measexpnew
  mcount += mexpcount

  ENDBOMB1:

  printlog,logfile,' '
Endfor ;; exposure loop

;; trim MEAS
if mcount lt n_elements(meas) then meas=meas[0:mcount-1]
;; trim OBJ
if ocount lt n_elements(obj) then obj=obj[0:ocount-1]

;stop


;; Step 2: Add in completely new exposures
;;----------------------------------------

printlog,logfile,' '
printlog,logfile,'Step 2: Adding ALLSTAR measurements for NEW exposures'
printlog,logfile,' '

;; Get leftover chips
MATCH,chstr.base,fmeta.base,ind1,ind2,/sort,count=nmatch
if nmatch eq n_elements(chstr) then begin
  printlog,logfile,'NO new exposures to add'
  goto,FINAL
endif else begin
  newmeta = chstr
  add_tag,newmeta,'nmeas',0L,newmeta
  REMOVE,ind1,newmeta
endelse
;; new exposures
ui = uniq(newmeta.expnum,sort(newmeta.expnum))
uexpnum = newmeta[ui].expnum
nuexpnum = n_elements(uexpnum)

printlog,logfile,strtrim(nuexpnum,2)+' new exposures to add'
printlog,logfile,' '

;; Loop over the new exposures
For e=0,nuexpnum-1 do begin
  eind = where(newmeta.expnum eq uexpnum[e],neind)
  printlog,logfile,strtrim(e+1,2)+' '+uexpnum[e]+'  '+strtrim(neind,2)+' chip files'

  measexpnew = replicate(meas_schema,100000L)
  mexpcount = 0LL

  ;; Loop over the chips in this exposure
  for i=0,neind-1 do begin
    meta1 = newmeta[eind[i]]
    chcat = LOAD_CHIPCAT(meta1)  ; forced=0
    if size(chcat,/type) ne 8 then goto,BOMB2
    nchcat = n_elements(chcat)
    printlog,logfile,'  chip '+strtrim(i+1,2)+' '+repstr(meta1.file,'.fits','.phot')+' '+strtrim(nchcat,2)
    ;; Only keep meaurements INSIDE the brick/tile area
    HEAD_ADXY,tilestr.head,chcat.ra,chcat.dec,bx,by,/deg
    gdcat = where(bx ge 0 and bx le (tilestr.nx-1) and by ge 0 and by le (tilestr.ny-1),ngdcat)
    if ngdcat gt 0 then begin
      printlog,logfile,'  '+strtrim(ngdcat,2)+' ALS measurements overlap the brick'
      chcat = chcat[gdcat]
      nchcat = n_elements(chcat)
      ;; Add measurements to MEASEXPNEW
      if nchcat gt 0 then begin
        printlog,logfile,'  Adding '+strtrim(nchcat,2)+' remaining measurements to MEASEXPNEW'
        if nchcat+mexpcount gt n_elements(measexpnew) then measexpnew=add_elements(measexpnew,100000L>nchcat)
        ;; Add all measurements
        temp = measexpnew[mexpcount:mexpcount+nchcat-1]
        struct_assign,chcat,temp,/nozero
        temp.brick = brick
        measexpnew[mexpcount:mexpcount+nchcat-1] = temp
        mexpcount += nchcat
        ;; Update meta
        newmeta[eind[i]].nmeas += nchcat
      endif else begin
        printlog,logfile,'  No measurements left to add'
      endelse

    endif else begin
      printlog,logfile,'  NO measurements overlap the brick'
    endelse
    printlog,logfile,' '
    BOMB2:
  endfor
  ;; No new measurments to add
  if mexpcount eq 0 then begin
    printlog,logfile,'NO new measurements from this exposure to add'
    goto,ENDBOMB2
  endif
  ;; trim MEASEXPNEW
  if mexpcount lt n_elements(measexpnew) then measexpnew=measexpnew[0:mexpcount-1]

  printlog,logfile,'  '+strtrim(mexpcount,2)+' measurements to add for this exposure'

  ;; ADD EXPOSURE MEASUREMENTS TO OBJECT CATALOG!!!
  SRCMATCH,obj[0:ocount-1].ra,obj[0:ocount-1].dec,measexpnew.ra,measexpnew.dec,1.0,ind1,ind2,/sph,count=nmatch

  ;; Some matches, add the photometry
  left = lindgen(mexpcount)
  nleft = mexpcount
  if nmatch gt 0 then begin
    printlog,logfile,'  '+strtrim(nmatch,2),' measurements matched to previous objects'
    tempobj = obj[ind1]
    tempmeas = measexpnew[ind2]
    ADD_MEAS2OBJ,tempobj,tempmeas
    tempobj.depthflag OR= 1
    obj[ind1] = tempobj
    ;; update meas OBJID
    measexpnew[ind2].objid = obj[ind1].objid
    ;; Remove these from the list
    if nmatch eq mexpcount then begin
      undefine,left
      nleft = 0
    endif else begin
      REMOVE,ind2,left
      nleft = n_elements(left)
    endelse
  endif

  ;; Some leftover, add new objects
  if nleft gt 0 then begin
    printlog,logfile,'  '+strtrim(nleft,2),' measurements added as new objects'
    newobj = replicate(obj_schema,nleft)
    tempmeas = measexpnew[left]
    struct_assign,tempmeas,newobj,/nozero
    ADD_MEAS2OBJ,newobj,tempmeas
    newobj.depthflag OR= 1
    ;; Assign new OBJID
    newobjid = brick+'.'+strtrim(lindgen(nleft)+1+objectidcount,2)
    objectidcount += nleft
    newobj.objid = newobjid
    measexpnew[left].objid = newobjid
    ;; Add to OBJ
    if ocount+nleft gt n_elements(obj) then obj=add_elements(obj,100000L>nleft)
    obj[ocount:ocount+nleft-1] = newobj
    ocount += nleft
  endif

  ;; Add to MEAS
  ;;  add elemeents
  if mexpcount+mcount gt n_elements(meas) then meas=add_elements(meas,100000L>mexpcount)
  meas[mcount:mcount+mexpcount-1] = measexpnew
  mcount += mexpcount

  ENDBOMB2:

  printlog,logfile,' '
Endfor ;; exposure loop

;; trim MEAS
if mcount lt n_elements(meas) then meas=meas[0:mcount-1]
;; trim OBJ
if ocount lt n_elements(obj) then obj=obj[0:ocount-1]

;; Add NEWMETA to META
meta = [meta,newmeta]


;; Step 3: Redo average object photometry, keep best measurements only
;;--------------------------------------------------------------------
;; alex suggested having two mags, one is the weighted of all, the
;; other is the best measurement
FINAL:

;; Exposure information
ui = uniq(meta.expnum,sort(meta.expnum))
expstr = meta[ui]
si = sort(expstr.expnum)   ;; sort by exposure
expstr = expstr[si]
dum = reform((strsplitter(meas.exposure,'_',/extract))[0,*])
expnum = reform((strsplitter(dum,'-',/extract))[1,*])
eindex = create_index(expnum)
match,expstr.expnum,eindex.value,ind1,ind2,/sort,count=nmatch
expstr[ind1].nmeas = eindex.num[ind2]


;; Average all of the measurements
;;  THIS CREATES A NEW OBJECT TABLE WITH THE FINAL SCHEMA!!!
oldobj = obj
undefine,obj
printlog,logfile,'--- AVERAGING THE PHOTOMETRY ---'
DELVERED_AVGMEAS,expstr,meas,obj


;; Copy over SExtractor information for the forced objects
;;   the measurement table only has ALF/ALS information
MATCH,obj.objid,fobj.objid,ind1,ind2,/sort,count=nmatch
obj[ind1].prob = fobj[ind2].prob
obj[ind1].mag_auto = fobj[ind2].mag_auto
obj[ind1].magerr_auto = fobj[ind2].magerr_auto
obj[ind1].asemi = fobj[ind2].asemi
obj[ind1].bsemi = fobj[ind2].bsemi
obj[ind1].theta = fobj[ind2].theta
obj[ind1].ellipticity = fobj[ind2].ellipticity
obj[ind1].fwhm = fobj[ind2].fwhm

;; Copy over other information
obj[ind1].nalfdetiter = fobj[ind2].nalfdetiter
if tag_exist(fobj,'NEIMERGED') then obj[ind1].neimerged = fobj[ind2].neimerged
obj.brick = brick

;; Calculate photometric variability metrics
printlog,logfile,'--- Calculating photometric variability metrics ---'
DELVERED_PHOTVAR,meas,obj

;; Fill in mlon/mlat
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
obj.mlon = mlon
obj.mlat = mlat


;; Fill in BRICKUNIQ
;; Getting objects that are in the UNIQUE brick area
if brickstr1.dec eq -90 then begin
  ;; the brick right at the pole does not have any RA limits
  ginside = where(obj.dec lt brickstr1.dec2,ninside)
endif else begin
  ginside = where(obj.ra ge brickstr1.ra1 and obj.ra lt brickstr1.ra2 and $
                  obj.dec ge brickstr1.dec1 and obj.dec lt brickstr1.dec2,ninside)
endelse
if ninside gt 0 then obj[ginside].brickuniq=1B

;; Get Gaia DR2 data
;; bricks are 0.25 x 0.25 deg, so half the diagonal is 0.177 deg
printlog,logfile,'Crossmatching with Gaia DR2'
gaia = DELVERED_GETREFCAT(brickstr1.ra,brickstr1.dec,0.2,'gaiadr2')
srcmatch,obj.ra,obj.dec,gaia.ra,gaia.dec,1.0,ind1,ind2,/sph,count=nmatch
printlog,logfile,strtrim(nmatch,2)+' Gaia DR2 matches'
if nmatch gt 0 then begin
  obj[ind1].gaia_match = 1
  obj[ind1].gaia_xdist = sphdist(obj[ind1].ra,obj[ind1].dec,gaia[ind2].ra,gaia[ind2].dec,/deg)*3600
  obj[ind1].gaia_sourceid = gaia[ind2].source
  obj[ind1].gaia_ra = gaia[ind2].ra
  obj[ind1].gaia_ra_error = gaia[ind2].ra_error
  obj[ind1].gaia_dec = gaia[ind2].dec
  obj[ind1].gaia_dec_error = gaia[ind2].dec_error
  obj[ind1].gaia_parallax = gaia[ind2].parallax
  obj[ind1].gaia_parallax_error = gaia[ind2].parallax_error
  obj[ind1].gaia_pmra = gaia[ind2].pmra
  obj[ind1].gaia_pmra_error = gaia[ind2].pmra_error
  obj[ind1].gaia_pmdec = gaia[ind2].pmdec
  obj[ind1].gaia_pmdec_error = gaia[ind2].pmdec_error
  obj[ind1].gaia_gmag = gaia[ind2].gmag
  gmag_error = 2.5*alog10(1.0+gaia[ind2].e_fg/gaia[ind2].fg)
  obj[ind1].gaia_gmag_error = gmag_error
  obj[ind1].gaia_bpmag = gaia[ind2].bp
  bpmag_error = 2.5*alog10(1.0+gaia[ind2].e_fbp/gaia[ind2].fbp)
  obj[ind1].gaia_bpmag_error = bpmag_error
  obj[ind1].gaia_rpmag = gaia[ind2].rp
  rpmag_error = 2.5*alog10(1.0+gaia[ind2].e_frp/gaia[ind2].frp)
  obj[ind1].gaia_rpmag_error = rpmag_error
endif else print,'NO Gaia DR2 matches'

;stop

;; Save JOINT files
;;-------------------
printlog,logfile,' '
printlog,logfile,'Saving joint catalogs'

;; Saving object photometry catalog
objfile = bdir+brick+'_joint_object.fits'
MWRFITS,obj,objfile,/create
spawn,['gzip','-f',objfile],/noshell

;; Saving measurement catalog
expfile = bdir+brick+'_joint_meas.fits'
MWRFITS,meas,expfile,/create
spawn,['gzip','-f',expfile],/noshell

;; Save metadata
metafile = bdir+brick+'_joint_meta.fits'
printlog,logfile,'Writing meta-data to '+metafile
MWRFITS,meta,metafile,/create

printlog,logfile,'dt = '+stringize(systime(1)-t0,ndec=1)+' sec.'

;stop

end
