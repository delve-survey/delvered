;+
;
; DELVERED_ZEROPOINT_EXPOSURE
;
; This program calculates NSC-like exposure-level zero-point
; for a single DELVERED_EXPOSURE.
;
; This needs to be run at the main directory level, like DELVERED_ZEROPOINT.PRO.
;
; By D. Nidever  Dec 2019
;-

pro delvered_zeropoint_exposure,expname,redo=redo,modeleqnfile=modeleqnfile

;; Not enough inputs
if n_elements(expname) eq 0 then begin
  print,'Syntax - delvered_zeropoint_zeropoint,expname,redo=redo,modeleqnfile=modeleqnfile'
  return
endif

COMMON photred,setup

CD,current=curdir

logfile = -1
t0 = systime(1)

; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

; TELESCOPE
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif

; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif

;; Model magnitude equation file
if n_elements(modeleqnfile) eq 0 then begin
  modeleqnfile = READPAR(setup,'MODELEQNFILE')
  if (modeleqnfile eq '0' or modeleqnfile eq '' or modeleqnfile eq '-1') then begin
    printlog,logfile,'NO MODELEQNFILE FOUND.  Please add to >>photred.setup<< file'
    return
  endif
endif

; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]

;; Output file, are we redoing
field = first_el(strsplit(file_basename(expname),'-',/extract))  ; F1
outfile = field+'/'+expname+'_zeropoint.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  printlog,logfile,outfile+' EXISTS and /redo NOT SET'
  return
endif



;; Get all of the files from DAOPHOT.success
lists = PHOTRED_GETINPUT('ZEROPOINT','DAOPHOT.success',redo=redo,ext=['fits','fits.fz'],/silent)
if lists.ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif
fitsfiles = lists.inputlines
nfitsfiles = n_elements(fitsfiles)
;; Get the unique exposures
allbase = PHOTRED_GETFITSEXT(fitsfiles,/basename)
;; Remove the ccdnum suffix
expbase = allbase
if thisimager.namps gt 1 then $
  for i=0,nfitsfiles-1 do expbase[i] = first_el(strsplit(expbase[i],thisimager.separator,/extract))
;; Keep the files for this exposure
gdbase = where(expbase eq expname,ngdbase)
if ngdbase eq 0 then begin
  printlog,logfile,'No files for exposure='+expname+' in DAOPHOT.success'
  return
endif
fitsfiles = fitsfiles[gdbase]
nfitsfiles = ngdbase
allbase = allbase[gdbase]
expbase = expbase[gdbase]

;; Look for other chips for this exposure that were previously
;; processed successfully
if lists.nsuccesslines gt 0 then begin
  sallbase = PHOTRED_GETFITSEXT(lists.successlines,/basename)
  sexpbase = sallbase
  if thisimager.namps gt 1 then $
    for i=0,n_elements(sallbase)-1 do sexpbase[i] = first_el(strsplit(sexpbase[i],thisimager.separator,/extract))
  ind = where(sexpbase eq expname,nmatch)
  if nmatch gt 0 then begin
    push,allbase,sallbase[ind]
    push,expbase,sexpbase[ind]
    push,fitsfiles,lists.successlines[ind]
    ;; Make sure they are unique
    ui = uniq(allbase,sort(allbase))
    allbase = allbase[ui]
    expbase = expbase[ui]
    fitsfiles = fitsfiles[ui]
    nfitsfiles = n_elements(fitsfiles)
  endif
endif

;; Load the apcor.lst file
apcor = IMPORTASCII('apcor.lst',fieldnames=['name','value'],/noprint)
add_tag,apcor,'file','',apcor
apcor.file = file_basename(apcor.name,'a.del')


;; Get exposure-level zero-point
expstr = {name:'',field:'',filter:'',exptime:0.0,ncat:0L,nref:0L,num:0L,zpterm:99.99,zptermerr:9.99,translines:strarr(3)}
expfiles = fitsfiles
hd = PHOTRED_READFILE(fitsfiles[0],/header)
exptime = PHOTRED_GETEXPTIME(fitsfiles[0])
filter = PHOTRED_GETFILTER(fitsfiles[0])
printlog,logfile,expname+' '+filter+' '+stringize(exptime,ndec=1)+' '+strtrim(nfitsfiles,2)
expstr.name = expname
expstr.field = field
expstr.filter = filter
expstr.exptime = exptime

;if filter eq 'u' then begin
;  print,'Skipping all u-band exposures'
;  goto,BOMB
;endif

;; Load all of the ALS files and add coordinates
undefine,cat
For j=0,nfitsfiles-1 do begin
  dir1 = FILE_DIRNAME(fitsfiles[j])
  base1 = allbase[j]
  alsfile1 = dir1+'/'+base1+'.als'
  ccdnum = PHOTRED_GETCHIPNUM(base1,thisimager)
  ifield = first_el(strsplit(base1,'-',/extract))
  als0 = PHOTRED_READFILE(alsfile1,count=nals)
  ;; Change ID to a string
  schema = {id:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,mag:0.0,err:0.0,sky:0.0,iter:0.0,chi:0.0,sharp:0.0,cmag:0.0}
  als = REPLICATE(schema,nals)
  STRUCT_ASSIGN,als0,als,/nozero
  ;; Get the coordinates
  if strpos(fitsfiles[j],'.fits.fz') ne -1 then exten=1 else exten=0
  head1 = PHOTRED_READFILE(fitsfiles[j],exten=exten,/header)
  ;; Converting to IDL X/Y convention, starting at (0,0)
  ;; DAOPHOT has X/Y start at (1,1)
  HEAD_XYAD,head1,als.x-1.0,als.y-1.0,ra,dec,/degree
  als.ra = ra
  als.dec = dec
  ;; Modify the IDs
  if thisimager.namps gt 1 then begin
    ; FIELD_EXT.IDNUMBER, i.e. 190L182a_5.17366
    ;---------------------------------------------
    id2 = ifield+'_'+strtrim(ccdnum,2)+'.'+strtrim(als.id,2)
    als.id = id2
  endif else begin
    ;; Updating the IDs
    ;; FIELD.IDNUMBER, i.e. 190L182a.17366
    ;;---------------------------------------------
    id2 = ifield+'.'+strtrim(als.id,2)
    als.id = id2
  endelse
  ;; Get calibrated magnitude, take exptime and aperture correction
  ;; into acount
  MATCH,apcor.file,base1,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    apcor1 = apcor[ind1].value
    if finite(apcor1) eq 0 then begin
      print,'NaN aperture correction.  Using 0.0 instead'
      apcor1 = 0.0
    endif
  endif else begin
    print,'No aperture correction for ',base1
    apcor1 = 0.0
  endelse
  ;; aperture correction is SUBTRACTIVE, makes it brighter
  als.cmag = als.mag + 2.5*alog10(exptime) - apcor1
  PUSH,cat,als
Endfor
ncat = n_elements(cat)
expstr.ncat = ncat
;; Get central RA and DEC
if range(cat.ra) gt 100 then begin
  ra = cat.ra
  bd = where(ra gt 180,nbd)
  if nbd gt 0 then ra[bd]-=360
  cenra = mean(minmax(ra))
  if cenra lt 0 then cenra+=360
endif else cenra=mean(minmax(cat.ra))
cendec = mean(minmax(cat.dec))

;; Load the reference data for this field
reffile = 'refcat/'+field+'_refcat.fits.gz'
if file_test(reffile) eq 0 then begin
  print,reffile,' REFERENCE FILE NOT FOUND'
  push,failurelist,expfiles
  goto,BOMB
endif
ref = MRDFITS(reffile,1,/silent)  
nref = n_elements(ref)
reftags = tag_names(ref)
expstr.nref = nref

;; Crossmatch
dcr = 1.0
SRCMATCH,ref.ra,ref.dec,cat.ra,cat.dec,dcr,ind1,ind2,/sph,count=nmatch
printlog,logfile,strtrim(nmatch,2),' matches to reference catalog'
if nmatch eq 0 then begin
  printlog,logfile,'No matches to reference catalog'
  push,failurelist,expfiles
  goto,BOMB
endif
;; Matched catalogs
cat1 = cat[ind2]
ref1 = ref[ind1]
  
;; Get the model magnitudes
instfilt = 'c4d-'+filter
mmags = DELVERED_GETMODELMAG(ref1,instfilt,cendec,modeleqnfile)
if n_elements(mmags) eq 1 and mmags[0] lt -1000 then begin
  printlog,logfile,'No good model mags'
  push,failurelist,expfiles
  goto,BOMB
endif
;; Get good stars
gdcat = where(cat1.mag lt 50 and cat1.err lt 0.05 and abs(cat1.sharp) lt 1 and cat.chi lt 3 and $
              mmags[*,0] lt 50 and mmags[*,1] lt 5,ngdcat)
if ngdcat lt 10 then $
  gdcat = where(cat1.mag lt 50 and cat1.err lt 0.08 and abs(cat1.sharp) lt 1 and cat.chi lt 3 and $
                mmags[*,0] lt 50 and mmags[*,1] lt 5,ngdcat)
if ngdcat eq 0 then begin
  printlog,logfile,'No stars that pass all of the quality/error cuts'
  push,failurelist,expfiles
  goto,BOMB
endif
ref2 = ref1[gdcat]
mmags2 = mmags[gdcat,*]
cat2 = cat1[gdcat]

;; Matched structure
;mag2 = cat2.mag_auto + 2.5*alog10(exptime) ; correct for the exposure time
;mstr = {col:float(mmags2[*,2]),mag:float(mag2),model:float(mmags2[*,0]),err:float(mmags2[*,1]),ccdnum:long(cat2.ccdnum)}
;; Measure the zero-point
;NSC_INSTCAL_CALIBRATE_FITZPTERM,mstr,expstr,chstr
;expstr.zptype = 1

;; Calculate the zeropoint
diff = mmags2[*,0] - cat2.cmag
err = sqrt(mmags2[*,1]^2 + cat2.err^2)
;; Make a sigma cut
med = median([diff])
sig = mad([diff])
gd = where(abs(diff-med) lt 3*sig,ngd)
x = fltarr(ngdcat)
undefine,zpterm,zptermerr
zpterm = dln_poly_fit(x[gd],diff[gd],0,measure_errors=err[gd],sigma=zptermerr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
zpterm = zpterm[0]
zptermerr = zptermerr[0]
printlog,logfile,'ZPTERM = '+stringize(zpterm,ndec=4)+' +/- '+stringize(zptermerr,ndec=4)

;; Add to the exposure structure
expstr.num = ngd
expstr.zpterm = zpterm
expstr.zptermerr = zptermerr
;; Lines for PHOTRED trans file
;;  ZPTERM is a SUBTRACTIVE constant offset
expstr.translines = [expname+'  '+filter+'  '+filter+'-'+filter+'  '+string(-zpterm,format='(f7.4)')+'    0.0000    0.0000   0.0000   0.0000',$
                        '                     '+string(zptermerr,format='(f7.4)')+'    0.0000    0.0000   0.0000   0.0000','']
;  F5-00517150_43  G  G-R  -0.4089    0.1713   -0.1193   0.0000   0.0000
;                           0.0040   -0.0000    0.0001   0.0000   0.0000


;; Write out the exposure structure
printlog,logfile,'Wrting exposure structure to '+outfile
MWRFITS,expstr,outfile,/create

printlog,logfile,'dt = '+strtrim(systime(1)-t0,2)+' sec.'

BOMB:

if keyword_set(stp) then stop

end
