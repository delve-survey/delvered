;+
;
; DELVERED_EXPOSURES
;
; This is a script to run PHOTRED on DELVE MC exposures to create
; single-image level catalogs and PSFs.
;
; INPUTS:
;  input         What nights to run PHOTRED on.  Either an array or
;                  a range such as 20160101-20160506.
;  /uselocal     Do all processing in a local drive.  Off by default.
;  /startfresh   Start completely fresh.  Remove old files and start
;                  from the very beginning.
;  =doapcor      Control whether APCOR is run.
;  =dozeropoint  Control whether ZEROPOINT is run.
;
; OUTPUTS:
;  PHOTRED will be run on each night and a final summary file
;  will be created for each night.
;
; USAGE:
;  IDL>delvered_exposures,'20160101'
;
; By D. Nidever  Feb 2019
;-

pro delvered_exposures,input,delvedir=delvedir,redo=redo,dozeropoint=dozeropoint,doapcor=doapcor,$
                       uselocal=uselocal,startfresh=startfresh,stp=stp,nmulti=nmulti

;; Defaults
;if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
workdir = '/data0/dnidever/delve/'  ;; temporary work directory

COMMON photred,setup

;; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - delvered_exposures,input,delvedir=delvedir,redo=redo,stp=stp'
  return
endif

; Start the logfile
;------------------
; format is delvered_exposures.host.DATETIME.log
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
logfile = logsdir+'delvered_exposures.'+hostname+'.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

;; Get all of the nights
dirs = FILE_SEARCH(expdir+'20??????',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)

;; Parse the input nights
for i=0,n_elements(input)-1 do begin
  input1 = input[i]
  ;; See if there is a -
  if strpos(input1,'-') ne -1 then begin
    nightrange = strsplit(input1,'-',/extract)
    ind_dirs = where(numdirs ge long(nightrange[0]) and numdirs le long(nightrange[1]),nind_dirs)
    if nind_dirs gt 0 then push,nights,dirs[ind_dirs] else print,'No directories found matching ',input1
  endif else begin
    MATCH,dirs,input1,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then push,nights,input1 else print,input1,' directory not found'
  endelse
endfor
nnights = n_elements(nights)
if nnights eq 0 then begin
  print,'No nights to process'
  return
endif

;; Which email program is available
mailprog = ''
spawn,'which sendmail',out,errout
if file_test(out[0]) eq 1 then mailprog='sendmail'
if mailprog eq '' then begin
  spawn,'which mail',out,errout
  if file_test(out[0]) eq 1 then mailprog='mail'
endif

; Print info
;-----------
print,''
print,'#######################################################'
print,'Starting DELVERED_EXPOSURES   ',systime(0)
print,'Running on ',host
print,'#######################################################'
print,''


; Make sure we have the right printlog.pro, not Markwardt's version
tempprogs = strsplit(!path,':',/extract)+'/printlog.pro'
test = file_test(tempprogs)
ind = where(test eq 1,nind)
bd = where(stregex(tempprogs[ind],'markwardt',/boolean) eq 1,nbd)
if nbd gt 0 then begin
  baddir = file_dirname(tempprogs[ind[bd]])
  print,"There is a version of Markwardt's PRINTLOG.PRO in "+baddir
  print,'Please rename this program (i.e. printlog.pro.orig)'
  return
endif

print,'Processing ',strtrim(nnights,2),' nights of data'


;#########################################
;#  STARTING THE PROCESSING
;#########################################
CD,current=origdir

;; Night loop
FOR i=0,nnights-1 do begin
  t0 = systime(1)
  inight = nights[i]
  CD,expdir+inight
  print,'' & print,'================================='
  print,' RUNNING PHOTRED on ',inight
  print,'=================================' & print,''

  if file_test('photred.setup') eq 0 then begin
    print,'NO photred.setup'
    goto,nightbomb
  endif
  ; LOAD THE SETUP FILE
  undefine,setup
  PHOTRED_LOADSETUP,setup,count=count
  if (count lt 1) then begin
    print,'Problem with photred.setup file'
    goto,nightbomb
  endif

  ;; Starting fresh
  ;;-----------------
  if keyword_set(startfresh) then begin
    print,'Starting fresh.  Deleting old files'
    ;; Delete old files
    ;;--------------------
    ;; zero-out fits files
    fitsfiles = file_search('F*/chip*/F*_??.fits',count=nfitsfiles)
    if nfitsfiles gt 0 then begin
      print,'Deleting image-level files'
      touchzero,fitsfiles
      base = repstr(fitsfiles,'.fits','')
      ext = ['_cat.dat','.fits.head','.opt','.als.opt','.coo','.ap','.grp','.lst','.cmn.lst','.psfini.ap','.lst1','.lst1.chi',$
             '.lst2','.lst2.chi','.plst','.psf','.psf.log','.als','.als.inp','.log','.nei','.nst',$
             'a.ap','a.als','a.als.inp','a.log','s.fits.fz','.mch','.tfr','.raw','.ast','.input','.phot','.ast']
      for j=0,n_elements(ext)-1 do file_delete,base+ext[j],/allow
      ;; remove HEADER from resource files
      print,'Updating resource files'
      for j=0,nfitsfiles-1 do begin
        dir1 = file_dirname(fitsfiles[j])
        base1 = file_basename(fitsfiles[j])
        rfile1 = dir1+'/.'+base1
        if file_test(rfile1) eq 1 then begin
          READLINE,rfile1,rlines
          bd = where(stregex(rlines,'HEADER',/boolean) eq 1,nbd)
          if nbd gt 0 then begin
            REMOVE,bd,rlines
            WRITELINE,rfile1,rlines
          endif
        endif
      endfor

      ;; chip directory-level files
      print,'Deleting chip directory-level files'
      fitsdirs = file_dirname(fitsfiles)
      udirs = fitsdirs[uniq(fitsdirs,sort(fitsdirs))]
      ;; scripts: daophot.sh, lstfilter.py, srcfilter.pro, goodpsf.pro, apcor.opt
      for j=0,n_elements(udirs)-1 do file_delete,udirs[j]+'/'+['daophot.sh','lstfilter.py','srcfilter.pro','goodpsf.pro','photo.opt','apcor.opt'],/allow
      ;; batch files, wfit, dopt, dao, initpsf, match,
      for j=0,n_elements(udirs)-1 do begin
        bfiles = file_search(udirs[j]+'/'+['wfit*.batch*','dopt*.batch*','dao*.sh*','initpsf*batch*','match*batch*'],count=nbfiles)
        if nbfiles gt 0 then file_delete,bfiles,/allow
      endfor
      ;; filters
      file_delete,udirs+'/filters',/allow
      ;; dangling rsrcXXXXXX/ directories
      for j=0,n_elements(udirs)-1 do begin
        rdirs = file_search(udirs[j]+'/rsrc*',/test_directory,count=nrdirs)
        if nrdirs gt 0 then begin
          rfiles = file_search(rdirs+'/*',count=nrfiles)
          if nrfiles gt 0 then file_delete,rfiles,/allow
          file_delete,rdirs,/allow
        endif
      endfor

      ;; field-level files: cmb, dered, final
      fdirs = file_dirname(udirs)  ; field directories
      fdirs = fdirs[uniq(fdirs,sort(fdirs))]
      for j=0,n_elements(fdirs)-1 do begin
        ffiles = file_search(fdirs[j]+'/'+['*.cmb','*.dered','*.final','*.dat','*.fits.gz','*_summary.fits','filters'],count=nffiles)
        if nffiles gt 0 then file_delete,ffiles,/allow
      endfor
    endif

    ;; daogrow files
    dfiles = file_search('daogrow*/*',count=ndfiles)
    if ndfiles gt 0 then begin
      print,'Deleting daogrow files'
      file_delete,dfiles,/allow
      ddirs = file_dirname(dfiles)
      ddirs = ddirs[uniq(ddirs,sort(ddirs))]
      file_delete,ddirs,/allow
    endif

    ;; summary files, delvered night summary file
    sfiles = file_search(['apcor*','delve.trans','*.dat','*.final','*.fits.gz','*_summary.fits','extinction','filters','idlbatch','runbatch'],count=nsfiles)
    if nsfiles gt 0 then file_delete,sfiles,/allow

    ;; zero-out lists and log files
    lfiles = file_search('logs/*',count=nlfiles)
    if nlfiles gt 0 then file_delete,lfiles,/allow
    ;; Put all fits files in WCS.inlist
    WRITELINE,'logs/WCS.inlist',fitsfiles
  endif


  ;; Copy everything to the work directory
  if keyword_set(uselocal) then begin
    if file_test(workdir+inight) eq 0 then FILE_MKDIR,workdir+inight
    print,''
    print,'Copying all files to temporary directory >>'+workdir+inight+'<<'
    SPAWN,['rsync','-av',expdir+inight+'/',workdir+inight+'/'],out1,errout1,/noshell
    if n_elements(errout1) gt 1 or errout1[0] ne '' then begin
      print,'There was a problem with the rsync '
      printline,errout1
      return
    endif
    print,strtrim(n_elements(out1)-4,2)+' files/directories copied'
    printline,out1
    ;; Go to the temporary directory
    CD,workdir+inight
  endif


  ;; Make sure the WCS.inlist files are relative
  READLINE,'logs/WCS.inlist',wcslist,count=nwcslist
  if nwcslist gt 0 then begin
    wcslist = repstr(wcslist,expdir+inight+'/','')
    WRITELINE,'logs/WCS.inlist',wcslist
  endif

  if READPAR(setup,'WCS') ne '0' then DELVERED_WCS,redo=redo,nmulti=nmulti
  if READPAR(setup,'DAOPHOT') ne '0' then DELVERED_DAOPHOT,redo=redo,nmulti=nmulti
  if READPAR(setup,'MATCH') ne '0' then DELVERED_MATCH,redo=redo,nmulti=nmulti

  if n_elements(doapcor) eq 1 then begin
    if keyword_set(doapcor) then PHOTRED_APCOR,redo=redo
  endif else if READPAR(setup,'APCOR') ne '0' then PHOTRED_APCOR,redo=redo

  if READPAR(setup,'ASTROM') ne '0' then PHOTRED_ASTROM,redo=redo

  if n_elements(dozeropoint) eq 1 then begin
    if keyword_set(dozeropoint) then DELVERED_ZEROPOINT,redo=redo,nmulti=nmulti
  endif else if READPAR(setup,'ZEROPOINT') ne '0' then DELVERED_ZEROPOINT,redo=redo,nmulti=nmulti

  if READPAR(setup,'CALIB') ne '0' then PHOTRED_CALIB,redo=redo
  if READPAR(setup,'COMBINE') ne '0' then PHOTRED_COMBINE,redo=redo
  if READPAR(setup,'DEREDDEN') ne '0' then PHOTRED_DEREDDEN,redo=redo
  if READPAR(setup,'SAVE') ne '0' then PHOTRED_SAVE,redo=redo,/sumquick,nmulti=nmulti

  print,'DELVERED FINISHED'
  stages = ['WCS','DAOPHOT','MATCH','APCOR','ASTROM','ZEROPOINT','CALIB','COMBINE','DEREDDEN','SAVE']
  PHOTRED_SUMMARY,outlines=outlines,stages=stages,/quick

  ;; Copy everything back to permanent directory
  if keyword_set(uselocal) then begin
    print,''
    print,'Copying all files back to permanent directory >>'+expdir+inight+'<<'
    print,''
    SPAWN,['rsync','-av',workdir+inight+'/',expdir+inight+'/'],out2,errout2,/noshell
    if n_elements(errout2) gt 1 or errout2[0] ne '' then begin
      print,'There was a problem with the rsync.'
      printline,errout2
      return
    endif
    printline,out2

    ;; SHOULD I CHECK HERE THAT EVERYTHING WAS COPIED CORRECTLY???

    CD,origdir

    ;; Delete temporary directory
    print,''
    print,'Deleting temporary directory >>'+workdir+inight+'<<'
    print,''
    SPAWN,['rm','-R',workdir+inight],/noshell

    ;; Fix the absolute paths in the summary files
    print,'Fixing absolute paths in summary files'
    sumfiles = file_search(expdir+inight+'/*_summary.fits',count=nsumfiles)
    for j=0,nsumfiles-1 do begin
      print,strtrim(j+1,2),' ',sumfiles[j]
      expstr = mrdfits(sumfiles[j],1,/silent)
      chstr = mrdfits(sumfiles[j],2,/silent)
      chstr.file = repstr(chstr.file,workdir,expdir)
      MWRFITS,expstr,sumfiles[j],/create
      MWRFITS,chstr,sumfiles[j],/silent
    endfor
  endif  ; /uselocal

  ;; Create the nightly summary file
  DELVERED_NIGHTSUMMARY,inight,delvedir=delvedir,redo=redo

  dt = systime(1)-t0
  print,'dt = ',strtrim(dt,2),' sec.'

  ;; Send email that this night is done
  CASE mailprog of
  'sendmail': begin
    undefine,elines
    push,elines,'From: dnidever@noao.edu'
    push,elines,'To: dnidever@noao.edu'
    push,elines,'Subject: delvered_exposures night='+inight+' FINISHED on '+hostname
    push,elines,'Content-Type: text/html'
    push,elines,'MIME-Version: 1.0'
    push,elines,'<pre>'
    push,elines,'delvered_exposures night='+inight+' FINISHED'
    push,elines,systime(0)
    push,elines,'HOST='+hostname
    push,elines,'dt='+strtrim(dt,2)+' sec.'
    push,elines,''
    push,elines,outlines
    push,elines,'</pre>'
    tempfile = mktemp('mail',/nodot)
    writeline,tempfile,elines
    cmd = 'cat '+tempfile+' | sendmail -t'
    spawn,cmd,out,errout
    file_delete,tempfile
  end
  'mail': begin
    undefine,elines
    push,elines,'delvered_exposures night='+inight+' FINISHED'
    push,elines,systime(0)
    push,elines,'HOST='+hostname
    push,elines,'dt='+strtrim(dt,2)+' sec.'
    push,elines,''
    push,elines,outlines
    tempfile = mktemp('mail',/nodot)
    writeline,tempfile,elines
    cmd = 'cat '+tempfile+' | mail  -s "delvered_exposures night='+inight+' FINISHED on '+hostname+'" dnidever@noao.edu'
    spawn,cmd,out,errout
    file_delete,tempfile
  end
  else: print,'No mail program available'
  ENDCASE

  NIGHTBOMB:
ENDFOR

; End logfile
;------------
JOURNAL

;; Send out final email that we are done
sinput = strtrim(input,2)
if n_elements(sinput) gt 1 then sinput='['+strjoin(sinput,',')+']'
body = 'delvered_exposures '+sinput+' FINISHED at '+systime(0)+' on HOST='+hostname
cmd = 'echo "'+body+'" |  mail -s "delvered_exposures FINISHED on '+hostname+'" dnidever@noao.edu'
spawn,cmd,out,errout

if keyword_set(stp) then stop

end
