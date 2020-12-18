pro check_resourcefiles_extnameproblem,nights,redo=redo

;; Check if there is any problem with extension number issues
;; with the flux and mask files.

;delvedir = '/dl1/users/dnidever/delve/'
;nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
;nights = file_basename(nights)
;cmd = "check_resourcefiles_extnameproblem,'"+nights+"'"
;cmddir = '/data0/dnidever/delve/'+strarr(nnights)
;pbs_daemon,cmd,cmddir,jobs=jobs,/idle,/hyper,prefix='rfix',nmulti=20,wait=1

if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)


nnights = n_elements(nights)
if nnights eq 0 then begin
  nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
  nights = file_basename(nights)
endif

;goto,summary

undefine,sumstr

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]
  print,strtrim(n+1,2),' ',inight

  outfile = delvedir+'exposures/rfile_check/'+inight+'_summary.fits'
  if file_test(outfile) eq 1 and not keyword_set(redo) then begin
    print,outfile,' EXISTS and /red NOT set'
    goto,NIGHTBOMB
  endif

  ;; Get the exposures file
  ;expfile = delvedir+'exposures/'+inight+'/'+inight+'_exposures.fits'
  ;if file_test(expfile) eq 0 then begin
  ;  print,expfile,' NOT FOUND'
  ;  goto,NIGHTBOMB
  ;endif
  ;expstr = mrdfits(expfile,1)
  ;expstr.expnum = strtrim(expstr.expnum,2)
  ;expstr.fluxfile = strtrim(expstr.fluxfile,2)
  ;expstr.maskfile = strtrim(expstr.maskfile,2)
  ;expstr.wtfile = strtrim(expstr.wtfile,2)

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

  undefine,sumstr

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])

    expfiles = file_search(delvedir+'exposures/'+inight+'/'+ifield+'/chip01/'+ifield+'-????????_01.fits',count=nexpfiles)
    expnum = file_basename(expfiles,'_01.fits')
    len = strlen(ifield)
    expnum = strmid(expnum,len+1)
    nexpnum = n_elements(expnum)

    ;; Loop over exposures
    For e=0,nexpnum-1 do begin
      iexpnum = expnum[e]

      ;;; Get the associated mass store file
      ;MATCH,expstr.expnum,iexpnum,ind1,ind2,/sort,count=nmatch
      ;if nmatch eq 0 then stop,'no match to exposure list'
      ;expstr1 = expstr[ind1]
      
      newstr = {night:'',field:'',expnum:'',fluxfile:'',maskfile:'',okay:0}
      newstr.night = inight
      newstr.field = ifield
      newstr.expnum = iexpnum
      newstr.okay = -1

      ;; Get flux and mask file from resource file
      rfile1 = delvedir+'exposures/'+inight+'/'+ifield+'/chip01/.'+ifield+'-'+iexpnum+'_01.fits'
      if file_test(rfile1) eq 0 then begin
        print,rfile1,' NOT FOUND.  Skipping this exposure'
        goto,EXPBOMB
      endif
      READLINE,rfile1,rlines
      arr = strtrim(strsplitter(rlines,' ',/extract),2)
      find = where(arr[0,*] eq 'fluxfile',nfind)
      fluxfile = arr[2,find[0]]
      hi = strpos(fluxfile,'[')
      fluxfile = strmid(fluxfile,0,hi)
      mind = where(arr[0,*] eq 'maskfile',nfind)
      maskfile = arr[2,mind[0]]
      hi = strpos(maskfile,'[')
      maskfile = strmid(maskfile,0,hi)

      ;; Get fluxfile and maskfile header/extension information
      ;;   use symlink to make fits_open think it's a normal
      ;;   FITS file
      tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
      tmpfile += '.fits'
      FILE_LINK,fluxfile,tmpfile
      FITS_OPEN,tmpfile,fcb & FITS_CLOSE,fcb
      FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink 
      tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
      tmpfile += '.fits'
      FILE_LINK,maskfile,tmpfile
      FITS_OPEN,tmpfile,mfcb & FITS_CLOSE,mfcb
      FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink 

      newstr.fluxfile = fluxfile
      newstr.maskfile = maskfile
      newstr.okay = 1  ;; good until found bad

      ;; Check if there are differences
      bd = where(fcb.extname ne mfcb.extname,nbd)
      if nbd gt 0 then begin
        print,ifield,' ',iexpnum,'  PROBLEM'
        newstr.okay = 0
      endif else print,ifield,' ',iexpnum,'  OK'

      PUSH,sumstr,newstr

      EXPBOMB:
    Endfor  ; exposure loop

    BOMB:
  Endfor ; field loop

  ;; Write out the summary information
  print,'Writing summary information to ',outfile
  MWRFITS,sumstr,outfile,/create

  NIGHTBOMB:
Endfor  ; night loop

stop

summary:

for i=0,nnights-1 do begin
  inight = nights[i]
  outfile = delvedir+'exposures/rfile_check/'+inight+'_summary.fits'
  if file_test(outfile) eq 1 then begin
    sumstr = mrdfits(outfile,1,/silent)
    bd = where(sumstr.okay eq 0,nbd)
    if nbd eq 0 then print,strtrim(i+1,2),' ',inight else print,strtrim(i+1,2),' ',inight,' ',strtrim(nbd,2),' bad'
    ;if nbd eq 0 then print,strtrim(i+1,2),' ',inight else print,strtrim(i+1,2),' ',inight,' ',strtrim(nbd,2),' bad: ',strjoin(sumstr[bd].expnum,' ')
  endif else print,outfile,' NOT FOUND'
endfor

;stop

end
