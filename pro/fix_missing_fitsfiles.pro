pro fix_missing_fitsfiles,nights

;; Some fits files are missing.  Make them at least 0B files.

;delvedir = '/dl1/users/dnidever/delve/'
;nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
;nights = file_basename(nights)
;cmd = "fix_resourcefiles,'"+nights+"'"
;cmddir = '/data0/dnidever/delve/'+strarr(nnights)
;pbs_daemon,cmd,cmddir,jobs=jobs,/idle,/hyper,prefix='rfix',nmulti=20,wait=1

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

nnights = n_elements(nights)
if nnights eq 0 then begin
  nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
  nights = file_basename(nights)
endif

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]
  print,strtrim(n+1,2),' ',inight

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])

    ;; Get exposure list
    ;expfile = delvedir+'exposures/'+inight+'/'+inight+'_exposures.fits'
    ;if file_test(expfile) eq 0 then stop,expfile,' NOT FOUND'
    ;expstr = MRDFITS(expfile,1)

    expfiles = file_search(delvedir+'exposures/'+inight+'/'+ifield+'/chip01/.'+ifield+'-????????_01.fits',count=nexpfiles)
    expnum = file_basename(expfiles,'_01.fits')
    expnum = strmid(expnum,1)  ; trim leading .
    len = strlen(ifield)
    expnum = strmid(expnum,len+1)
    nexpnum = n_elements(expnum)

    ;; Loop over exposures
    For e=0,nexpnum-1 do begin
      iexpnum = expnum[e]

      ;; Loop over the chips
      For c=0,61 do begin
        ccdnum = c+1
        schip = string(ccdnum,format='(i02)')
        rfile = delvedir+'exposures/'+inight+'/'+ifield+'/chip01/.'+ifield+'-'+iexpnum+'_'+schip+'.fits'
        if file_test(rfile) eq 0 then goto,CHIPBOMB
        fitsfile = delvedir+'exposures/'+inight+'/'+ifield+'/chip01/'+ifield+'-'+iexpnum+'_'+schip+'.fits'
        if file_test(fitsfile) eq 0 then begin
          print,fitsfile,' missing'
          touchzero,fitsfile
        endif
        CHIPBOMB:
      Endfor  ; chip loop
    Endfor  ; exposure loop

    BOMB:
  Endfor ; field loop

  NIGHTBOMB:
Endfor  ; night loop

;stop

end
