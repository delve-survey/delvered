pro fix_resourcefiles_extname,nights

;; Fix resources files to use extension names instead of extension
;; numbers, because the extension numbers can be different in the flux
;; and mask files.

;delvedir = '/dl1/users/dnidever/delve/'
;nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
;nights = file_basename(nights)
;cmd = "fix_resourcefiles,'"+nights+"'"
;cmddir = '/data0/dnidever/delve/'+strarr(nnights)
;pbs_daemon,cmd,cmddir,jobs=jobs,/idle,/hyper,prefix='rfix',nmulti=20,wait=1

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)


nnights = n_elements(nights)
if nnights eq 0 then begin
  nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
  nights = file_basename(nights)
endif

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]
  print,strtrim(n+1,2),' ',inight

  ;; Get the exposures file
  expfile = delvedir+'exposures/'+inight+'/'+inight+'_exposures.fits'
  if file_test(expfile) eq 0 then begin
    print,expfile,' NOT FOUND'
    goto,NIGHTBOMB
  endif
  expstr = mrdfits(expfile,1)
  expstr.expnum = strtrim(expstr.expnum,2)
  expstr.fluxfile = strtrim(expstr.fluxfile,2)
  expstr.maskfile = strtrim(expstr.maskfile,2)
  expstr.wtfile = strtrim(expstr.wtfile,2)

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

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

      ;; Get the associated mass store file
      MATCH,expstr.expnum,iexpnum,ind1,ind2,/sort,count=nmatch
      if nmatch eq 0 then stop,'no match to exposure list'
      expstr1 = expstr[ind1]
      
      ;; Get number of extensions
      ;;   use symlink to make fits_open think it's a normal
      ;;   FITS file
      ;tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
      ;tmpfile += '.fits'
      ;FILE_LINK,expstr1.fluxfile,tmpfile
      ;FITS_OPEN,tmpfile,fcb & FITS_CLOSE,fcb
      ;FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink 

      ;; Loop over the chips
      For c=0,61 do begin
        ccdnum = c+1
        ;if ccdnum eq 61 then goto,CHIPBOMB
        schip = string(ccdnum,format='(i02)')
        ;chipdir = fielddir+'chip'+schip+'/'
        ind = where(decam.ccdnum eq ccdnum,nind)
        extname = decam[ind[0]].name

        file = delvedir+'exposures/'+inight+'/'+ifield+'/chip'+schip+'/'+ifield+'-'+string(iexpnum,format='(i08)')+'_'+schip+'.fits'
        rfile = delvedir+'exposures/'+inight+'/'+ifield+'/chip'+schip+'/.'+ifield+'-'+string(iexpnum,format='(i08)')+'_'+schip+'.fits'
        if file_test(rfile) eq 0 then begin
          if ccdnum ne 2 and ccdnum ne 61 then print,file,' NOT FOUND'
          goto,CHIPBOMB
        endif

        ;; Load resource file
        READLINE,rfile,rlines
        ;; Remove blank spaces
        arr = strsplitter(rlines,' ',/extract)
        sz = size(arr)
        tags = ['fluxfile','wtfile','maskfile']
        for k=0,n_elements(tags)-1 do begin
          tind = where(arr[0,*] eq tags[k],ntind)
          if ntind eq 0 then stop,'cannot find '+tags[k]
          arr1 = reform(arr[*,tind])
          if n_elements(arr1) gt 3 then stop,'more than 3 elements'
          hi = strpos(arr1[2],'[')
          rline1 = arr1[0]+' = '+strmid(arr1[2],0,hi)+'['+extname+']'
          rlines[tind] = rline1
        endfor
        ;print,'  Fixing ',rfile
        WRITELINE,rfile,rlines

        CHIPBOMB:
      Endfor  ; chip loop

    Endfor  ; exposure loop

    BOMB:
  Endfor ; field loop

  NIGHTBOMB:
Endfor  ; night loop

;stop

end
