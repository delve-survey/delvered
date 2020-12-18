pro fix_resourcefiles,nights

;; Fix resources files for extra spaces

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

    expfiles = file_search(delvedir+'exposures/'+inight+'/'+ifield+'/chip01/'+ifield+'-????????_01.fits',count=nexpfiles)
    expnum = file_basename(expfiles,'_01.fits')
    len = strlen(ifield)
    expnum = strmid(expnum,len+1)
    nexpnum = n_elements(expnum)

    ;; Loop over exposures
    For e=0,nexpnum-1 do begin
      iexpnum = expnum[e]

      ;; Loop over the chips
      For c=0,61 do begin
        ccdnum = c+1
        if ccdnum eq 61 then goto,CHIPBOMB
        schip = string(ccdnum,format='(i02)')
        chipdir = delvedir+'exposures/'+inight+'/'+ifield+'/chip'+schip+'/'
        outfile1 = chipdir+ifield+'-'+iexpnum+'_'+schip+'.fits'
        if file_test(outfile1) eq 0 and ccdnum ne 2 then begin
          print,outfile1,' NOT FOUND'
          goto,CHIPBOMB
        endif

        ;; Check resource file
        ;;   some have extra spaces in them
        rfile = file_dirname(outfile1)+'/.'+file_basename(outfile1)
        if file_test(rfile) eq 1 then begin
          READLINE,rfile,rlines
          ;; Remove blank spaces
          arr = strsplitter(rlines,' ',/extract)
          sz = size(arr)
          if sz[1] gt 3 then begin
            tags = ['fluxfile','wtfile','maskfile']
            for k=0,n_elements(tags)-1 do begin
              tind = where(arr[0,*] eq tags[k],ntind)
              if ntind eq 0 then stop,'cannot find '+tags[k]
              arr1 = reform(arr[*,tind])
              if n_elements(arr1) gt 4 then stop,'more than 4 elements'
              rline1 = strjoin(arr1[0:2],' ')+arr1[3]
              rlines[tind] = rline1
            endfor
            ;print,'  Fixing ',rfile
            WRITELINE,rfile,rlines
          endif
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
