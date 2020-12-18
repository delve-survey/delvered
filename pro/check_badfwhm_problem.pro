pro check_badfwhm_problem,nights,redo=redo

;; Check if there are any images that have bad FWHM compared to the
;; other chips of the same exposure

if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)


nnights = n_elements(nights)
if nnights eq 0 then begin
  nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
  nights = file_basename(nights)
endif

undefine,sumstr

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]
  print,strtrim(n+1,2),' ',inight

  outfile = delvedir+'exposures/'+inight+'/'+inight+'_badfwhm.fits'
  if file_test(outfile) eq 1 and not keyword_set(redo) then begin
    print,outfile,' already EXISTS and /redo NOT set'
    goto,NIGHTBOMB
  endif

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

  undefine,expstr,chstr

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])

    files = file_search(delvedir+'exposures/'+inight+'/'+ifield+'/chip??/'+ifield+'-????????_??.fits',count=nfiles)
    base = file_basename(files,'.fits')
    dum = strsplitter(base,'-',/extract)
    dum2 = strsplitter(reform(dum[1,*]),'_',/extract)
    expnum = reform(dum2[0,*])
    index = create_index(expnum)
    nexpnum = n_elements(index.value)    

    ;; Loop over exposures
    For e=0,nexpnum-1 do begin
      iexpnum = index.value[e]
      ind = index.index[index.lo[e]:index.hi[e]-1]
      nind = n_elements(ind)
      expfiles = files[ind] 
      chstr1 = replicate({file:'',expnum:'',chip:0L,exists:0B,fwhm:99.99,bad:0B},nind)
      chstr1.expnum = iexpnum

      ;; Loop over chips
      for c=0,nind-1 do begin
        chstr1[c].file = expfiles[c]
        chstr1[c].chip = photred_getchipnum(expfiles[c],{separator:'_',namps:62L})
        optfile = repstr(expfiles[c],'.fits','.opt')
        if file_test(optfile) eq 1 then begin
          chstr1[c].exists = 1B
          READLINE,optfile,olines
          fwind = where(strmid(olines,0,2) eq 'FW',nfwind)
          fwhm = float(first_el(strsplit(olines[fwind[0]],'=',/extract),/last))
          chstr1[c].fwhm = fwhm
        endif
      endfor

      ;; Do the opt files exist
      gd = where(chstr1.exists eq 1,ngd)
      if ngd eq 0 then begin
        print,'No good opt files for ',iexpnum
        goto,EXPBOMB
      endif

      ;; Check for bad FWHM
      medfwhm = median(chstr1[gd].fwhm)
      sigfwhm = mad(chstr1[gd].fwhm)
      bdfwhm = where(abs(chstr1.fwhm-medfwhm) gt (4*sigfwhm>0.5),nbdfwhm)
      if nbdfwhm gt 0 then chstr1[bdfwhm].bad = 1B
      print,strtrim(nbdfwhm,2),' bad chips for ',iexpnum

      expstr1 = {expnum:iexpnum,nfiles:nind,nbad:nbdfwhm,medfwhm:medfwhm,sigfwhm:sigfwhm}
      PUSH,expstr,expstr1
      PUSH,chstr,chstr1

      EXPBOMB:
    Endfor  ; exposure loop

    BOMB:
  Endfor ; field loop


  if n_elements(expstr) eq 0 then begin
    print,'No exposures'
    goto,NIGHTBOMB
  endif

  ;; Add these to the lists
  bdfwhm = where(chstr.bad eq 1,nbdfwhm)
  if nbdfwhm gt 0 then begin
    print,'Adding ',strtrim(nbdfwhm,2),' files to DAOPHOT.inlist and removing from DAOPHOT.success'
    ;; Add to DAOPHOT.inlist
    badfiles = chstr[bdfwhm].file
    len = strlen(delvedir+'exposures/'+inight+'/')
    badfiles = strmid(badfiles,len)
    WRITELINE,delvedir+'exposures/'+inight+'/logs/DAOPHOT.inlist',badfiles,/append
    ;; Remove from DAOPHOT.success
    READLINE,delvedir+'exposures/'+inight+'/logs/DAOPHOT.success',slines,count=nslines
    if nslines gt 0 then begin
      MATCH,file_basename(slines),file_basename(badfiles),ind1,ind2,count=nmatch
      if nmatch gt 0 then begin
        if nmatch lt nslines then begin
          REMOVE,ind1,slines
          WRITELINE,delvedir+'exposures/'+inight+'/logs/DAOPHOT.success',slines
        endif else begin
          TOUCHZERO,delvedir+'exposures/'+inight+'/logs/DAOPHOT.success'
        endelse
      endif
    endif
    ;; Remove the bad opt files
    FILE_DELETE,repstr(chstr[bdfwhm].file,'.fits','.opt'),/allow
  endif

  ;; Write out the summary information
  print,'Writing summary information to ',outfile
  MWRFITS,expstr,outfile,/create
  MWRFITS,chstr,outfile,/silent

  NIGHTBOMB:
Endfor  ; night loop

;stop

end
