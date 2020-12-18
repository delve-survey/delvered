pro fix_refcat

;; Fix reference catalogs near the RA=0/360 rollover

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

nights = file_search(delvedir+'exposures/201?????',/test_directory,count=nnights)
nights = file_basename(nights)

;; Night loop
For n=0,nnights-1 do begin
;For n=251,710 do begin
  inight = nights[n]
  print,strtrim(n+1,2),' ',inight

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])
    files = file_search(fdir[f]+'/chip28/*_refcat.fits.gz',count=nfiles)  ; central chip
    if nfiles gt 0 then begin
      refcat1 = mrdfits(files[0],1,/silent)
      if mean(refcat1.ra) gt 3 and mean(refcat1.ra) lt 357 then goto,BOMB
    endif else print,'no refcat.fits.gz files in '+fdir[f]

    print,fdir[f],' needs to be corrected'
    refcat = MRDFITS(delvedir+'exposures/'+inight+'/refcat/'+ifield+'_refcat.fits.gz',1)

    ;; coordinates relative to the center of the field
    cendec = mean(minmax(refcat.dec))
    if range(refcat.ra) gt 100 then begin
      ra = refcat.ra
      bd = where(ra gt 180,nbd)
      if nbd gt 0 then ra[bd]-=360
      cenra = mean(minmax(ra))
      if cenra lt 0 then cenra+=360
    endif else cenra=mean(minmax(refcat.ra))
    ROTSPHCEN,refcat.ra,refcat.dec,cenra,cendec,reflon,reflat,/gnomic

    ;; Get exposure list
    expfile = delvedir+'exposures/'+inight+'/'+inight+'_exposures.fits'
    if file_test(expfile) eq 0 then stop,expfile,' NOT FOUND'
    expstr = MRDFITS(expfile,1)

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
        if file_test(outfile1) eq 0 then begin
          if ccdnum ne 2 then print,outfile1,' NOT FOUND'
          goto,CHIPBOMB
        endif

        ;; Save the reference catalog for this chip
        hd = PHOTRED_READFILE(outfile1,/header)
        ;; Temporarily fix NAXIS1/2 values
        if sxpar(hd,'ZNAXIS1') ne 0 then begin
          sxaddpar,hd,'NAXIS1',sxpar(hd,'ZNAXIS1')
          sxaddpar,hd,'NAXIS2',sxpar(hd,'ZNAXIS2')
        endif
        nx = sxpar(hd,'naxis1')
        ny = sxpar(hd,'naxis2')
        HEAD_XYAD,hd,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
        ROTSPHCEN,vra,vdec,cenra,cendec,vlon,vlat,/gnomic
        ;gdrefcat = where(refcat.ra ge min(vra)-0.02 and refcat.ra le max(vra)+0.02 and $
        ;                 refcat.dec ge min(vdec)-0.02 and refcat.dec le max(vdec)+0.02,ngdrefcat)
        gdrefcat = where(reflon ge min(vlon)-0.02 and reflon le max(vlon)+0.02 and $
                         reflat ge min(vlat)-0.02 and reflat le max(vlat)+0.02,ngdrefcat)
        refcat1 = refcat[gdrefcat]
        refcatfile = chipdir+ifield+'-'+iexpnum+'_'+schip+'_refcat.fits'
        print,'  Writing ',refcatfile
        MWRFITS,refcat1,refcatfile,/create
        SPAWN,['gzip','-f',refcatfile],/noshell
        CHIPBOMB:
      Endfor  ; chip loop

    Endfor  ; exposure loop


    BOMB:
  Endfor ; field loop
  NIGHTBOMB:
Endfor  ; night loop

;stop

end
