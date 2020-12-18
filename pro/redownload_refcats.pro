pro redownload_refcats,inight

;; Redownload reference catlaogs

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

if n_elements(inight) eq 0 then begin
  print,'Need to input night number'
  return
endif

  print,inight

  nightdir = delvedir+'exposures/'+inight+'/'
  fdir = file_search(nightdir+'F*',/test_directory,count=nfdir)

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])
    files = file_search(fdir[f]+'/chip28/*_refcat.fits.gz',count=nfiles)  ; central chip
    if nfiles eq 0 then stop,'No '+fdir[f]+'/chip28/*_refcat.fits.gz file found'
    refcat1 = mrdfits(files[0],1,/silent)

    ;; coordinates relative to the center of the field
    cendec = mean(minmax(refcat1.dec))
    if range(refcat1.ra) gt 100 then begin
      ra = refcat1.ra
      bd = where(ra gt 180,nbd)
      if nbd gt 0 then ra[bd]-=360
      cenra = mean(minmax(ra))
      if cenra lt 0 then cenra+=360
    endif else cenra=mean(minmax(refcat1.ra))

    ;; Redownload the referencec catalog
    savefile = nightdir+'refcat/'+ifield+'_refcat.fits'
    refcat = DELVERED_GETREFDATA('c4d-'+['u','g','r','i','z','Y','VR'],cenra,cendec,1.5,savefile=savefile)
    SPAWN,['gzip','-f',savefile],/noshell

    ;; coordinates relative to the center of the field for the FINAL CATALOG
    cendec = mean(minmax(refcat.dec))
    if range(refcat.ra) gt 100 then begin
      ra = refcat1.ra
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
        off = 0.2
        gdrefcat = where(reflon ge min(vlon)-off and reflon le max(vlon)+off and $
                         reflat ge min(vlat)-off and reflat le max(vlat)+off,ngdrefcat)
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

;stop

end
