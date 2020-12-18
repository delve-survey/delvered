pro delvered_prep_night,expfile,delvedir=delvedir,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo

;; This pre-processing script gets DELVE and community MC data ready
;; from the NOAO mass store to be processed with PHOTRED.

;; Defaults
;if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
if FILE_TEST(delvedir,/directory) eq 0 then FILE_MKDIR,delvedir
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
if n_elements(scriptsdir) gt 0 then scriptsdir=trailingslash(scriptsdir) else scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
if n_elements(irafdir) gt 0 then irafdir=trailingslash(irafdir) else irafdir='/home/dnidever/iraf/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir='/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Model magnitudes equation file
modeleqnfile = delvereddir+'params/modelmag_equations.txt'

;; Print information
print,'--------------------------------------------------------------------'
print,' PREPARING DELVE MC EXPOSURES FOR PHOTRED EXPOSURE-LEVEL PROCESSING'
print,'--------------------------------------------------------------------'
print,'DELVEDIR = ',delvedir
print,'SCRIPTSDIR = ',scriptsdir
print,'IRAFDIR = ',irafdir
print,'EXPOSUREDIR = ',expdir
print,'WORKINGDIR = ',workdir
print,'HOST = ',host

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; The photed.setup file
setup = ['##### REQUIRED #####',$
         'scriptsdir  '+scriptsdir,$
         'irafdir     '+irafdir,$
         'telescope   Blanco',$
         'instrument  DECAM',$
         'observatory CTIO',$
         'nmulti      10',$
         'nmulti_wcs       20',$
         'nmulti_daophot   20',$
         'nmulti_allframe  10',$
         'filtref     g,i,r,z,u,Y',$
         'modeleqnfile '+modeleqnfile,$
         'trans       delve.trans',$
         '##### OPTIONAL #####',$
         'sepfielddir  1',$
         'sepchipdir   1',$
         'keepmef      0',$
         'catformat    FITS',$
         'workdir      '+workdir,$
         'clean        1',$
         'skipcheck    1',$
         'redo         0',$
         'wcsrefname   GAIADR2',$
         'searchdist   20',$
         '#wcsrmslim   1.0',$
         'hyperthread   1',$
         'daopsfva      1',$
         'daofitradfwhm 1.0',$
         'psfstars      1',$
         'psfcomsrc     0',$
         'psfcomglobal  0',$
         'psfcomgauss   0',$
         '#mchmaxshift  50.0',$
         'finditer      1',$     ; 2->1 on 7/22/20
         'alfdetprog  sextractor',$
         '#alfnocmbimscale 0',$
         'alftrimcomb   0',$
         '#ddo51radoffset  1',$
         'cmbforce      1',$
         'keepinstr     1',$
         'avgmag        1',$
         'avgonlymag    0',$
         'todered       u,g,r,i,z,g-i',$
         'sumquick      1',$
         '##### STAGES #####',$
         '#rename',$
         '#split',$
         ' wcs',$
         ' daophot',$
         ' match',$
         '#allframe',$
         ' apcor',$
         ' astrom',$
         ' zeropoint',$
         ' calib',$
         ' combine',$
         ' deredden',$
         ' save',$
         '#html']

  ;; not enough inputs
  if n_elements(expfile) eq 0 then begin
    print,'Syntax - delvered_prep_night,expfile,delvedir=delvedir,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo'
    return
  endif
  ;; File not found
  if file_test(expfile) eq 0 then begin
    print,expfile,' NOT FOUND'
    return
  endif

  expstr1 = MRDFITS(expfile,1,/silent)
  expstr1.fluxfile = strtrim(expstr1.fluxfile,2)
  expstr1.maskfile = strtrim(expstr1.maskfile,2)
  expstr1.wtfile = strtrim(expstr1.wtfile,2)
  nexp = n_elements(expstr1)
  inight = file_basename(file_dirname(expfile))

  ;; Check that all of the exposures are for the same night
  observatory,'ctio',obs
  nightarr = strarr(nexp)
  for i=0,nexp-1 do begin
    dateobs = expstr1[i].date_obs
    dateobs = strjoin(strsplit(dateobs,' ',/extract),'T')
    ;; Get night number
    jd = date2jd(dateobs)
    jd -= obs.tz/24.0  ;; correct for time zone difference
    caldat,double(long(jd)),Month, Day, Year, Hour, Minute, Second  ;; get date for beginning of the night
    ;; Convert to YYYYMMDD format
    night1 = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
    nightarr[i] = night1
  endfor
  ui = uniq(nightarr,sort(nightarr))
  if n_elements(ui) gt 1 then begin
    print,'More than 1 night. '+strjoin(nightarr[ui],' ')
    return
  endif

  print,'' & print,'----------------------------------'
  print,'Night ',strtrim(inight,2),' - ',strtrim(nexp,2),' exposures'
  print,'=================================='

  ;; Get the WCS.inlist and WCS.success files and see if we already
  ;; have these exposures
  nightdir = expdir+strtrim(inight,2)+'/'
  setupfile = nightdir+'photred.setup'
  if file_test(nightdir,/directory) eq 1 then begin
    if file_test(setupfile) eq 0 then WRITELINE,setupfile,setup
    undefine,files
    READLIST,nightdir+'logs/WCS.inlist',inlist,/unique,setupdir=nightdir,count=ninlist,/silent
    if ninlist gt 0 then PUSH,files,inlist
    READLIST,nightdir+'logs/WCS.success',successlist,/unique,setupdir=nightdir,count=nsuccess,/silent
    if nsuccess gt 0 then PUSH,files,successlist
    ;; Some previous files found
    if n_elements(files) gt 0 then begin
      ;; Make these relative paths
      for k=0,n_elements(files)-1 do begin
        pos = strpos(files[k],inight)
        if pos ge 0 then begin
          len = strlen(inight)
          files[k] = strmid(files[k],pos+len+1)
        endif
        ;; I needed this for SMASH, but not for new data
        ;if strmid(files[k],4,5,/reverse_offset) eq '.fits' and file_test(files[k]) eq 0 then files[k]+='.fz'
      endfor
      files = files[uniq(files,sort(files))]  ;; only want unique ones
      ;; Make sure they exist
      gdf = where(file_test(nightdir+files) eq 1,ngdf,comp=bdf,ncomp=nbdf)
      if ngdf gt 0 then files=files[gdf] else undefine,files
    endif
    nfiles = n_elements(files)
    ;; Remove the existing exposures
    if n_elements(files) gt 0 then begin
      fbase = PHOTRED_GETFITSEXT(files,/basename)
      base = reform((strsplitter(fbase,'-',/extract))[0,*])
      MATCH,base,expstr1.expnum,ind1,ind2,/sort,count=nmatch
      if nmatch gt 0 then begin
        exptoadd = expstr1
        if nmatch lt nexp then begin
          print,'Removing ',strtrim(nmatch,2),' exposures that already exist'
          REMOVE,ind2,exptoadd
        endif else undefine,exptoadd
      endif
    endif else exptoadd=expstr1
    ;; Load the existing "fields" file 
    if file_test(nightdir+'fields') eq 1 then begin
      oldfieldstr = IMPORTASCII(nightdir+'fields',fieldname=['shname','name'],fieldtypes=[7,7],/silent)
      add_tag,oldfieldstr,'fieldnum',0L,oldfieldstr
      oldfieldstr.fieldnum = long(strmid(oldfieldstr.shname,1))
    endif
  ;; No directory yet
  endif else begin
    print,strtrim(inight,2),' not started yet'
    FILE_MKDIR,[nightdir,nightdir+'logs']
    WRITELINE,setupfile,setup
    exptoadd = expstr1
  endelse
  ntoadd = n_elements(exptoadd)
  print,strtrim(ntoadd,2),' exposures to add'
  if ntoadd eq 0 then goto,NIGHTBOMB

  ;; Group the exposures
  fields = strarr(ntoadd)
  fieldcount = 1L
  if n_elements(oldfieldstr) gt 0 then fieldcount=max(oldfieldstr.fieldnum)+1
  undefine,newfieldstr
  if ntoadd gt 1 then begin
    matches = MATCHALL_SPH(exptoadd.ra,exptoadd.dec,exptoadd.ra,exptoadd.dec,0.1,nmatches)
    for j=0,ntoadd-1 do begin
      if matches[j+1] ne matches[j] then begin
        find = matches[matches[j]:matches[j+1]-1]
        ;; New field
        ;;   not already in a field
        if fields[j] eq '' then begin
          fieldname = 'F'+strtrim(fieldcount,2)
          fieldstr1 = {shname:fieldname,name:'',fieldnum:fieldcount}
          push,newfieldstr,fieldstr1
          gf = where(fields[find] eq '',ngf)  ; only get exposures not already in a group
          fields[find[gf]] = fieldname
          fieldcount++
        endif
      endif
    endfor
  endif else begin
    fields = 'F'+strtrim(fieldcount,2)
    fieldstr1 = {shname:fields,name:'',fieldnum:fieldcount}
    push,newfieldstr,fieldstr1
  endelse
  field_index = CREATE_INDEX(fields)
  nfields = n_elements(field_index.num)
  print,strtrim(nfields,2),' fields found'
  print,''


  ;; Loop over the fields
  undefine,outfiles
  outfiles = strarr(62*ntoadd)
  ocount = 0LL
  For f=0,nfields-1 do begin
    find = field_index.index[field_index.lo[f]:field_index.hi[f]]
    nfind = n_elements(find)
    ifield = field_index.value[f]
    print,ifield,' - ',strtrim(nfind,2),' exposures'
    fexptoadd = exptoadd[find]
    ;; Create the field directory
    fielddir = nightdir+ifield+'/'
    if FILE_TEST(fielddir,/directory) eq 0 then FILE_MKDIR,fielddir

    ;; Get Gaia DR2 and other reference data for this field
    if FILE_TEST(nightdir+'refcat/',/directory) eq 0 then FILE_MKDIR,nightdir+'refcat/'
    savefile = nightdir+'refcat/'+ifield+'_refcat.fits'
    undefine,refcat
    if file_test(savefile+'.gz') eq 1 then begin
      refcat = MRDFITS(savefile+'.gz',1,/silent)
      tra = mean(minmax(refcat.ra))
      tdec = mean(minmax(refcat.dec))
      if sphdist(fexptoadd[0].ra,fexptoadd[0].dec,tra,tdec,/deg) gt 0.05 then undefine,refcat
    endif
    if (file_test(savefile+'.gz') eq 0 and n_elements(refcat) eq 0) or keyword_set(redo) then begin
      refcat = DELVERED_GETREFDATA('c4d-'+['u','g','r','i','z','Y','VR'],fexptoadd[0].ra,fexptoadd[0].dec,1.5,savefile=savefile)
      SPAWN,['gzip','-f',savefile],/noshell
    endif

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

    ;; Loop over exposures
    For e=0,nfind-1 do begin
      filename = fexptoadd[e].fluxfile
      filebase = PHOTRED_GETFITSEXT(filename,/basename)
      ;; Get number of extensions
      ;;   use symlink to make fits_open think it's a normal FITS file
      tmpfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
      tmpfile += '.fits'
      FILE_LINK,filename,tmpfile
      FITS_OPEN,tmpfile,fcb & FITS_CLOSE,fcb
      ;; Get the field longname
      if e eq 0 then begin
        hd0 = HEADFITS(tmpfile,exten=0)
        object = sxpar(hd0,'OBJECT')
        if strtrim(object,2) eq '' then object='Exposure'+strtrim(fexptoadd[e].expnum,2)
        MATCH,newfieldstr.shname,ifield,indfield
        newfieldstr[indfield].name = strcompress(object,/remove_all)
      endif

      ;; Get the CCDNUM for the extensions
      MATCH,decam.name,fcb.extname,ind1,ind2,/sort,count=nmatch
      ;extnum = ind2
      ccdnum = decam[ind1].ccdnum
      extname = fcb.extname[ind2]
      ;; Loop over the extensions and create the resource files
      For c=0,fcb.nextend-1 do begin
        schip = string(ccdnum[c],format='(i02)')
        ;; Create chip directory if needed
        chipdir = fielddir+'chip'+schip+'/'
        if FILE_TEST(chipdir,/directory) eq 0 then FILE_MKDIR,chipdir
        outfile1 = chipdir+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'
        WRITELINE,outfile1,''
        routfile1 = chipdir+'.'+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'
        rlines = ['fluxfile = '+strtrim(fexptoadd[e].fluxfile,2)+'['+strtrim(extname[c],2)+']',$
                  'wtfile = '+strtrim(fexptoadd[e].wtfile,2)+'['+strtrim(extname[c],2)+']',$
                  'maskfile = '+strtrim(fexptoadd[e].maskfile,2)+'['+strtrim(extname[c],2)+']']
        ;rlines = ['fluxfile = '+strtrim(fexptoadd[e].fluxfile,2)+'['+strtrim(extnum[c],2)+']',$
        ;          'wtfile = '+strtrim(fexptoadd[e].wtfile,2)+'['+strtrim(extnum[c],2)+']',$
        ;          'maskfile = '+strtrim(fexptoadd[e].maskfile,2)+'['+strtrim(extnum[c],2)+']']
        WRITELINE,routfile1,rlines
        outfiles[ocount] = ifield+'/chip'+schip+'/'+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'  ; relative path
        ocount++

        ;; Save the reference catalog for this chip
        hd = HEADFITS(tmpfile,exten=extname[c],/silent)
        ;; Temporarily fix NAXIS1/2 values
        sxaddpar,hd,'NAXIS1',sxpar(hd,'ZNAXIS1')
        sxaddpar,hd,'NAXIS2',sxpar(hd,'ZNAXIS2')
        nx = sxpar(hd,'naxis1')
        ny = sxpar(hd,'naxis2')
        HEAD_XYAD,hd,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
        ROTSPHCEN,vra,vdec,cenra,cendec,vlon,vlat,/gnomic
        offset = 0.02
        if abs(cendec) gt 70 then offset=0.2
        if abs(cendec) gt 80 then offset=0.4
        ;gdrefcat = where(refcat.ra ge min(vra)-offset and refcat.ra le max(vra)+offset and $
        ;                 refcat.dec ge min(vdec)-offset and refcat.dec le max(vdec)+offset,ngdrefcat)
        gdrefcat = where(reflon ge min(vlon)-offset and reflon le max(vlon)+offset and $
                         reflat ge min(vlat)-offset and reflat le max(vlat)+offset,ngdrefcat)
        refcat1 = refcat[gdrefcat]
        refcatfile = chipdir+ifield+'-'+fexptoadd[e].expnum+'_'+string(ccdnum[c],format='(i02)')+'_refcat.fits'
        MWRFITS,refcat1,refcatfile,/create
        SPAWN,['gzip','-f',refcatfile],/noshell
      Endfor  ; chip loop
      FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink
    Endfor  ; exposure loop
  Endfor  ; field loop
  outfiles = outfiles[0:ocount-1]  ;; trim outfiles

  ;; Create the WCS.inlist file
  if FILE_TEST(nightdir+'logs/',/directory) eq 0 then FILE_MKDIR,nightdir+'logs'
  WRITELINE,nightdir+'logs/WCS.inlist',outfiles

  ;; Write the "fields" file
  if n_elements(oldfieldstr) gt 0 then newfieldstr=[oldfieldstr,newfieldstr]

  ;; Make sure the field names are unique
  findex = create_index(newfieldstr.name)
  bd = where(findex.num gt 1,nbd)
  for j=0,nbd-1 do begin
    find = findex.index[findex.lo[bd[j]]:findex.hi[bd[j]]]
    nfind = n_elements(find)
    newfieldstr[find[1:*]].name = newfieldstr[find[0]].name + '_'+strtrim(lindgen(nfind-1)+2,2)
  endfor

  WRITELINE,nightdir+'fields',newfieldstr.shname+'   '+newfieldstr.name

  ;; Copy scripts into the directory
  ;;  I DON'T THINK I NEED TO DO THIS
  NIGHTBOMB:


;stop

end
