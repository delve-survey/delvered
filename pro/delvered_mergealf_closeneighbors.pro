;+
;
; DELVERED_MERGEALF_CLOSENEIGHBORS
;
; This detects "spurious" close neighbors stars and merges their measurements.
; These close neighbors are an artifact of the iterative detection algorithm
; in ALLFRAME.PRO/ALLFPREP.PRO and can happen when the PSF is not a good
; match (from lots of images with different seeing/PSF being combined
; and producing systematic residuals which get detected as sources in the
; second detection iteration.
;
; INPUTS:
;  combfile  The absolute path to the combined image.
;  sex       SExtractor structure from ALLFPREP.pro with the NDETITER column.
;  obj       Structure of objects produced by DELVERED_FORCEBRICK.PRO.
;  meas      Structure of measurements produced by DELVERED_FORCEBRICK.PRO.
;
; OUTPUTS:
;  newobj   New object structure with close neighbors merged.
;  newmeas  New measurement structure with close neighbors merged.
;
; USAGE:
;  IDL>delvered_mergealf_closeneighbors,combfile,sex,obj,meas,newobj,newmeas,logfile=logfile
;
; By D. Nidever  July 2020
;-

pro delvered_mergealf_closeneighbors,combfile,sex,obj,meas,newobj,newmeas,logfile=logfile

;; Not enough inputs
nobj = n_elements(obj)
nmeas = n_elements(meas)
undefine,newobj,newmeas
if n_elements(combfile) eq 0 or n_elements(sex) eq 0 or nobj eq 0 or nmeas eq 0 then begin
  print,'Syntax - delvered_mergealf_closeneighbors,combfile,sex,obj,meas,newobj,newmeas,logfile=logfile'
  return
endif

; Defaults
if n_elements(logfile) eq 0 then logfile=-1

;; Make sure the names don't have trailing spaces
obj.objid = strtrim(obj.objid,2)
meas.objid = strtrim(meas.objid,2)

;; Load the combined image
dir = file_dirname(combfile)+'/'
combbase = file_basename(combfile,'.fits.fz')
FITS_READ,combfile,cim,chead

;; Load NMG file
nmgfile = dir+combbase+'.nmg'
LOADALS,nmgfile,nmg
;; match to OBJ, some objects might already have been removed
objnum = long(reform((strsplitter(obj.objid,'.',/extract))[1,*]))
MATCH,nmg.id,objnum,ind1,ind2,/sort
si = sort(ind2)  ; sort by obj
ind1 = ind1[si]
ind2 = ind2[si]
nmg = nmg[ind1]

;; Load OPT file
LOADOPT,dir+combbase+'.opt',optstr
fwhm = optstr.fw

;; Add columns to NMG for easy access
temp = nmg
add_tags,temp,['ndet','mag_auto','ndetiter','imval','flux','gmaxflux','mmaxflux','mrad'],$
              ['0L','99.99','0L','0.0','0.0','0.0','0.0','0.0'],nmg  & undefine,temp
otags = tag_names(obj)
detind = where(strmid(otags,0,4) eq 'NDET',ndetind)
ndet = lonarr(n_elements(obj))
for i=0,ndetind-1 do ndet += obj.(detind[i])
nmg.ndet = ndet
match,nmg.id,sex.number,ind1,ind2,/sort
nmg[ind1].mag_auto = sex[ind2].mag_auto
nmg[ind1].ndetiter = sex[ind2].ndetiter


; Moffat and Gaussian peak values
; beta = 2.5 or 3.0
; profile = (beta-1)/(pi*alpha^2)*(1+(r/alpha)^2)^(-beta)
; max flux = flux * (beta-1)/(pi*alpha^2)
beta = 3.0d0
alpha = double( fwhm/(2*sqrt(2^(1.0/beta)-1)) )
flux = 10^((25.0-nmg.mag)/2.5)  ; in counts
mmaxflux = flux * (beta-1)/(!dpi*alpha^2)   ;; moffat peak value
gmaxflux = flux/(1.138*fwhm^2)              ;; 2D Gaussian peak value
;; Corrections, for some reason the NMG flux values are off
;;  when compared to the actual pixel values.
imval = interpolate(cim,nmg.x-1,nmg.y-1)
nmg.imval = imval
gcorr = median(imval/gmaxflux)
gmaxflux *= gcorr
mcorr = median(imval/mmaxflux)
mmaxflux *= mcorr
nmg.flux = flux
nmg.gmaxflux = gmaxflux
nmg.mmaxflux = mmaxflux

;; 2D Moffat profile
r = findgen(5*ceil(fwhm))
r2 = scale_vector(findgen(500),0,max(r))
thresh = median(cim)+3*mad(cim)
rad = fltarr(nobj)+0.5   ; 0.5 is the minimum
for i=0,nobj-1 do begin
  profile = (beta-1)/(!dpi*alpha^2)*(1+(r/alpha)^2)^(-beta)*flux[i]
  ;; correct the profile
  profile *= mcorr
  indhi = max(where(profile gt thresh))
  if indhi gt -1 then begin
    lo = (indhi-1) > 0
    hi = (indhi+1) < (n_elements(r)-1)
    r2 = scale_vector(findgen(21),r[lo],r[hi])
    interp,r,profile,r2,profile2
    indhi2 = max(where(profile2 gt thresh))
    rad[i] = r2[indhi2]
  endif
endfor
nmg.mrad = rad


;; Find the close matches
dcr = fwhm
matches = MATCHALL_2D(nmg.x,nmg.y,nmg.x,nmg.y,dcr,nmatches)
dupind = where(nmatches gt 1,ndupind)

;; We can't just choose the iter=1 objects that have iter=2
;; close neighbors, because sometimes ALLSTAR eliminated the original
;; iter=1 detection.


;; Loop over the objects that have close matches
mstr = replicate({num:0L,objid:'',primind:-1LL,primdetiter:-1,ngroup:0L,objind:lonarr(10)-1,nmeas:lonarr(10)-1,$
                  mag:fltarr(10)+99.99,roff:fltarr(10)+99.99,fratio:fltarr(10)+999999.0,$
                  inside:lonarr(10)-1,ninside:0L,ngrouporig:0L,objindorig:lonarr(10)-1},ndupind)
mstr.num = lindgen(ndupind)+1
mstr.objid = obj[dupind].objid
objtorem = lonarr(nobj)
prim = lonarr(nobj)
nei = lonarr(nobj)
count = 0LL
for i=0,ndupind-1 do begin
  j = dupind[i]
  if matches[j+1] ne matches[j] then begin
    ind = matches[matches[j]:matches[j+1]-1]
    ind = ind[sort(nmg[ind].mag)]  ;; sort by magnitude
    nind = n_elements(ind)

    ;; Sometimes the magnitudes order is not the same
    ;; as the order of the image values at the source positions

    ;; Primary, first detection iteration source
    ;;   if one of the sources was detected in the 1st iteration, then
    ;;   pick that one
    primind = where(nmg[ind].ndetiter eq 1,nprimind)
    if nprimind gt 0 then begin
      mstr[i].primdetiter = 1
      ;; If more than one then pick the brightest one
      if nprimind gt 1 then begin
        si = sort(nmg[ind[primind]].mag)
        primind = primind[si[0]]
      endif
      ;; Put primary in the 1st element
      if primind ne 0 then begin
        ind = [ind[primind],ind]
        REMOVE,primind+1,ind
      endif

    ;; All secondary iteration sources
    ;;   pick one with highest flux value at position in combined image
    endif else begin
      mstr[i].primdetiter = 2
      si = reverse(sort(nmg[ind].imval))
      ind = ind[si]
    endelse


    ;; Only merge now if this is the right primary
    ;if nmg[j].mag eq min(nmg[ind].mag) then begin
    if j eq ind[0] then begin
      mstr[i].primind = j
      mstr[i].ngroup = nind
      mstr[i].objind[0:nind-1] = ind
      mstr[i].mag[0:nind-1] = nmg[ind].mag
      mstr[i].nmeas[0:nind-1] = nmg[ind].ndet
      ;mstr[i].fluxratio[0:nind-1] = 10^((nmg[ind].mag-nmg[ind[0]].mag)/2.5)
      ;mstr[i].snr[0:nind-1] = 1.087/nmg[ind].err

      mstr[i].ngrouporig = mstr[i].ngroup
      mstr[i].objindorig = mstr[i].objind

      ;; More in-depth checking
      ;print,strtrim(i+1,2),' ',strtrim(mstr[i].ngroup,2),' ',nmg[mstr[i].primind].x,' ',nmg[mstr[i].primind].y,' ',nmg[mstr[i].primind].mag
      primind = mstr[i].primind
      x1 = nmg[primind].x
      y1 = nmg[primind].y
      off = 15
      pl = 0
      if pl eq 1 then begin
        displayc,cim,xr=[-off,off]+x1-1,yr=[-off,off]+y1-1,tit=strtrim(i+1,2)+' '+strtrim(mstr[i].primind,2),/z
        oplot,nmg[ind].x-1,nmg[ind].y-1,ps=1,sym=3.0,thick=2
      endif
      ;print,nmg[ind].mag


      ;; Comparing the primary profile at the positions of the neighbors
      roff = sqrt((nmg[ind].x-nmg[primind].x)^2+(nmg[ind].y-nmg[primind].y)^2)
      mstr[i].roff = roff
      profile = (beta-1)/(!dpi*alpha^2)*(1+(roff/alpha)^2)^(-beta)*nmg[primind].flux
      profile *= mcorr
      fratio = profile / nmg[ind].gmaxflux
      mstr[i].fratio[0:nind-1] = fratio
      inside = where(fratio gt 1.0,ninside)
      ;print,'Roff = ',roff
      ;print,'fratio = ',fratio
      ;print,'Ninside = ',strtrim(ninside,2)

      ;; This works well for fainter stars, but not as well for brighter
      ;; ones where the detections in the 2nd detection pass can be quite
      ;; bright as well, almost as bright as the primary
      ;; ~13th mag or brighter
      ;; sometimes the fratio test works on the bright stars and sometimes
      ;; it doesn't

      ;; if there are two nearly-equally bright sources near each other
      ;; and the fratio test fails, then maybe check the image to see
      ;; if there are two peaks there.

      ;; the actual flux values at the source positions
      pixval = interpolate(cim,round(nmg[ind].x-1),round(nmg[ind].y-1))

      inside = lonarr(mstr[i].ngroup)
      for k=0,mstr[i].ngroup-1 do begin
        ;; Check flux ratio
        if fratio[k] ge 0.99 then inside[k]=1
        ;; Check cross-sectional profile
        if inside[k] eq 0 then begin
          rdiff = roff[k]
          np = ceil(rdiff)+1
          xx = scale_vector(findgen(np),nmg[primind].x-1,nmg[ind[k]].x-1)
          yy = scale_vector(findgen(np),nmg[primind].y-1,nmg[ind[k]].y-1)
          xprof = interpolate(cim,xx,yy)
          ;; check for a dip, if none, then this is not a separate object
          dip = where(xprof lt xprof[np-1]*0.8,ndip)
          if ndip eq 0 then inside[k]=2
        endif    
      endfor
      ;print,'inside = ',inside
      dum = where(inside ge 1,ninside)
      mstr[i].inside[0:nind-1] = inside
      mstr[i].ninside = ninside

      ;; Remove neighbors that are "outside"
      bd = where(inside le 0,nbd)
      if nbd gt 0 then begin
        REMOVE,bd,ind,roff,fratio
        nind = n_elements(ind)       
        mstr[i].ngroup = nind
        mstr[i].objind = -1
        mstr[i].objind[0:nind-1] = ind
        mstr[i].nmeas = -1
        mstr[i].nmeas[0:nind-1] = nmg[ind].ndet
        mstr[i].mag = 99.99
        mstr[i].mag[0:nind-1] = nmg[ind].mag
        mstr[i].roff = 99.99
        mstr[i].roff[0:nind-1] = roff
        mstr[i].fratio = 999999.00
        mstr[i].fratio[0:nind-1] = fratio
      endif

      ;; Some neighbors to remove
      if ninside gt 1 then begin
        objtorem[ind[1:*]] = 1
        prim[ind[0]] = 1
        nei[ind[1:*]] += 1
      endif
      ;stop
    endif  ; this is the primary
  endif  ; multiple matches
endfor  ; groups loop

;; Keep good ones
gd = where(mstr.ngroup gt 1 and mstr.ninside gt 1,ngd)
mstr = mstr[gd]
nmstr = ngd
printlog,logfile,strtrim(nmstr,2),' objects with close neighbors within ',stringize(dcr,ndec=2),' pixels'


;; Deal with neighbors claimed by multiple objects
multi = where(nei gt 1,nmulti)
printlog,logfile,'Resolving ',strtrim(nmulti,2),' neighbors that are claimed by multiple objects'
for i=0,nmulti-1 do begin
  indmulti = multi[i]
  bdind1 = where(mstr.objind eq indmulti,nbd)
  bdind2 = array_indices(mstr.objind,bdind1)
  sind = reform(bdind2[0,*])  ;; source column index
  mind = reform(bdind2[1,*])  ;; mstr row index

  ;; Pick the "winner"
  choice = intarr(nbd)  ; which one do we choose
  ;; 1) If one has PRIMDETITER=1 and other 2, pick detiter=1
  primdetiter = mstr[mind].primdetiter
  first = where(primdetiter eq 1,nfirst)
  if nfirst eq 1 then choice[first] = 1
  ;; 2) Use the flux ratio (uses mag and distance).
  if max(choice) eq 0 then begin
    fratio = fltarr(nbd)
    for k=0,nbd-1 do fratio[k]=mstr[mind[k]].fratio[sind[k]]
    si = sort(fratio)
    if fratio[si[0]] lt fratio[si[1]]*0.8 then choice[si[0]] = 1
  endif
  ;; Pick brighter one
  if max(choice) eq 0 then begin
    mag = mstr[mind].mag[0]
    si = sort(mag)
    choice[si[0]] = 1
  endif
  
  ;; Fix up the structure for the "losers"
  lost = where(choice eq 0,nlost)
  for k=0,nlost-1 do begin
    mind1 = mind[lost[k]]
    sind1 = sind[lost[k]]

    ngroup = mstr[mind1].ngroup
    objind = mstr[mind1].objind[0:ngroup-1]
    nmeas = mstr[mind1].nmeas[0:ngroup-1]
    mag = mstr[mind1].mag[0:ngroup-1]
    roff = mstr[mind1].roff[0:ngroup-1]
    fratio = mstr[mind1].fratio[0:ngroup-1]
    inside = mstr[mind1].inside[0:ngroup-1]
    ninside = mstr[mind1].ninside
    REMOVE,sind1,objind,nmeas,mag,roff,fratio,inside
    ngroup -= 1
    ninside -= 1
    nobjind = n_elements(objind)

    mstr[mind1].ngroup = ngroup
    mstr[mind1].objind = -1
    mstr[mind1].objind[0:nobjind-1] = objind
    mstr[mind1].nmeas = -1
    mstr[mind1].nmeas[0:nobjind-1] = nmg[objind].ndet
    mstr[mind1].mag = 99.99
    mstr[mind1].mag[0:nobjind-1] = nmg[objind].mag
    mstr[mind1].roff = 99.99
    mstr[mind1].roff[0:nobjind-1] = roff
    mstr[mind1].fratio = 999999.00
    mstr[mind1].fratio[0:nobjind-1] = fratio
    mstr[mind1].inside = -1
    mstr[mind1].inside[0:nobjind-1] = inside
    mstr[mind1].ninside = ninside
  endfor

  ;; Fix the NEI array
  nei[multi[i]] = 1
endfor
;; Remove any that don't have close neighbors anymore
gd = where(mstr.ngroup gt 1 and mstr.ninside gt 1,ngd)
mstr = mstr[gd]
nmstr = ngd


;; Merging flux in the measurements
;;----------------------------------
;;  exposure by exposure
printlog,logfile,'Merging flux from close neighbors exposure-by-exposure in the measurements table'
objindex = create_index(meas.objid)
MATCH,objindex.value,obj.objid,ind1,ind2,/sort,count=nmatch
si = sort(ind2)  ; sort by obj index
ind1 = ind1[si]
ind2 = ind2[si]
newmeas = meas
nmeas = n_elements(meas)
meastorem = lonarr(nmeas)
;; Loop over the groups
for i=0,nmstr-1 do begin
  ;print,strtrim(i+1,2),' ',mstr[i].objid,' ',strtrim(mstr[i].ngroup,2)
  ;; Make an array of all the measurments of all objects of this group
  ngrp = mstr[i].ngroup
  ngrpmeas = total(mstr[i].nmeas[0:ngrp-1],/int)
  mind = lonarr(ngrpmeas)
  cnt = 0LL
  for k=0,ngrp-1 do begin
    oind = mstr[i].objind[k]  ;; index into OBJ for this object
    mind1 = objindex.index[objindex.lo[ind1[oind]]:objindex.hi[ind1[oind]]]
    nmind1 = n_elements(mind1)
    mind[cnt:cnt+nmind1-1] = mind1   ;; add to meas index for this group
    cnt += nmind1
  endfor  
  ;; Index the measurements on exposure
  expindex = create_index(meas[mind].exposure)
  nexp = n_elements(expindex.value)
  ;; Loop over exposures
  for k=0,nexp-1 do begin
    eind = mind[expindex.index[expindex.lo[k]:expindex.hi[k]]]
    neind = n_elements(eind)
    ;print,'  exposure ',expindex.value[k],' ',meas[eind[0]].filter 
    ;print,'  ',strjoin(meas[eind].objid,', ')
    ;print,'  ',strjoin(meas[eind].mag,', ')
    if neind gt 1 then begin
      ;; Find the primary
      pind = where(meas[eind].objid eq mstr[i].objid,npind)
      if npind gt 0 then begin
        primind = eind[pind]
        meas1 = meas[primind]
      ;; If the primary is not detected in this exposure, then use the brightest source
      endif else begin
        ;print,'Primary not found'
        si = sort(meas[eind].mag)
        pind = si[0]
        primind = eind[pind]
        meas1 = meas[primind]
        meas1.objid = mstr[i].objid   ; make sure it has the right objid
      endelse
      ;; Now combine the flux,  leave all the other parameters
      ;;  m = -2.5*alog10( 10.^(-m1/2.5) + 10.^(-m2/2.5) )
      imag = meas[eind].imag
      meas1.imag = -2.5*alog10( total( 10^(-imag/2.5) ))
      mag = meas[eind].mag
      meas1.mag = -2.5*alog10( total( 10^(-mag/2.5) ))
      ;; Update in NEWMEAS
      newmeas[primind] = meas1
      ;; Measurments to remove
      meastorem[eind] = 1     ; remove all except for the primary one 
      meastorem[primind] = 0  
    ;; Only ONE measurements, make sure it has the right OBJID
    endif else begin
      ;if newmeas[eind].objid ne mstr[i].objid then print,'Single meas. Primary not found'
      newmeas[eind].objid = mstr[i].objid
    endelse
  endfor
endfor



;; Removing neighbor measurements
indnei = where(nei eq 1,nindnei)
newobjindex = create_index(newmeas.objid)
MATCH,newobjindex.value,obj[indnei].objid,ind1,ind2,/sort,count=nmatch
;; remove measurements belonging to the removed neighbors
meastorem2 = lonarr(nmeas)
for i=0,nmatch-1 do begin
  ind = newobjindex.index[newobjindex.lo[ind1[i]]:newobjindex.hi[ind1[i]]]
  meastorem2[ind] = 1
endfor
torem2 = where(meastorem2 eq 1,ntorem)
printlog,logfile,strtrim(ntorem,2)+' neighbor measurements to remove'
if ntorem gt 0 then REMOVE,torem2,newmeas
;; sanity check, meastorem and meastorem2 should be the same

;; Add MERGE column to NEWOBJ
newobj = obj
add_tag,newobj,'neimerged',0,newobj
newobj[mstr.primind].neimerged = 1

;; Remove close neighbor objects from NEWOBJ
indnei = where(nei eq 1,nindnei)
printlog,logfile,strtrim(nindnei,2)+' neighbor objects to remove'
REMOVE,indnei,newobj

;; Average the photometry
DELVERED_SIMPLEAVGMEAS,newmeas,avgobj
;; copy the values to the final object structure
MATCH,newobj.objid,avgobj.objid,ind1,ind2,/sort,count=nmatch
temp = newobj[ind1]
STRUCT_ASSIGN,avgobj[ind2],temp,/nozero
newobj[ind1] = temp
otags = tag_names(newobj)
atags = tag_names(avgobj)
ufilter = meas[uniq(meas.filter,sort(meas.filter))].filter
for i=0,n_elements(ufiles)-1 do begin
  scatind = where(otags eq strupcase(ufilter[i])+'SCATTER',nscatind)
  rmsind = where(atags eq strupcase(ufilter[i])+'RMS',nrmsind)
  newobj[ind1].(scatind[0]) = avgobj[ind2].(rmsind[0])
endfor
;; copy over other columns from the old object structure
MATCH,newobj.objid,obj.objid,ind1,ind2,/sort,count=nmatch
cols = ['PROB','EBV','MAG_AUTO','MAGERR_AUTO','ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM','BRICKUNIQ']
for i=0,n_elements(cols)-1 do begin
  colind = where(otags eq cols[i],ncolind)
  newobj[ind1].(colind[0]) = obj[ind2].(colind[0])
endfor

;stop

end
