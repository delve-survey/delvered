;+
;
; DELVERED_AVGMEAS
;
; Average measurements for objects. OBJID must have been set.
;
; INPUTS:
;  expstr  Catalog of exposures and their metadata.
;  meas    Catalog of all measurements with OBJID already filled in.
;
; OUTPUTS:
;  obj     Catalog of unique objects with average photometry,
;            morphology, coordinates, etc.
;
; USAGE:
;  IDL>delvered_avgmeas,meas,obj
;
; By D.Nidever   June 2020
;
;-

pro delvered_avgmeas,expstr,meas,obj

;; Not enough inputs
undefine,obj
nmeas = n_elements(meas)
if nmeas eq 0 then begin
  print,'Syntax - delvered_avgmeas,meas,obj'
  return
endif
print,strtrim(nmeas,2),' measurements'

;; Create index on OBJID
oindex = create_index(meas.objid)
nobj = n_elements(oindex.value)
print,strtrim(nobj,2),' objects'

;; Index from MEAS to OBJ
measobj_index = lon64arr(nmeas)
for i=0,nobj-1 do measobj_index[oindex.index[oindex.lo[i]:oindex.hi[i]]] = i


;; Initialize the final object table, with ALL BANDS
obj_schema = {objid:'',brick:'',depthflag:0,nalfdetiter:0,neimerged:0,ra:0.0d0,dec:0.0d0,rarms:0.0,decrms:0.0,ndet:0,ndetall:0,$
              mlon:0.0d0,mlat:0.0d0,$
              umag:99.99,uerr:9.99,urms:99.99,ndetu:0,umagall:99.99,uerrall:9.99,urmsall:99.99,ndetallu:0,$
              gmag:99.99,gerr:9.99,grms:99.99,ndetg:0,gmagall:99.99,gerrall:9.99,grmsall:99.99,ndetallg:0,$
              rmag:99.99,rerr:9.99,rrms:99.99,ndetr:0,rmagall:99.99,rerrall:9.99,rrmsall:99.99,ndetallr:0,$
              imag:99.99,ierr:9.99,irms:99.99,ndeti:0,imagall:99.99,ierrall:9.99,irmsall:99.99,ndetalli:0,$
              zmag:99.99,zerr:9.99,zrms:99.99,ndetz:0,zmagall:99.99,zerrall:9.99,zrmsall:99.99,ndetallz:0,$
              ymag:99.99,yerr:9.99,yrms:99.99,ndety:0,ymagall:99.99,yerrall:9.99,yrmsall:99.99,ndetally:0,$
              chi:99.99,sharp:99.99,prob:99.99,ebv:99.99,mag_auto:99.99,magerr_auto:9.99,$
              asemi:999999.0,bsemi:999999.0,theta:999999.0,ellipticity:999999.0,fwhm:999999.0,$
              rmsvar:999999.0,madvar:999999.0,iqrvar:999999.0,etavar:999999.0,jvar:999999.0,$
              kvar:999999.0,chivar:999999.0,romsvar:999999.0,variable10sig:-1,nsigvar:999999.0,$
              gaia_match:0,gaia_xdist:999999.0,gaia_sourceid:0LL,gaia_ra:999999.0d0,gaia_ra_error:999999.0,$
              gaia_dec:999999.0d0,gaia_dec_error:999999.0,gaia_parallax:999999.0,gaia_parallax_error:999999.0,gaia_pmra:999999.0,$
              gaia_pmra_error:999999.0,gaia_pmdec:999999.0,gaia_pmdec_error:999999.0,gaia_gmag:999999.0,$
              gaia_gmag_error:999999.0,gaia_bpmag:999999.0,gaia_bpmag_error:999999.0,gaia_rpmag:999999.0,$
              gaia_rpmag_error:999999.0,$
              brickuniq:0B}
; depthflag: 1-allstar, single processing; 2-forced photometry; 3-both
obj = replicate(obj_schema,nobj)
obj.objid = oindex.value
otags = tag_names(obj)


;; Exposure names
;;   F23-00912982_15, trim off ccdnum
dum = reform((strsplitter(meas.exposure,'_',/extract))[0,*])
expnum = reform((strsplitter(dum,'-',/extract))[1,*])
eindex = create_index(expnum)
;; match EXPSTR to EINDEX
;;  also some exposures have 0 measurements
MATCH,eindex.value,expstr.expnum,ind1,ind2,/sort
mexpstr = expstr[ind2]   ;; matched expstr, ind1 is in the right order
nexposure = n_elements(mexpstr)
print,strtrim(nexposure,2),' exposures with measurements'
ufilter = mexpstr[uniq(mexpstr.filter,sort(mexpstr.filter))].filter
nufilter = n_elements(ufilter)
print,strtrim(nufilter,2),' unique filters'

;; Loop over unique filters
For f=0,nufilter-1 do begin
  filtind = where(mexpstr.filter eq ufilter[f],nfiltind)
  
  ;; Indices for the magnitude and errors in OBJ
  magind = where(otags eq strupcase(ufilter[f])+'MAGALL')
  errind = where(otags eq strupcase(ufilter[f])+'ERRALL')
  rmsind = where(otags eq strupcase(ufilter[f])+'RMSALL')
  nobsind = where(otags eq 'NDETALL'+strupcase(ufilter[f]))

  ;; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin

    ;; Now copy in the values, MEAS only has "good" detections
    mind = eindex.index[eindex.lo[filtind]:eindex.hi[filtind]]
    obj[measobj_index[mind]].(magind) = meas[mind].mag
    obj[measobj_index[mind]].(errind) = meas[mind].err
    obj[measobj_index[mind]].(nobsind) = 1
    ;; Depthflag
    ;;  OR combine to depthflag, 1-single-exposure, 2-forced, 3-both
    ;;  within a single exposure we can have some forced and some
    ;;  not, need to do this object by object
    ;; forced=0 -> depthflag=1, forced=1 -> depthflag=2, depthflag=forced+1
    obj[measobj_index[mind]].depthflag OR= (meas[mind].forced+1)

  ;; Multiple exposures for this filter to average
  endif else begin

    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(nobj)
    totalfluxwt = dblarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      totalwt[measobj_index[mind]] += 1.0d0/meas[mind].err^2
      totalfluxwt[measobj_index[mind]] += 2.5118864d^meas[mind].mag * (1.0d0/meas[mind].err^2)
      ;; Depthflag
      obj[measobj_index[mind]].depthflag OR= (meas[mind].forced+1)
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif

    ;; Measure rms, RMS
    ;;  sqrt(mean(diff^2))
    totaldiff = dblarr(nobj)
    numobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      totaldiff[measobj_index[mind]] += (newmag[measobj_index[mind]] - meas[mind].mag)^2
      numobs[measobj_index[mind]]++
    endfor
    newrms = sqrt( totaldiff/(numobs>1) )
    if nbdmag gt 0 then newrms[bdmag]=99.99

    ;; Set rms=99.99 for numobs=1
    oneobs = where(numobs eq 1,noneobs)
    if noneobs gt 0 then newrms[oneobs]=99.99

    obj.(magind) = newmag
    obj.(errind) = newerr
    obj.(rmsind) = newrms
    obj.(nobsind) = numobs
  endelse  ; combine multiple exposures for this filter

Endfor  ;; unique filter loop


;; Compute "best" photometry
;;--------------------------------

;; Loop over unique filters
For f=0,nufilter-1 do begin
  filtind = where(mexpstr.filter eq ufilter[f],nfiltind)
  
  ;; Indices for the magnitude and errors in OBJ
  magind = where(otags eq strupcase(ufilter[f])+'MAG')
  errind = where(otags eq strupcase(ufilter[f])+'ERR')
  rmsind = where(otags eq strupcase(ufilter[f])+'RMS')
  nobsind = where(otags eq 'NDET'+strupcase(ufilter[f]))

  ;; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin

    ;; Now copy in the values, MEAS only has "good" detections
    mind = eindex.index[eindex.lo[filtind]:eindex.hi[filtind]]
    obj[measobj_index[mind]].(magind) = meas[mind].mag
    obj[measobj_index[mind]].(errind) = meas[mind].err
    obj[measobj_index[mind]].(nobsind) = 1

  ;; Multiple exposures for this filter to average
  endif else begin

    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(nobj)
    totalfluxwt = dblarr(nobj)
    numobs = lonarr(nobj)
    numgdobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      wt = 1.0d0/meas[mind].err^2
      ;; Severely downweight shallow non-forced photometry with err>0.2 mag
      if mexpstr[filtind[k]].exptime lt 90 then begin
        bd = where(meas[mind].forced eq 0 and meas[mind].err gt 0.2,nbd,comp=gd,ncomp=ngd)
        if nbd gt 0 then wt[bd] *= 1e-4
      endif else gd=lindgen(n_elements(mind))
      numgdobs[measobj_index[mind[gd]]]++
      numobs[measobj_index[mind]]++
      totalwt[measobj_index[mind]] += wt
      totalfluxwt[measobj_index[mind]] += 2.5118864d^meas[mind].mag * wt
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif
    obj.(nobsind) = numgdobs
    ;; Some objects only have poor shallow exposures, drop them
    ;;   they are not very reliable, and we have the ALL phot for them
    badobs = where(numgdobs eq 0 and newmag lt 50,nbadobs)
    if nbadobs gt 0 then begin
      newmag[badobs] = 99.99
      newerr[badobs] = 9.99
    endif

    ;if ufilter[f] eq 'g' then stop,'gmag 1'

    ;; Measure rms, RMS
    ;;  sqrt(mean(diff^2))
    totaldiff = dblarr(nobj)
    numobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      if mexpstr[filtind[k]].exptime lt 90 then begin
        bd = where(meas[mind].forced eq 0 and meas[mind].err gt 0.2,nbd,comp=gd,ncomp=ngd)
      endif else begin
        gd = lindgen(n_elements(mind))
        ngd = n_elements(gd)
      endelse
      if ngd gt 0 then begin
        totaldiff[measobj_index[mind[gd]]] += (newmag[measobj_index[mind[gd]]] - meas[mind[gd]].mag)^2
        numobs[measobj_index[mind[gd]]]++
      endif
    endfor
    newrms = sqrt( totaldiff/(numobs>1) )
    if nbdmag gt 0 then newrms[bdmag]=99.99

    ;; Set rms=99.99 for numobs=1
    oneobs = where(numobs eq 1,noneobs)
    if noneobs gt 0 then newrms[oneobs]=99.99

    obj.(magind) = newmag
    obj.(errind) = newerr
    obj.(rmsind) = newrms

    ;if ufilter[f] eq 'g' then stop,'gmag 2'

  endelse  ; combine multiple exposures for this filter

Endfor  ;; unique filter loop


;; NDET and NDETALL
obj.ndet = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz+obj.ndety
obj.ndetall = obj.ndetallu+obj.ndetallg+obj.ndetallr+obj.ndetalli+obj.ndetallz+obj.ndetally



;; Coordinates, morphology, ebv, etc
;;-----------------------------------

;; Get average chi, sharp, prob
totchi = fltarr(nobj) & numchi = lon64arr(nobj)
totsharp = fltarr(nobj) & totwtsharp = fltarr(nobj)
;totprob = fltarr(nobj) & numprob = lon64arr(nobj)

;; Exposure loop
For i=0,nexposure-1 do begin
  ;; All measurements for this exposure
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]
  nmind = n_elements(mind)
  ;; CHI
  chi1 = fltarr(nobj)+!values.f_nan
  chi1[measobj_index[mind]] = meas[mind].chi
  gdchi = where(finite(chi1) eq 1 and chi1 lt 1e5,ngdchi)
  if ngdchi gt 0 then begin
    totchi[gdchi] += chi1[gdchi]
    numchi[gdchi]++
  endif
  ;; SHARP
  ;; use "weights" with wt=1 for normal values
  ;; and wt=0.0001 for short, allstar and S/N<5
  ;; Using weights will still use the low S/N
  ;; sharp value if it's the ONLY detection, which
  ;;  is what I want.
  sharp1 = fltarr(nobj)+!values.f_nan
  sharp1[measobj_index[mind]] = meas[mind].sharp
  gdsharp = where(finite(sharp1) eq 1 and sharp1 lt 1e5,ngdsharp)
  if ngdsharp gt 0 then begin
    wtsharp1 = fltarr(nobj)
    wtsharp1[gdsharp] = 1.0
    ;; Set wt to 1e-4 for ALLSTAR short exposures with S/N<5
    snr = fltarr(nobj)
    snr[measobj_index[mind]] = 1.087/meas[mind].err
    forced = intarr(nobj)-1
    forced[measobj_index[mind]] = meas[mind].forced
    lowsnrind = where(finite(sharp1) eq 1 and abs(sharp1) lt 1e5 and snr lt 5 and forced ne 1,nlowsnrind)
    if nlowsnrind gt 0 then wtsharp1[lowsnrind]=1e-4
    ;; weighted sum
    totsharp[gdsharp] += wtsharp1[gdsharp]*sharp1[gdsharp]
    totwtsharp[gdsharp] += wtsharp1[gdsharp]
  endif
Endfor
;; Make average CHI
gdchi = where(numchi gt 0,ngdchi)
avgchi = fltarr(nobj)+99.99
if ngdchi gt 0 then avgchi[gdchi]=totchi[gdchi]/numchi[gdchi]
obj.chi = avgchi
;; Make average SHARP
gdsharp = where(totwtsharp gt 0,ngdsharp)
avgsharp = fltarr(nobj)+99.99
if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/totwtsharp[gdsharp]
obj.sharp = avgsharp


;; Weighted RA/DEC
totalwt = dblarr(nobj)
totalrawt = dblarr(nobj)
totaldecwt = dblarr(nobj)
for i=0,nexposure-1 do begin
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]  
  totalwt[measobj_index[mind]] += 1.0d0/meas[mind].err^2
  totalrawt[measobj_index[mind]] += meas[mind].ra * (1.0d0/meas[mind].err^2)
  totaldecwt[measobj_index[mind]] += meas[mind].dec * (1.0d0/meas[mind].err^2)
endfor
newra = totalrawt/totalwt
newdec = totaldecwt/totalwt
obj.ra = newra
obj.dec = newdec


; RA/DEC RMS
;  sqrt(mean(diff^2))
totalradiff = dblarr(nobj)
totaldecdiff = dblarr(nobj)
for i=0,nexposure-1 do begin
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]  
  totalradiff[measobj_index[mind]] += (obj[measobj_index[mind]].ra - meas[mind].ra)^2
  totaldecdiff[measobj_index[mind]] += (obj[measobj_index[mind]].dec - meas[mind].dec)^2
endfor
newrarms = sqrt( totalradiff/(obj.ndetall>1) ) * 3600 * cos(obj.dec/!radeg)
newdecrms = sqrt( totaldecdiff/(obj.ndetall>1) ) * 3600
bd = where(obj.ndetall eq 0,nbd)
if nbd gt 0 then newrarms[bd]=99.99
if nbd gt 0 then newdecrms[bd]=99.99
; Set rms=99.99 for numobs=1
oneobs = where(obj.ndetall eq 1,noneobs)
if noneobs gt 0 then newrarms[oneobs]=99.99
if noneobs gt 0 then newdecrms[oneobs]=99.99
obj.rarms = newrarms
obj.decrms = newdecrms

;; EBV
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
obj.ebv = dust_getval(glon,glat,/noloop,/interp)

end
