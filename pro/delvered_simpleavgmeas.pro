;+
;
; DELVERED_SIMPLEAVGMEAS
;
; Average measurements for objects. OBJID must have been set.
; Like DELVERED_AVGMEAS.PRO, nothing fancy.
;
; INPUTS:
;  meas    Catalog of all measurements with OBJID already filled in.
;
; OUTPUTS:
;  obj     Catalog of unique objects with average photometry,
;            morphology, coordinates, etc.
;
; USAGE:
;  IDL>delvered_simpleavgmeas,meas,obj
;
; By D.Nidever   July 2020
;
;-

pro delvered_simpleavgmeas,meas,obj

;; Not enough inputs
undefine,obj
nmeas = n_elements(meas)
if nmeas eq 0 then begin
  print,'Syntax - delvered_simpleavgmeas,meas,obj'
  return
endif

;; Create index on OBJID
oindex = create_index(meas.objid)
nobj = n_elements(oindex.value)

;; Index from MEAS to OBJ
measobj_index = lon64arr(nmeas)
for i=0,nobj-1 do measobj_index[oindex.index[oindex.lo[i]:oindex.hi[i]]] = i


;; Initialize the final object table, with ALL BANDS
obj_schema = {objid:'',depthflag:0,ra:0.0d0,dec:0.0d0,rarms:0.0,decrms:0.0,ndet:0,$
              umag:99.99,uerr:9.99,urms:99.99,ndetu:0,$
              gmag:99.99,gerr:9.99,grms:99.99,ndetg:0,$
              rmag:99.99,rerr:9.99,rrms:99.99,ndetr:0,$
              imag:99.99,ierr:9.99,irms:99.99,ndeti:0,$
              zmag:99.99,zerr:9.99,zrms:99.99,ndetz:0,$
              ymag:99.99,yerr:9.99,yrms:99.99,ndety:0,$
              chi:99.99,sharp:99.99}
obj = replicate(obj_schema,nobj)
obj.objid = oindex.value
otags = tag_names(obj)
hasforced = tag_exist(meas,'FORCED')

;; Exposure names
;;   F23-00912982_15, trim off ccdnum
exposure = reform((strsplitter(meas.exposure,'_',/extract))[0,*])
eindex = create_index(exposure)
nexposure = n_elements(eindex.value)
;; Create exposure structure
expstr = replicate({exposure:'',filter:''},nexposure)
expstr.exposure = eindex.value
expstr.filter = meas[eindex.index[eindex.lo]].filter
;print,strtrim(nexposure,2),' exposures with measurements'
ufilter = expstr[uniq(expstr.filter,sort(expstr.filter))].filter
nufilter = n_elements(ufilter)
;print,strtrim(nufilter,2),' unique filters'

;; Loop over unique filters
For f=0,nufilter-1 do begin
  filtind = where(expstr.filter eq ufilter[f],nfiltind)
  
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
    ;; Depthflag
    ;;  OR combine to depthflag, 1-single-exposure, 2-forced, 3-both
    ;;  within a single exposure we can have some forced and some
    ;;  not, need to do this object by object
    ;; forced=0 -> depthflag=1, forced=1 -> depthflag=2, depthflag=forced+1
    if hasforced eq 1 then $
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
      if hasforced eq 1 then $
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

;; NDET
obj.ndet = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz+obj.ndety


;; Coordinates, morphology, ebv, etc
;;-----------------------------------

;; Get average chi, sharp
totchi = fltarr(nobj) & numchi = lon64arr(nobj)
totsharp = fltarr(nobj) & numsharp = lon64arr(nobj)

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
    totsharp[gdsharp] += sharp1[gdsharp]
    numsharp[gdsharp]++
  endif
Endfor
;; Make average CHI
gdchi = where(numchi gt 0,ngdchi)
avgchi = fltarr(nobj)+99.99
if ngdchi gt 0 then avgchi[gdchi]=totchi[gdchi]/numchi[gdchi]
obj.chi = avgchi
;; Make average SHARP
gdsharp = where(numsharp gt 0,ngdsharp)
avgsharp = fltarr(nobj)+99.99
if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/numsharp[gdsharp]
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
newrarms = sqrt( totalradiff/(obj.ndet>1) ) * 3600 * cos(obj.dec/!radeg)
newdecrms = sqrt( totaldecdiff/(obj.ndet>1) ) * 3600
bd = where(obj.ndet eq 0,nbd)
if nbd gt 0 then newrarms[bd]=99.99
if nbd gt 0 then newdecrms[bd]=99.99
; Set rms=99.99 for numobs=1
oneobs = where(obj.ndet eq 1,noneobs)
if noneobs gt 0 then newrarms[oneobs]=99.99
if noneobs gt 0 then newdecrms[oneobs]=99.99
obj.rarms = newrarms
obj.decrms = newdecrms

end
