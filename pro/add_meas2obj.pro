pro add_meas2obj,obj,meas

;; Add measurements to OBJ, update OBJID in MEAS
;; The structures should be matched already
;; MEAS should be of a single exposure and band

nobj = n_elements(obj)
nmeas = n_elements(meas)
if nobj ne nmeas then begin
  print,'OBJ and MEAS just have same number of elements'
  return
endif

ui = uniq(meas.filter,sort(meas.filter))
if n_elements(ui) ne 1 then begin
  print,'MEAS must only be one exposure'
  return
endif
filter = meas[0].filter
otags = tag_names(obj)
magind = where(otags eq strupcase(filter)+'MAG',nmagind)
if nmagind eq 0 then begin
  print,'NO ',strupcase(filter)+'MAG column'
  return
endif
errind = where(otags eq strupcase(filter)+'ERR',nerrind)
if nerrind eq 0 then begin
  print,'NO ',strupcase(filter)+'ERR column'
  return
endif
scatind = where(otags eq strupcase(filter)+'SCATTER',nscatind)
if nscatind eq 0 then begin
  print,'NO ',strupcase(filter)+'SCATTER column'
  return
endif
detind = where(otags eq 'NDET'+strupcase(filter),ndetind)
if ndetind eq 0 then begin
  print,'NO NDET'+strupcase(filter)+' column'
  return
endif

origobj = obj

;; Leave RA/DEC unchanged


ndet = obj.(detind)
ndetall = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz+obj.ndety
nomeas = where(obj.(magind) gt 50,nnomeas,comp=hasmeas,ncomp=nhasmeas)
oldmag = double(obj.(magind))

;; Objects with no good measurements in this band
if nnomeas gt 0 then begin
  obj[nomeas].(magind) = meas[nomeas].mag
  obj[nomeas].(errind) = meas[nomeas].err
  obj[nomeas].(scatind) = 99.99
  obj[nomeas].(detind) = 1
  ;; Update chi/sharp
  ;; makemag.pro takes the median, we'll just average this in
  chisum = obj[nomeas].chi * ndetall[nomeas]
  newchi = (chisum+meas[nomeas].chi)/float(ndetall[nomeas]+1)
  obj[nomeas].chi = newchi
  sharpsum = obj[nomeas].sharp * ndetall[nomeas]
  newsharp = (sharpsum+meas[nomeas].sharp)/float(ndetall[nomeas]+1)
  obj[nomeas].sharp = newsharp
endif
;; Objects that already have some good measurements
if nhasmeas gt 0 then begin

  ;; ONLY TAKE BEST MAGNITUDES??

  ;; Take flux-weighed mean
  totalwt = 1/obj[hasmeas].(errind)^2
  totflux = (10d0)^(obj[hasmeas].(magind)/2.5d0)
  totalfluxwt = totflux*totalwt
  ;; add new measurements
  totalwt += 1.0d0/(meas[hasmeas].err)^2
  totalfluxwt += 2.5118864d^meas[hasmeas].mag * (1.0d0/(meas[hasmeas].err)^2)
  ;; get new mags/errs
  newflux = totalfluxwt/totalwt
  newmag = 2.50*alog10(newflux)
  newerr = sqrt(1.0/totalwt)
  obj[hasmeas].(magind) = newmag
  obj[hasmeas].(errind) = newerr
  obj[hasmeas].(detind) += 1
  ;; Update scatter
  ;; Measure scatter, RMS
  ;;  rms = sqrt(mean(diff^2))
  ;;      = sqrt(mean(x^2)-<x>^2)
  ;; no previous scatter
  hasmeas1 = where(obj.(magind) lt 50 and ndet eq 1,nhasmeas1)
  if nhasmeas1 gt 0 then begin
    newrms = sqrt( ( double(oldmag[hasmeas1]-obj[hasmeas1].(magind))^2 + $
                     double(meas[hasmeas1].mag-obj[hasmeas1].(magind))^2 )/2.0 )
    obj[hasmeas1].(scatind) = newrms
  endif
  ;; a previous scatter
  hasmeas2 = where(obj.(magind) lt 50 and ndet gt 1,nhasmeas2)
  if nhasmeas2 gt 0 then begin
    ;; rms = sqrt(mean(x^2)-<x>^2)
    meanx2 = ( (double(obj[hasmeas2].(scatind)))^2 + oldmag[hasmeas2]^2 )*ndet[hasmeas2]
    ;; add new measurement
    meanx2 += (double(meas[hasmeas2].mag))^2
    ;; true mean magnitude of previous value and this new magnitude
    meanmag = (oldmag[hasmeas2]*ndet[hasmeas2]+meas[hasmeas2].mag)/double(ndet[hasmeas2]+1)
    ;; rms2 = 1/n*sum(x^2) - 2*<x>*xwt + xwt^2
    ;; need the true mean in there, even if we want to use a new midpoint/wt
    newrms = sqrt( meanx2/double(ndet[hasmeas2]+1) - 2*meanmag*double(obj[hasmeas2].(magind)) + double(obj[hasmeas2].(magind))^2 )
    obj[hasmeas2].(scatind) = newrms
  endif

  ;; Update chi/sharp
  ;; makemag.pro takes the median, we'll just average this in
  chisum = obj[hasmeas].chi * ndetall[hasmeas]
  newchi = (chisum+meas[hasmeas].chi)/float(ndetall[hasmeas]+1)
  obj[hasmeas].chi = newchi
  sharpsum = obj[hasmeas].sharp * ndetall[hasmeas]
  newsharp = (sharpsum+meas[hasmeas].sharp)/float(ndetall[hasmeas]+1)
  obj[hasmeas].sharp = newsharp
endif

;stop

end
