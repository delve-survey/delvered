;+
;
; DELVERED_PHOTVAR
;
; A program to calculate photometric variability indices.
;
; INPUTS:
;  meas   Structure of all of the measurements.
;  obj    Structure for the objects.  Variability metric columsn must
;           already exist.
;
; OUTPUTS:
;  The variability metric columns in OBJ are updated.
;
; USAGE:
;  IDL>delvered_photvar,meas,obj
;
; By D.Nidever  copied from NSC DR2 combination code
;-

pro delvered_photvar,meas,obj

;; Not enough inputs
nobj = n_elements(obj)
nmeas = n_elements(meas)
if nmeas eq 0 or nobj eq 0 then begin
  print,'Syntax - delvered_photvar,meas,obj'
  return
endif
print,strtrim(nmeas,2),' measurements'

otags = tag_names(obj)
mtags = tag_names(meas)

;; Create index on OBJID
oindex = create_index(meas.objid)
nobj = n_elements(oindex.value)
print,strtrim(nobj,2),' objects'

;; All variability metrics are bad to start with
obj.rmsvar = !values.f_nan
obj.madvar = !values.f_nan
obj.iqrvar = !values.f_nan
obj.etavar = !values.f_nan
obj.jvar = !values.f_nan
obj.kvar = !values.f_nan
obj.chivar = !values.f_nan
obj.romsvar = !values.f_nan


;; Loop over the objects
fidmag = fltarr(nobj)+!values.f_nan  ; fiducial magnitude
for i=0L,nobj-1 do begin
  ind = oindex.index[oindex.lo[i]:oindex.hi[i]]
  nmeas = n_elements(ind)
  meas1 = meas[ind]

  filtindex = create_index(meas1.filter)
  nfilters = n_elements(filtindex.value)
  resid = dblarr(nmeas)+!values.f_nan     ; residual mag
  relresid = dblarr(nmeas)+!values.f_nan  ; residual mag relative to the uncertainty
  for f=0,nfilters-1 do begin
    filt = strupcase(filtindex.value[f])
    findx = filtindex.index[filtindex.lo[f]:filtindex.hi[f]]
    magind = where(otags eq filt+'MAG',nmagind)
    errind = where(otags eq filt+'ERR',nerrind)
    gph = where(meas1[findx].mag lt 50,ngph)
    if ngph gt 1 then begin
      ;; Residual mag
      resid[findx[gph]] = meas1[findx[gph]].mag-obj[i].(magind)
      ;; Residual mag relative to the uncertainty
      ;;  set a lower threshold of 0.02 in the uncertainty
      relresid[findx[gph]] = sqrt(ngph/(ngph-1.0)) * (meas[findx[gph]].mag-obj[i].(magind))/(meas[findx[gph]].err > 0.02)
    endif
  endfor ; filter loop

  ;; Calculate variability indices
  gdresid = where(finite(resid) eq 1,ngdresid)
  if ngdresid gt 0 then begin
    resid2 = resid[gdresid]
    sumresidsq = total(resid2^2)
    tsi = sort(meas1[gdresid].mjd)
    resid2tsi = resid2[tsi]
    quartiles = cgpercentiles(resid2,percentiles=[0.25,0.50,0.75])
    ;; RMS
    rms = sqrt(sumresidsq/ngdresid)
    ;; MAD
    madvar = 1.4826*median(abs(resid2-quartiles[1]))
    ;; IQR
    iqrvar = 0.741289*(quartiles[2]-quartiles[0])
    ;; 1/eta
    etavar = sumresidsq / total((resid2tsi[1:*]-resid2tsi[0:ngdresid-2])^2)
    obj[i].rmsvar = rms
    obj[i].madvar = madvar
    obj[i].iqrvar = iqrvar
    obj[i].etavar = etavar
  endif

  ;; Calculate variability indices wrt to uncertainties
  gdrelresid = where(finite(relresid) eq 1,ngdrelresid)
  if ngdrelresid gt 0 then begin
    relresid2 = relresid[gdrelresid]
    pk = relresid2^2-1
    jvar = total( sign(pk)*sqrt(abs(pk)) )/ngdrelresid
    ;avgrelvar = np.mean(np.abs(relresid2))    ; average of absolute relative residuals
    chivar = sqrt(total(relresid2^2))/ngdrelresid
    kdenom = sqrt(total(relresid2^2)/ngdrelresid)
    if kdenom ne 0 then begin
      kvar = (total(abs(relresid2))/ngdrelresid) / kdenom
    endif else begin
      kvar = !values.f_nan
    endelse
    ;; RoMS
    romsvar = total(abs(relresid2))/(ngdrelresid-1)
    obj[i].jvar = jvar
    obj[i].kvar = kvar
    ;;obj['avgrelvar'] = avgrelvar
    obj[i].chivar = chivar
    obj[i].romsvar = romsvar
  endif


  ;; Fiducial magnitude, used to select variables below
  ;;  order of priority: r,g,i,z,Y,VR,u
  if obj[i].ndet gt 0 then begin
    magarr = [obj[i].rmag,obj[i].gmag,obj[i].imag,obj[i].zmag,obj[i].ymag,obj[i].umag]
    gfid = where(magarr lt 50,ngfid)
    if ngfid gt 0 then fidmag[i]=magarr[gfid[0]]
  endif
endfor  ; object loop


;; Select variables using photometric variability indices.
    
;; Select Variables
;;  1) Construct fiducial magnitude (done in loop above)
;;  2) Construct median VAR and sigma VAR versus magnitude
;;  3) Find objects that Nsigma above the median VAR line
si = sort(fidmag)   ; NaNs are at end
varcol = 'madvar'
gdvar = where(finite(obj.madvar) eq 1 and finite(fidmag) eq 1,ngdvar,comp=bdvar,ncomp=nbdvar)
obj.variable10sig = 0
if ngdvar gt 0 then begin
  binsize = 0.25
  BINDATA,fidmag[gdvar],fidmag[gdvar],fidmagmed_bins,fidmagmed,binsize=binsize,/med
  BINDATA,fidmag[gdvar],fidmag[gdvar],xbin2,numhist,binsize=binsize,/hist
  ;; Fix NaNs in fidmagmed
  bdfidmagmed = where(finite(fidmagmed) eq 0,nbdfidmagmed)
  if nbdfidmagmed gt 0 then fidmagmed[bdfidmagmed] = fidmagmed_bins[bdfidmagmed]
  ;; Median metric
  BINDATA,fidmag[gdvar],obj[gdvar].madvar,xbin,varmed,binsize=binsize,/med
  ;; Smooth, it handles NaNs well
  smlen = 5
  smvarmed = gsmooth(varmed,smlen)
  bdsmvarmed = where(finite(smvarmed) eq 0,nbdsmvarmed)
  if nbdsmvarmed gt 0 then smvarmed[bdsmvarmed] = median(smvarmed)
  ;; Interpolate to all the objects
  gv = where(finite(smvarmed) eq 1,ngv,comp=bv,ncomp=nbv)
  if ngv le 1 then begin
    print,'Not enough good MADVAR values to detect variables'
    obj.variable10sig = 0
    obj.nsigvar = !values.f_nan
    return
  endif
  INTERP,fidmagmed[gv],smvarmed[gv],fidmag[gdvar],objvarmedgd
  objvarmed = dblarr(nobj)
  objvarmed[gdvar] = objvarmedgd
  objvarmed[gdvar] = min(smvarmed[gv]) > objvarmed[gdvar]   ; lower limit
  if nbdvar gt 0 then objvarmed[bdvar]=smvarmed[gv[ngv-1]]   ; objects with bad fidmag, set to last value
  ;; Scatter in metric around median
  ;;  calculate MAD ourselves so that it's around our computed median metric line
  BINDATA,fidmag[gdvar],abs(obj[gdvar].madvar-objvarmed[gdvar]),xbin3,varsig,binsize=binsize,/med
  varsig *= 1.4826   ; scale MAD to stddev
  ;; Fix values for bins with few points
  bdhist = where(numhist lt 3,nbdhist,comp=gdhist,ncomp=ngdhist)
  if nbdhist gt 0 then begin
    if ngdhist gt 0 then begin
      varsig[bdhist] = median(varsig[gdhist])
    endif else begin
      varsig[*] = 0.02
    endelse
  endif
  ;; Smooth
  smvarsig = gsmooth(varsig,smlen)
  ;; Interpolate to all the objects
  gv = where(finite(smvarsig) eq 1,ngv,comp=bv,ncomp=nbv)
  if ngv le 1 then begin
    print,'Not enough good MADVAR values to detect variables'
    obj.variable10sig = 0
    obj.nsigvar = !values.f_nan
    return
  endif
  INTERP,fidmagmed[gv],smvarsig[gv],fidmag[gdvar],objvarsiggd
  objvarsig = dblarr(nobj)
  objvarsig[gdvar] = objvarsiggd
  objvarsig[gdvar] = min(smvarsig[gv]) > objvarsig[gdvar]   ; lower limit
  if nbdvar gt 0 then objvarsig[bdvar]=smvarsig[gv[ngv-1]]   ; objects with bad fidmag, set to last value
  ;; Detect positive outliers
  nsigvarthresh = 10.0
  nsigvar = (obj.madvar-objvarmed)/objvarsig
  obj[gdvar].nsigvar = nsigvar[gdvar]
  isvar = where(nsigvar[gdvar] gt nsigvarthresh,nisvar)
  print,strtrim(nisvar,2)+' variables detected'
  if nisvar gt 0 then obj[gdvar[isvar]].variable10sig = 1
endif

end
