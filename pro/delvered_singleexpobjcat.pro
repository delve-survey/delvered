pro delvered_singleexpobjcat,brick,meas,obj

;; Create measurement and object catalogs from the single-exposure
;; level processing ALS files.

;brick = '0933m587'
dir = '/net/dl2/dnidever/delve/bricks/'+strmid(brick,0,4)+'/'+brick+'/'
meta = mrdfits(dir+brick+'_joint_meta.fits',1)
nmeta = n_elements(meta)

;; Add ALF filename
add_tag,meta,'alffile','',meta
meta.alffile = dir+strtrim(meta.base,2)+'.alf'


;; Measurement schema
meas_schema = {id:'',objid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,forced:0B,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,$
               imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0}
mtags = tag_names(meas_schema)

;; Initialize the final object table, with ALL BANDS
obj_schema = {objid:'',forced:0B,ra:0.0d0,dec:0.0d0,ndet:0L,$
              umag:99.99,uerr:9.99,uscatter:99.99,ndetu:0L,$
              gmag:99.99,gerr:9.99,gscatter:99.99,ndetg:0L,$
              rmag:99.99,rerr:9.99,rscatter:99.99,ndetr:0L,$
              imag:99.99,ierr:9.99,iscatter:99.99,ndeti:0L,$
              zmag:99.99,zerr:9.99,zscatter:99.99,ndetz:0L,$
              ymag:99.99,yerr:9.99,yscatter:99.99,ndety:0L,$
              chi:99.99,sharp:99.99,depthflag:0,brickuniq:0B}
otags = tag_names(obj_schema)

obj = replicate(obj_schema,1e6)
ocount = 0LL
meas = replicate(meas_schema,1e6)
mcount = 0LL

print,strtrim(nmeta,2),' chip files'

;; Loop over the chips
for i=0,nmeta-1 do begin
  print,strtrim(i+1,2),' ',meta[i].base
  meta1 = meta[i]
  usealf = 0  ; 1
  cat1 = LOAD_CHIPCAT(meta1,usealf=usealf)
  ncat1 = n_elements(cat1)
  if size(cat1,/type) ne 8 then goto,BOMB
  
  ;; Apply S/N cut
  ;snr = 1.087/cat1.err
  ;snrthresh = 5.0
  ;bd = where(snr lt snrthresh,nbd,comp=gd,ncomp=ngd)
  ;print,'Applying S/N>'+stringize(snrthresh,ndec=2)+' cut'
  ;if ngd eq 0 then goto,BOMB
  ;REMOVE,bd,cat1
  ;ncat1 = n_elements(cat1)

  magind = where(otags eq strupcase(meta1.filter)+'MAG',nmagind)
  errind = where(otags eq strupcase(meta1.filter)+'ERR',nerrind)
  detind = where(otags eq 'NDET'+strupcase(meta1.filter),ndetind)

  ;; First catalog
  if ocount eq 0 then begin
    ;; Update the object catalog
    print,'Adding ',strtrim(ncat1,2),' objects'
    newobj = replicate(obj_schema,ncat1)
    struct_assign,cat1,newobj,/nozero
    ;newobj.(magind) = cat1.mag
    ;newobj.(errind) = cat1.err
    ;newobj.(detind) = 1
    newobj.ndet = 1
    newobj.objid = strtrim(brick,2)+'.'+strtrim(lindgen(ncat1)+1,2)
    obj[0:ncat1-1] = newobj
    ocount += ncat1
    ;; Update the measurement catalog
    cat1.objid = newobj.objid
    meas[0:ncat1-1] = cat1
    mcount += ncat1

  ;; Second and up
  endif else begin
   ;; Cross-matching to existing objects
   SRCMATCH,obj[0:ocount-1].ra,obj[0:ocount-1].dec,cat1.ra,cat1.dec,1.0,ind1,ind2,/sph,count=nmatch
   print,strtrim(nmatch,2),' matches'
   if nmatch gt 0 then begin
     obj[ind1].ndet += 1
     cat1[ind2].objid = obj[ind1].objid
     if n_elements(meas) lt mcount+ncat1 then meas=add_elements(meas,ncat1>1e5)
     meas[mcount:mcount+nmatch-1] = cat1[ind2]
     mcount += nmatch
     if nmatch eq ncat1 then begin
       undefine,cat1
       ncat1 = 0
     endif else begin
       REMOVE,ind2,cat1
       ncat1 = n_elements(cat1)
     endelse
   endif
   ;; Some left, add them
   if ncat1 gt 0 then begin
      ;; Update the object catalog
      print,'Adding ',strtrim(ncat1,2),' objects'
      newobj = replicate(obj_schema,ncat1)
      struct_assign,cat1,newobj,/nozero
      newobj.ndet = 1
      newobj.objid = strtrim(brick,2)+'.'+strtrim(lindgen(ncat1)+1+ocount,2)
      if n_elements(obj) lt ocount+ncat1 then obj=add_elements(obj,ncat1>1e5)
      obj[ocount:ocount+ncat1-1] = newobj
      ocount += ncat1
      ;; Update the measurement catalog
      if n_elements(meas) lt mcount+ncat1 then meas=add_elements(meas,ncat1>1e5)
      cat1.objid = newobj.objid
      meas[mcount:mcount+ncat1-1] = cat1
      mcount += ncat1
   endif
  endelse
  BOMB:
endfor  ;; file loop
;; Trim extra elements
if n_elements(obj) gt ocount then obj=obj[0:ocount-1]
if n_elements(meas) gt mcount then meas=meas[0:mcount-1]

;; Exposure information
ui = uniq(meta.expnum,sort(meta.expnum))
expstr = meta[ui]
add_tag,expstr,'exposure','',expstr
expstr.exposure = reform((strsplitter(expstr.base,'_',/extract))[0,*])
si = sort(expstr.exposure)   ;; sort by exposure
expstr = expstr[si]
exposure = reform((strsplitter(meas.exposure,'_',/extract))[0,*])
eindex = create_index(exposure)
match,expstr.exposure,eindex.value,ind1,ind2,/sort,count=nmatch
expstr[ind1].nmeas = eindex.num[ind2]

;; Average the photometry
print,'Averaging the photometry'
oldobj = obj
DELVERED_AVGMEAS,expstr,meas,obj

;; Calculate photometric variability metrics
print,'Calculating photometric variability metrics'
DELVERED_PHOTVAR,meas,obj

;stop

end
