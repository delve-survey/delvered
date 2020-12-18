;+
;
; DELVERED_NIGHTSUMMARY
;
; This creates a summary file of all the taken on one night.
;
; INPUTS:
;  night    The night name, e.g. 20160101.
;
; OUTPUTS:
;  A night summary file will be created in the night directory
;  with name NIGHT_summary.fits.
;
; USAGE:
;  IDL>delvered_nightsummary,'20160101'
;
; By D. Nidever  Aug 2019
;-

pro delvered_nightsummary,night,delvedir=delvedir,redo=redo,stp=stp

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'

;; Not enough inputs
if n_elements(night) eq 0 then begin
  print,'Syntax - delvered_nightsummary,night,delvedir=delvedir,redo=redo,stp=stp'
  return
endif
inight = strtrim(night,2)

CD,current=origdir

;; Check if directory exists
if file_test(expdir+inight) eq 0 then begin
  print,'Directory '+expdir+inight+' NOT FOUND'
  return
endif
;; Check if output file already exists
nightsumfile = expdir+inight+'/'+inight+'_summary.fits'
info = file_info(nightsumfile)
if info.exists eq 1 then nexp=sxpar(headfits(nightsumfile),'naxis2') else nexp=0
if (info.exists eq 1 and nexp gt 0) and not keyword_set(redo) then begin
  print,nightsumfile+' EXISTS and /redo NOT set'
  return
endif

CD,expdir+inight

;; Create the nightly summary file
sumfiles = FILE_SEARCH('*_summary.fits',count=nsumfiles)
undefine,expstr,chipstr
fields = IMPORTASCII('fields',fieldnames=['shname','name'],fieldtypes=[7,7],/silent)
nfields = n_elements(fields)
print,strtrim(nfields,2),' fields'
For j=0,nfields-1 do begin
  ;print,strtrim(j+1,2),' ',fields[j].name
  sumfile = fields[j].name+'_summary.fits'
  if file_test(sumfile) eq 1 then begin
    name = fields[j].name
    shname = fields[j].shname
    expstr1 = MRDFITS(sumfile,1,/silent)
    add_tag,expstr1,'fieldname',name,expstr1
    add_tag,expstr1,'field',shname,expstr1
    PUSH,expstr,expstr1
    chipstr1 = MRDFITS(sumfile,2,/silent)
    add_tag,chipstr1,'fieldname',name,chipstr1
    PUSH,chipstr,chipstr1
  endif else print,sumfile,' NOT FOUND'
Endfor
if n_elements(expstr) eq 0 then begin
  print,'No exposures for night=',inight
endif
print,'Writing nightly summary file to ',nightsumfile
MWRFITS,expstr,nightsumfile,/create
MWRFITS,chipstr,nightsumfile,/silent

if keyword_set(stp) then stop

end
