;+
;
; DELVERED_ZEROPOINT
;
; This program calculates NSC-like exposure-level zero-points
; for DELVERED_EXPOSURES.
;
; By D. Nidever  Feb 2019
;-

pro delvered_zeropoint,redo=redo,nmulti=nmulti,modeleqnfile=modeleqnfile

COMMON photred,setup

print,''
print,'############################'
print,'RUNNING DELVERED_ZEROPOINT'
print,'############################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

logfile = -1

; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting DELVERED_ZEROPOINT  ',systime(0)

; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

;; TELESCOPE
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
;; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif
;; NMULTI
if n_elements(nmulti) eq 0 then begin
  nmulti_setup = READPAR(setup,'NMULTI')
  if nmulti_setup ne '0' and nmulti_setup ne '' and nmulti_setup ne '-1' then nmulti=long(nmulti_setup)   
  ;; Use NMULTI_ZEROPOINT if set
  nmultizp = READPAR(setup,'NMULTI_ZEROPOINT')
  if nmultizp ne '0' and nmultizp ne '' and nmultizp ne '-1' then nmulti=long(nmultizp)
  ;; Default is 5
  if n_elements(nmulti) eq 0 then nmulti=5
  nmulti = nmulti > 2       ; must be >=2
  ;; for some reason running it with nmulti=1 causes problems
endif
;; Hyperthread
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=1   ;; default is to use /hyperthread

; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif

;; Model magnitude equation file
if n_elements(modeleqnfile) eq 0 then begin
  modeleqnfile = READPAR(setup,'MODELEQNFILE')
  if (modeleqnfile eq '0' or modeleqnfile eq '' or modeleqnfile eq '-1') then begin
    printlog,logfile,'NO MODELEQNFILE FOUND.  Please add to >>photred.setup<< file'
    return
  endif
endif

; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''


;###################
; GETTING INPUTLIST
;###################

;; Get all of the files from DAOPHOT.success
;READLIST,curdir+'/logs/DAOPHOT.success',fitsfiles,/unique,/fully,setupdir=curdir,count=nfitsfiles,logfile=logfile,/silent
lists = PHOTRED_GETINPUT('ZEROPOINT','DAOPHOT.success',redo=redo,ext=['fits','fits.fz'])
if lists.ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif
fitsfiles = lists.inputlines
nfitsfiles = n_elements(fitsfiles)
;; Get the unique exposures
allbase = PHOTRED_GETFITSEXT(fitsfiles,/basename)
;; Remove the ccdnum suffix
expbase = allbase
if thisimager.namps gt 1 then $
  for i=0,nfitsfiles-1 do expbase[i] = first_el(strsplit(expbase[i],thisimager.separator,/extract))
expindex = CREATE_INDEX(expbase)
nexp = n_elements(expindex.value)
printlog,logfile,strtrim(nexp,2),' unique exposures to process'

;; Loop for other chips for these exposures that were previously
;; processed successfully
if lists.nsuccesslines gt 0 then begin
  printlog,logfile,'Looking for previous success for these exposures'
  sallbase = PHOTRED_GETFITSEXT(lists.successlines,/basename)
  sexpbase = sallbase
  if thisimager.namps gt 1 then $
    for i=0,n_elements(sallbase)-1 do sexpbase[i] = first_el(strsplit(sexpbase[i],thisimager.separator,/extract))
  ;; Loop over the exposures
  for i=0,nexp-1 do begin
    MATCH,expindex.value[i],sexpbase,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then begin
      push,allbase,sallbase[ind2]
      push,expbase,sexpbase[ind2]
      push,fitsfiles,lists.successlines[ind2]
    endif
  endfor
  ;; Some added
  if n_elements(expbase) gt n_elements(expindex.index) then begin
    ;; Make sure they are unique
    ui = uniq(allbase,sort(allbase))
    allbase = allbase[ui]
    expbase = expbase[ui]
    fitsfiles = fitsfiles[ui]
    expindex = CREATE_INDEX(expbase)
    nexp = n_elements(expindex.value)
  endif
endif

;; Load the apcor.lst file
apcor = IMPORTASCII('apcor.lst',fieldnames=['name','value'],/noprint)
add_tag,apcor,'file','',apcor
apcor.file = file_basename(apcor.name,'a.del')

;########################################
;#  PROCESSING THE FILES
;########################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

;; Run DELVERED_ZEROPOINT_EXPOSURE.PRO
cmd = "delvered_zeropoint_exposure,'"+expindex.value+"',modeleqnfile='"+modeleqnfile+"'"
if keyword_set(redo) then cmd+=',/redo'
field = reform((strsplitter(expindex.value,'-',/extract))[0,*])
cmddir = field
PBS_DAEMON,cmd,cmddir,jobs=jobs,/idle,/hyperthread,prefix='zpt',nmulti=nmulti,wait=1,scriptsdir=scriptsdir


;; Load the outputs, exposure loop
expstr = replicate({name:'',field:'',filter:'',exptime:0.0,ncat:0L,nref:0L,num:0L,zpterm:99.99,zptermerr:9.99,translines:strarr(3)},nexp)
undefine,outlist,successlist,failurelist
FOR i=0,nexp-1 do begin
  expname = expindex.value[i]
  field = first_el(strsplit(file_basename(expname),'-',/extract))  ; F1   
  ind = expindex.index[expindex.lo[i]:expindex.hi[i]]
  nind = n_elements(ind)
  expfiles = fitsfiles[ind]
  outfile = field+'/'+expname+'_zeropoint.fits'
  if file_test(outfile) eq 1 then begin
    expstr1 = MRDFITS(outfile,1,/silent)
    expstr[i] = expstr1 
    push,outlist,expstr1.name
    push,successlist,expfiles
  endif else begin
    printlog,logfile,outfile+' NOT FOUND'
    push,failurelist,expfiles
  endelse
ENDFOR

;; Check for any exposures that failed
bdexp = where(expstr.num le 0,nbdexp)
if nbdexp gt 0 then printlog,logfile,'Found '+strtrim(nbdexp,2)+' exposures with no zero-point'
for i=0,nbdexp-1 do begin
  printlog,logfile,strtrim(i+1,2)+' '+expstr[bdexp[i]].name
  ind = expindex.index[expindex.lo[bdexp[i]]:expindex.hi[bdexp[i]]]
  nind = n_elements(ind)
  expfiles = fitsfiles[ind]
  ;; See if there any other exposures for this filter
  gdexp = where(expstr.num gt 1 and expstr.filter eq expstr[bdexp[i]].filter,ngdexp)
  ;; Getting mean zero-point for this filter
  if ngdexp gt 0 then begin
    expname = expstr[bdexp[i]].name
    filter = expstr[bdexp[i]].filter
    mnzpterm = mean(expstr[gdexp].zpterm)
    mnzptermerr = mean(expstr[gdexp].zptermerr)
    expstr[bdexp[i]].zpterm = mnzpterm
    expstr[bdexp[i]].zptermerr = mnzptermerr
    expstr[bdexp[i]].num = 1
    expstr[bdexp[i]].translines = [expname+'  '+filter+'  '+filter+'-'+filter+'  '+string(-mnzpterm,format='(f7.4)')+'    0.0000    0.0000   0.0000   0.0000',$
                          '                     '+string(mnzptermerr,format='(f7.4)')+'    0.0000    0.0000   0.0000   0.0000','']
    push,outlist,expstr[bdexp[i]].name
    push,successlist,expfiles
  endif else begin
    printlog,logfile,'No good zero-points for filter='+expstr[bdexp[i]].filter
  endelse
endfor

;; Write out the transformation equations
gdexp = where(expstr.num gt 0,ngdexp,comp=bdexp,ncomp=nbdexp)
if ngdexp gt 0 then begin
  printlog,logfile,'Writing transformation equations to >>delve.trans<<'
  undefine,tlines
  ;; Check if it already exists
  if file_test('delve.trans') eq 1 and not keyword_set(redo) then READLINE,'delve.trans',tlines
  for i=0,ngdexp-1 do push,tlines,expstr[gdexp[i]].translines
  WRITELINE,'delve.trans',tlines
endif else begin
  printlog,logfile,'No good transformation equation to write out'
endelse


;##########################################
;#  UPDATING LIST FILES
;##########################################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir

printlog,logfile,'DELVERED_ZEROPOINT Finished  ',systime(0)

if keyword_set(stp) then stop

end
