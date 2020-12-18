pro dlvexpstat,delvedir=delvedir,redo=redo,all=all

;; DELVERED exposure status

;; Defaults
;if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'

;; Get all of the nights
dirs = FILE_SEARCH(expdir+'20??????',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)
print,strtrim(ndirs,2),' DELVE nights'

schema = {night:'',dir_times:lon64arr(3)-1,dir_size:-1L,logsdir:-1,logsdir_times:lon64arr(3)-1,logsdir_size:-1L,fieldsfile:-1,$
          setupfile:-1,exposurefile:-1,nightsumfile:-1,wcs_success:-1L,wcs_failure:0L,$
          daophot_success:-1L,daophot_failure:0L,match_success:-1L,match_failure:0L,apcor_success:-1L,apcor_failure:0L,$
          astrom_success:-1L,astrom_failure:0L,zeropoint_success:-1L,zeropoint_failure:0L,calib_success:-1L,calib_failure:0L,$
          combine_success:-1L,combine_failure:0L,deredden_success:-1L,deredden_failure:0L,save_success:-1L,save_failure:0L,$
          success:lonarr(10)-1,failure:lonarr(10),nexpall:-1L,nexp:-1L,done:0B}
tags = tag_names(schema)

;; Start with last summary
sumfiles = file_search(expdir+'summary/delve_expsummary_*.fits',count=nsumfiles)
last = replicate(schema,ndirs)
last.night = dirs
if nsumfiles gt 0 then begin
  si = reverse(sort(sumfiles))
  str0 = mrdfits(sumfiles[si[0]],1,/silent)
  nstr0 = n_elements(str0)
  match,last.night,str0.night,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    ;; start with the last one
    temp = last[ind1]
    struct_assign,str0[ind2],temp,/nozero
    last[ind1] = temp
    ;last[ind1] = str0[ind2]
  endif
  str = last    ;; start with the last one
endif else begin
  str = replicate(schema,ndirs)
  str.night = dirs
endelse
; REDO, start fresh
if keyword_set(redo) then begin
  str = replicate(schema,ndirs)
  str.night = dirs
endif

stages = ['WCS','DAOPHOT','MATCH','APCOR','ASTROM','ZEROPOINT','CALIB','COMBINE','DEREDDEN','SAVE']
nstages = n_elements(stages)
;; Night loop
ngivesum = 0
For i=0,ndirs-1 do begin
  inight = dirs[i]
  idir = expdir+inight+'/'
  str[i].night = inight
  ;; Checking directory times and size against previous one
  check = 1
  dir_info = file_info(expdir+inight)
  str[i].dir_times = [dir_info.atime,dir_info.ctime,dir_info.mtime]
  str[i].dir_size = dir_info.size
  str[i].logsdir = file_test(idir+'logs/',/directory)
  ;; Check logs directory now
  logs_info = file_info(idir+'logs/')
  str[i].logsdir_times = [logs_info.atime,logs_info.ctime,logs_info.mtime]
  str[i].logsdir_size = logs_info.size
  ;; Something changed 
  if max(abs([str[i].dir_times,str[i].dir_size]-[last[i].dir_times,last[i].dir_size])) gt 0 or keyword_set(redo) or $
     max(abs([str[i].logsdir_times,str[i].logsdir_size]-[last[i].logsdir_times,last[i].logsdir_size])) gt 0 or keyword_set(all) then begin

    str[i].fieldsfile = file_test(idir+'fields')
    str[i].setupfile = file_test(idir+'photred.setup')
    str[i].exposurefile = file_test(idir+inight+'_exposures.fits')
    str[i].nightsumfile = file_test(idir+inight+'_summary.fits')
    if str[i].exposurefile eq 1 then str[i].nexpall=sxpar(headfits(idir+inight+'_exposures.fits',exten=1),'naxis2')
    if str[i].nightsumfile eq 1 then str[i].nexp=sxpar(headfits(idir+inight+'_summary.fits',exten=1),'naxis2')
    if str[i].logsdir eq 1 then begin
      for j=0,nstages-1 do begin
        sind = where(tags eq stages[j]+'_SUCCESS',nsind)
        sinfo = file_info(idir+'/logs/'+stages[j]+'.success')
        if sinfo.exists eq 1 and sinfo.size gt 0 then str[i].(sind[0])=file_lines(sinfo.name) else str[i].(sind[0])=0
        find = where(tags eq stages[j]+'_FAILURE',nfind)
        finfo = file_info(idir+'/logs/'+stages[j]+'.failure')
        if finfo.exists eq 1 and finfo.size gt 0 then str[i].(find[0])=file_lines(finfo.name) else str[i].(find[0])=0
        str[i].success[j] = str[i].(sind[0])
        str[i].failure[j] = str[i].(find[0])
      endfor
      ;if str[i].nexp lt 0 then str[i].nexp = max(str[i].success) / 61
    endif
  endif

  if (str[i].logsdir eq 1 and max(str[i].success) gt 0) or keyword_set(all) then begin
    if ngivesum mod 40 eq 0 then print,'NIGHT      NEXP     WCS     DAOPHOT    MATCH     APCOR   ASTROM    ZEROPT  CALIB   COMBINE  DERED   SAVE  NIGHTSUM    DONE'

    comment = ''
    done = 0
    if min(str[i].failure) eq 0 and max(str[i].failure) eq 0 and min(str[i].success) gt 0 and $
       str[i].fieldsfile eq 1 and str[i].setupfile eq 1 and str[i].exposurefile eq 1 and str[i].nightsumfile eq 1 then str[i].done=1
    if str[i].done eq 1 then comment = '       FINISHED'
    if total(str[i].failure) gt 0 then comment = '       !!!!'

    format = '(A-11,I4, I7,A1,I-4, I5,A1,I-4, I5,A1,I-4, I5,A1,I-4, I4,A1,I-4, I4,A1,I-3, I4,A1,I-4, I4,A1,I-4, I3,A1,I-3, I3,A1,I-3, I5,A-15)'
    if str[i].nexp gt 0 then nexp=str[i].nexp else nexp=str[i].nexpall
    print,inight,nexp,str[i].wcs_success,'/',str[i].wcs_failure,str[i].daophot_success,'/',str[i].daophot_failure,$
          str[i].match_success,'/',str[i].match_failure,str[i].apcor_success,'/',str[i].apcor_failure,$
          str[i].astrom_success,'/',str[i].astrom_failure,str[i].zeropoint_success,'/',str[i].zeropoint_failure,$
          str[i].calib_success,'/',str[i].calib_failure,str[i].combine_success,'/',str[i].combine_failure,$
          str[i].deredden_success,'/',str[i].deredden_failure,str[i].save_success,'/',str[i].save_failure,$
          str[i].nightsumfile,comment,format=format
    ngivesum++
  endif  
Endfor
done = where(str.done eq 1,ndone)
print,strtrim(ndone,2),'/',strtrim(ndirs,2),' finished'
if ndone gt 0 then nexpdone = long(total(str[done].nexp)) else nexpdone=0
print,strtrim(nexpdone,2),'/',strtrim(long(total(str.nexpall>0)),2),' exposures finished'

;; Save summary file
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
time = string(year,month,day,hour,minute,second,format='(I4,I02,I02,I02,I02,I02)')
outfile = expdir+'/summary/delve_expsummary_'+time+'.fits'
print,'Writing summary file to ',outfile
MWRFITS,str,outfile,/create

;stop

end
