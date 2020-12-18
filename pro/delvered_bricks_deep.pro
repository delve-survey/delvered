;+
;
; DELVERED_BRICKS_DEEP
;
; This is a script to run PHOTRED ALLFRAME on DELVE MC exposures to create
; forced-photometry catalogs.
;
; INPUTS:
;  input    What bricks to run PHOTRED/ALLFRAME on.  Either an array or
;            a range such as 100-110.
;  =nmulti  Number of simultaneous jobs to run.
;  /redo    Redo bricks that were previously done.
;  /update  Check for updates and rerun if there are.
;  =delvedir  The main delve directory.  Default is /net/dl2/dnidever/delve/
;  =delvereddir  The main delvered software directory.  Default is /home/dnidever/projects/delvered/'.
;
; OUTPUTS:
;  PHOTRED_ALLFRAME will be run on each brick and a final catalog and
;  summary file created for each brick.
;
; USAGE:
;  IDL>delvered_bricks,100
;
; By D. Nidever  Aug 2019
;-

pro delvered_bricks_deep,input,nmulti=nmulti,redo=redo,update=update,delvedir=delvedir,delvereddir=delvereddir,stp=stp

t0 = systime(1)
  
;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/deep/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
tempdir = '/tmp/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir = '/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
if n_elements(nmulti) eq 0 then nmulti=10

;; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - delvered_bricks,input,nmulti=nmulti,redo=redo,update=update,stp=stp'
  return
endif

; Start the logfile
;------------------
; format is delvered_bricks.host.DATETIME.log
host = GETENV('HOST')
hostname = first_el(strsplit(host,'.',/extract))
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logtime = smonth+sday+syear+shour+sminute+ssecond
logfile = logsdir+'delvered_bricks.'+hostname+'.'+logtime+'.log'
JOURNAL,logfile

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delve_deep_bricks.fits',1)

;; Parse the input nights
for i=0,n_elements(input)-1 do begin
  input1 = input[i]
  ;; See if there is a -
  ;; Use RA-DEC coordinate as a BOX selection
  ;;   could also think of them as an ordered list and take all bricks
  ;;   between the two (inclusive)
  if strpos(input1,'-') ne -1 then begin
    brickrange = strsplit(input1,'-',/extract)
    br = float(strsplitter(brickrange,'pm',/extract)) / 10.  ;; RA and DEC components
    br = transpose(br)
    ind_bricks = where(brickstr.ra ge br[0,0] and brickstr.ra le br[0,1] and $
                       brickstr.dec ge br[0,1] and brickstr.dec le br[1,1],nind_bricks)
    if nind_bricks gt 0 then push,bricks,brickstr[ind_bricks] else print,'No brick found matching ',input1
  endif else begin
    if input1[0] eq '*' then push,bricks,brickstr.brickname else begin
      MATCH,brickstr.brickname,input1,ind1,ind2,/sort,count=nmatch
      if nmatch gt 0 then push,bricks,input1 else print,input1,' brick not found'
    endelse
  endelse
endfor
nbricks = n_elements(bricks)
if nbricks eq 0 then begin
  print,'No bricks to process'
  return
endif
bricks = strtrim(bricks,2)
bricks = bricks[uniq(bricks,sort(bricks))]
nbricks = n_elements(bricks)

; Print info
;-----------
print,''
print,'############################################'
print,'Starting DELVERED_BRICKS   ',systime(0)
print,'Running on ',host
print,'############################################'
print,''


; Make sure we have the right printlog.pro, not Markwardt's version
tempprogs = strsplit(!path,':',/extract)+'/printlog.pro'
test = file_test(tempprogs)
ind = where(test eq 1,nind)
bd = where(stregex(tempprogs[ind],'markwardt',/boolean) eq 1,nbd)
if nbd gt 0 then begin
  baddir = file_dirname(tempprogs[ind[bd]])
  print,"There is a version of Markwardt's PRINTLOG.PRO in "+baddir
  print,'Please rename this program (i.e. printlog.pro.orig)'
  return
endif

;; Check if bricks were previously done
if not keyword_set(redo) and not keyword_set(update) then begin
  print,'Checking for any bricks were previously done'
  outfiles = brickdir+bricks+'/'+bricks+'.fits.gz'
  bd = where(file_test(outfiles) eq 1,nbd,comp=gd,ncomp=ngd)
  if ngd eq 0 then begin
    print,'All bricks were previously processed.  Nothing to do.'
    return
  endif
  if nbd gt 0 then begin
    print,'REDO not set and ',strtrim(nbd,2),' bricks were previously processed. Removing them from current list to process.'
    bricks = bricks[gd]
    nbricks = ngd
  endif
endif

;; Check if we need to update the exposures database
dbfile = '/net/dl2/dnidever/delve/deep/bricks/db/delvered_summary.db'
lockfile = dbfile+'.lock'
print,'Checking if exposures database needs to be updated'
;; Check for lockfile
while (file_test(lockfile) eq 1) do begin
  print,'Lock file found. Waiting 60 seconds.'
  wait,60
endwhile
dbinfo = file_info(dbfile)
nightdirs = file_search(delvedir+'exposures/20??????',count=nnightdirs)
sumfiles = nightdirs+'/'+file_basename(nightdirs)+'_summary.fits'
suminfo = file_info(sumfiles)  
gsum = where(suminfo.exists eq 1 and suminfo.size gt 0,ngsum)
if ngsum gt 0 then begin
  if max(suminfo[gsum].mtime) gt dbinfo.mtime then begin
    print,'Need to update the database'
    touchzero,lockfile 
    print,'Waiting 10 sec to let current queries finish'
    wait,10
    ;; This takes about 2.5 min to run
    SPAWN,delvereddir+'bin/make_delvered_deep_summary_table',/noshell
    file_delete,lockfile,/allow
  endif
endif

print,'Processing ',strtrim(nbricks,2),' brick(s)'

;#########################################
;#  STARTING THE PROCESSING
;#########################################
cmd = "delvered_forcebrick_deep,'"+bricks+"',delvedir='"+delvedir+"'"
if keyword_set(redo) then cmd += ',/redo'
if keyword_set(update) then cmd += ',/update'
cmddirs = strarr(nbricks)+workdir
PBS_DAEMON,cmd,cmddirs,jobs=jobs,prefix='dlvbrcks',/idle,/hyperthread,nmulti=nmulti,wait=5

print,'DELVERED_BRICKS done after ',strtrim(systime(1)-t0,2),' sec'

; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
