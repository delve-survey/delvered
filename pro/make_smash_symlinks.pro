pro make_smash_symlinks,redo=redo

;; Make symlinks to the SMASH data that has already been processedd

delvedir = '/net/dl1/users/dnidever/delve/exposures/'
smashdir = '/net/dl1/users/dnidever/smash/cp/red/photred/'
nmulti = 3 ;5

nights = FILE_SEARCH(smashdir+'20??????',/test_directory,count=nnights)
CD,current=origdir

scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
irafdir = '/home/dnidever/iraf/'
workdir = '/data0/dnidever/delve/'

print,'######################################'
print,strtrim(nnights,2),' SMASH nights to create symlinks for'
print,'######################################'

;; Night loop
undefine,cmd,cmddir
For i=0,nnights-1 do begin
  inight = FILE_BASENAME(nights[i])
  ;; Skip if the directory exists
  ;if file_test(delvedir+inight,/directory) eq 0 then begin
  ;if inight ne '20160101' and inight ne '20121116' then begin
    cmd1 = "make_smash_symlinks_night,'"+inight+"'"
    if keyword_set(redo) then cmd1+=",/redo"
    push,cmd,cmd1
  ;endif else print,nights[i],' exists already.  Skipping.'
Endfor
ncmd = n_elements(cmd)
cmddir = strarr(ncmd)+'/data0/dnidever/delve/'

stop

;; Run PBS_DAEMON
print,'Starting the JOBS'
PBS_DAEMON,cmd,cmddir,jobs=jobs,nmulti=nmulti,/hyperthread,/idle,prefix='smash',wait=5

stop


end
