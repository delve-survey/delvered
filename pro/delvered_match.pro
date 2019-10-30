;+
;
; DELVERED_MATCH
;
; This does the daomatch/daomaster portion for photred
; See daomatch_dir.pro and daomatch.pro
;
; INPUTS:
;  /redo Redo files that were already done.
;  /fake Run for artificial star tests.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The DAOMASTER MCH/RAW output files.
;
; By D.Nidever  Feb 2008
;-

pro delvered_match,redo=redo,fake=fake,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING DELVERED_MATCH'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'MATCH'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting DELVERED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getinput','photred_updatelists','photred_loadsetup',$
         'photred_getfilter','photred_getexptime','maxloc','remove','daomatch','printlog','push',$
         'matchstars','loadals','importascii','minloc','psf_gaussian','slope','mad','scale_vector',$
         'roi_cut','array_indices2','range','srcmatch','robust_linefit','mpfit','loadmch','mktemp',$
         'undefine','first_el','strsplitter','touchzero','writeline','scale','stringize','signs','stress',$
         'strep','strmult','strtrim0','sign']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that the DAOMATCH/DAOMASTER programs exist
SPAWN,'which daomatch',out,errout
daomatchfile = FILE_SEARCH(out,count=ndaomatchfile)
if (ndaomatchfile eq 0) then begin
  print,'DAOMATCH PROGRAM NOT AVAILABLE'
  return
endif
SPAWN,'which daomaster',out,errout
daomasterfile = FILE_SEARCH(out,count=ndaomasterfile)
if (ndaomasterfile eq 0) then begin
  print,'DAOMASTER PROGRAM NOT AVAILABLE'
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif


; LOAD information from the "photred.setup" file
;-----------------------------------------------
; REDO
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1
; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0
; TELESCOPE
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif
; FILTREF
filtref = READPAR(setup,'FILTREF')
filtref = strtrim(filtref,2)
if (filtref eq '0' or filtref eq '' or filtref eq '-1') then begin
  printlog,logfile,'NO REFERENCE FILTER.  Please add to >>photred.setup<< file'
  return
endif
; check if filtref is a comma-delimited list of filters in order of
; preference
filtref = strsplit(filtref,',',/extract)
nfiltref = n_elements(filtref)
; MAXSHIFT
maxshift = READPAR(setup,'MCHMAXSHIFT')
if maxshift ne '0' and maxshift ne '' and maxshift ne '-1' then $
   maxshift=float(maxshift) else undefine,maxshift
if n_elements(maxshift) gt 0 then if maxshift le 0. then undefine,maxshift
; MCHUSETILES
mchusetiles = READPAR(setup,'MCHUSETILES')
if mchusetiles eq '0' or mchusetiles eq '' or mchusetiles eq '-1' then undefine,mchusetiles
tilesep = '+'
;tilesep = '.'

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif



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
; INLIST         ALS files
; OUTLIST        MCH files
; SUCCESSLIST    ALS files


; Get input
;-----------
precursor = 'DAOPHOT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='als')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines


alsdirlist = FILE_DIRNAME(inputlines)
alsbaselist = FILE_BASENAME(inputlines)
nalsbaselist = n_elements(alsbaselist)

; Unique directories
uidirs = uniq(alsdirlist,sort(alsdirlist))
uidirs = uidirs[sort(uidirs)]
dirs = alsdirlist[uidirs]
ndirs = n_elements(dirs)


;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,systime(0)

;successarr = intarr(ninputlines)-1            ; 0-bad, 1-good
undefine,outlist,successlist,failurelist

;#############################################
; Loop through all of the unique directories
;#############################################
;; Start new field structure for DAOMATCH jobs
undefine,mstr
FOR i=0,ndirs-1 do begin

  printlog,logfile,''
  printlog,logfile,'GROUPING FILES IN '+dirs[i]

  ; CD to the directory
  CD,dirs[i]

  ; Getting the files in this directory
  gdals = where(alsdirlist eq dirs[i],ngdals)
  alsfiles = alsbaselist[gdals]

  ; Get each file's field information
  ;-----------------------------------
  sep = '-'
  dum = strsplitter(alsfiles,sep,/extract)
  fieldarr = reform(dum[0,*])

  ; Unique fields
  ;--------------
  ui = uniq(fieldarr,sort(fieldarr))
  fields = fieldarr[ui]
  nfields = n_elements(fields)


  ;############################
  ; Loop through the fields
  ;############################
  FOR j=0,nfields-1 do begin

    thisfield = fields[j]

    printlog,logfile,''
    printlog,logfile,'=============================='
    printlog,logfile,'Grouping frames for Field='+thisfield
    printlog,logfile,'=============================='

    ; Get the files for this field
    gdfield = where(fieldarr eq thisfield,ngdfield)
    base = FILE_BASENAME(alsfiles[gdfield],'.als')
    nbase = n_elements(base)


    ; If MEF then match each chip/amplifier
    ; If SINGLE CHIP then match all frames

    ; Check that there are associated FITS and PSF files
    ;---------------------------------------------------
    foundarr = intarr(nbase)+1
    for k=0,nbase-1 do begin
      fitsfile = FILE_SEARCH(base[k]+['.fits','.fits.fz'],count=nfitsfile)
      if (nfitsfile eq 0) then begin
        printlog,logfile,'NO ASSOCIATED FITS FILE FOR ',base
        PUSH,failurelist,dirs[i]+'/'+base[k]+'.als'
        foundarr[k] = 0
      endif
      psffile = FILE_SEARCH(base[k]+'.psf',count=npsffile)
      if (npsffile eq 0) then begin
        printlog,logfile,'NO ASSOCIATED PSF FILE FOR ',base
        PUSH,failurelist,dirs[i]+'/'+base[k]+'.als'
        foundarr[k] = 0
      endif
    endfor
    ; Are any left?
    gd = where(foundarr eq 1,ngd)
    if (ngd gt 0) then begin
      base = base[gd]
      nbase = n_elements(base)
    endif else begin
      printlog,logfile,'NO FILES WITH ASSOCIATED FITS FILES'
      goto,BOMB
    endelse    

    
    ;;##########################
    ;;  MAKING THE GROUPS
    ;;##########################


      ;;-----------------------
      ;; Multi-chip/amp imagers
      ;;----------------------- 
      if (thisimager.namps gt 1) then begin
       
        ; Getting all of the bases with the right extensions, e.g. _1 to _16 for Blanco+MOSAIC
        ; also getting the array of amplifier numbers
        ; Assuming that there are no leading 0s on the numbers, i.e. _1 and NOT _01
        undefine,gdbase,gdbase1,amp
        for f=1,thisimager.namps do begin
          gdbase1 = where(stregex(base,thisimager.separator+strtrim(f,2)+'$',/boolean) eq 1,ngdbase1)
          iamp = strtrim(f,2)
          ; try leading 0
          if ngdbase1 eq 0 and f lt 10 then begin
            gdbase1 = where(stregex(base,thisimager.separator+'0'+strtrim(f,2)+'$',/boolean) eq 1,ngdbase1)
            iamp = '0'+strtrim(f,2)
          endif
          if ngdbase1 gt 0 then begin
            PUSH,gdbase,gdbase1
            PUSH,amp,replicate(iamp,ngdbase1)
          endif
        endfor
        ngdbase = n_elements(gdbase)
        ; Some good ones
        if (ngdbase gt 0) then begin
          printlog,logfile,strtrim(ngdbase,2)+' of the '+strtrim(nbase,2)+' files have the correct '+$
                           thisimager.instrument+' extension'

          ; Adding the ones that didn't match to the failure list
          if (ngdbase lt nbase) then begin
            bad = lindgen(nbase)
            REMOVE,gdbase,bad
            nbad = n_elements(bad)
            printlog,logfile,strtrim(nbad,2)+' files do not have the correct extension and will be added to the FAILURE list'
            PUSH,failurelist,dirs[i]+'/'+base[bad]+'.als'
          endif

          ; ONLY keep the bases that have the correct extensions
          orig_base = base
          base = base[gdbase]
          nbase = ngdbase

        ; No good ones
        endif else begin
          printlog,logfile,'NO SPLIT '+thisimager.instrument+' FILES with the correct extensions.  Skipping this field'
          PUSH,failurelist,dirs[i]+'/'+base+'.als'           ; ALL fail
          goto,BOMB
        endelse

        ; Getting unique amplifier numbers
        ui = uniq(amp,sort(amp))
        amps = amp[ui]
        si = sort(long(amps))  ; sort them
        amps = amps[si]
        namps = n_elements(amps)
    
        ; Less amps than expected
        if (namps lt thisimager.namps) then $
          printlog,logfile,'ONLY '+strtrim(namps,2)+' AMP(S). '+strtrim(thisimager.namps,2)+' EXPECTED.'

        ; Check that we have ALL files for this GROUP
        ;--------------------------------------------
        ; Some previous successes to check
        nsuccess = lists.nsuccesslines
        if (nsuccess gt 0) then begin
          printlog,logfile,''
          printlog,logfile,'Some previous successes.  Making sure we have all files in this group'

          ; Loop through the amps
          for k=0,namps-1 do begin
            printlog,logfile,'Checking amp='+amps[k]

            ; Get previously successful amp files that exist
            successbase = FILE_BASENAME(lists.successlines,'.als')
            alstest = FILE_TEST(lists.successlines)      
            matchind = where(stregex(successbase,'^'+thisfield+'-',/boolean) eq 1 and $
                             stregex(successbase,thisimager.separator+amps[k]+'$',/boolean) eq 1 and $
                             alstest eq 1,nmatchind)

            ; Found some matches
            if (nmatchind gt 0) then begin            

              ; Check if these are already in the INLIST
              undefine,ind1,ind2,num_alreadyinlist
              MATCH,successbase[matchind],base,ind1,ind2,count=num_alreadyinlist
              num_notinlist = nmatchind - num_alreadyinlist

              ; Some not in INLIST yet
              if (num_notinlist gt 0) then begin            
                printlog,logfile,'Found '+strtrim(num_notinlist,2)+' previously successful file(s) for this '+$
                                 'group NOT YET in the INLIST.  Adding.'
                indtoadd = matchind
                if num_alreadyinlist gt 0 then REMOVE,ind1,indtoadd
                PUSH,base,successbase[indtoadd]
                PUSH,amp,replicate(amps[k],num_notinlist)

                ; Setting REDO=1 so the files in the success list will be redone.
                if not keyword_set(redo) then begin
                  printlog,logfile,'Setting REDO=1'
                  redo = 1
                endif

              endif  ; some not in inlist yet
            endif  ; some files from this group in success file
          endfor ; amp loop

          ; Make sure they are unique
          ui = uniq(base,sort(base))
          ui = ui[sort(ui)]
          base = base[ui]
          amp = amp[ui]
          nbase = n_elements(base)

          printlog,logfile,''

        endif ; some successes


        ; Getting the REFERENCE Image
        ;---------------------------------
        ;  use photred_pickreferenceframe.pro
        if not keyword_set(fake) then begin
          ; Get filters for first amp
          ind1 = where(amp eq amps[0],nind1)
          base1 = base[ind1]
          ; Get fits/fits.fz file name
          base1fits = strarr(nind1)
          for l=0,nind1-1 do $
            if FILE_TEST(base1[l]+'.fits') eq 1 then base1fits[l]=base1[l]+'.fits' else base1fits[l]=base1[l]+'.fits.fz'
          filters = PHOTRED_GETFILTER(base1fits)
          exptime = PHOTRED_GETEXPTIME(base1fits)
          rexptime = round(exptime*10)/10.  ; rounded to nearest 0.1s
          utdate = PHOTRED_GETDATE(base1fits)
          uttime = PHOTRED_GETUTTIME(base1fits)
          dateobs = utdate+'T'+uttime
          jd = dblarr(nind1)
          for l=0,nind1-1 do jd[l]=DATE2JD(dateobs[l])
      
          ; Find matches to the reference filter in priority order
          ngdref=0 & refind=-1
          repeat begin
            refind++
            gdref = where(filters eq filtref[refind],ngdref)
          endrep until (ngdref gt 0) or (refind eq nfiltref-1)
          if ngdref gt 0 then usefiltref=filtref[refind]
          ;gdref = where(filters eq filtref,ngdref)
          ; No reference filters
          if ngdref eq 0 then begin
            printlog,logfile,'NO IMAGES IN REFERENCE FILTER - '+filtref
            printlog,logfile,'MODIFY >>photred.setup<< file parameter FILTREF'
            printlog,logfile,'FILTERS AVAILABLE: '+filters[uniq(filters,sort(filters))]
            printlog,logfile,'Failing field '+thisfield+' and going to the next'
            PUSH,failurelist,dirs[i]+'/'+base+'.als'
            goto,BOMB
          endif

          ; More than one exposure in reference filter
          ; Use image with LONGEST exptime
          if ngdref gt 1 then begin
            ; Getting image with longest exptime
            refbase = base1[gdref]
            maxind = maxloc(rexptime[gdref])
            if n_elements(maxind) gt 1 then begin  ; pick chronological first image
              si = sort(jd[gdref[maxind]])
              maxind = maxind[si[0]]  
            endif 
            ;exptime2 = exptime[gdref]
            ;maxind = first_el(maxloc(exptime2))
            refimbase = refbase[maxind[0]]
            refexptime = exptime[gdref[maxind[0]]]

            printlog,logfile,'Two images in reference filter.'
            printlog,logfile,'Picking the image with the longest exposure time'

          ; Single frame
          endif else begin
            refimbase = base1[gdref[0]]
            refexptime = exptime[gdref[0]]
          endelse

          ; Getting just the base, without the extension, e.g. "_1"
          len = strlen(refimbase)
          lenend = strlen(thisimager.separator+amps[0])
          refimbase = strmid(refimbase,0,len-lenend)

        ; FAKE, pick reference image of existing MCH file
        ;  This ensures that we use exactly the same reference frame.
        endif else begin
          ind1 = where(amp eq amps[0],nind1)
          base1 = base[ind1]
          ; Get fits/fits.fz file name
          base1fits = strarr(nind1)
          for l=0,nind1-1 do $
            if FILE_TEST(base1[l]+'.fits') eq 1 then base1fits[l]=base1[l]+'.fits' else base1fits[l]=base1[l]+'.fits.fz'          
          filters = PHOTRED_GETFILTER(base1fits)
          exptime = PHOTRED_GETEXPTIME(base1fits)
          gdref = where(file_test(base1+'.mch') eq 1,ngdref)
          if ngdref eq 0 then begin
            printlog,logfile,'/FAKE, no existing MCH file for '+amps[0]
            goto,BOMB
          endif
          if ngdref gt 1 then begin
            printlog,logfile,'/FAKE, '+strtrim(ngdref,2)+' MCH files for '+amps[0]+'. Too many!'
            goto,BOMB
          endif
          refimbase = base1[gdref[0]]
          usefiltref = filters[gdref[0]]
          refexptime = exptime[gdref[0]]

          ; Getting just the base, without the extension, e.g. "_1"
          len = strlen(refimbase)
          lenend = strlen(thisimager.separator+amps[0])
          refimbase = strmid(refimbase,0,len-lenend)
        endelse
        
        ; Reference image information
        printlog,logfile,'REFERENCE IMAGE = '+refimbase+' Filter='+usefiltref+' Exptime='+strtrim(refexptime,2)

        ;--------------
        ; AMP loop
        ;--------------
        For k=0,namps-1 do begin

          ampind = where(amp eq amps[k],nampind)
          ampfiles = base[ampind]
          printlog,logfile,''
          printlog,logfile,'AMP='+amps[k]+'  '+strtrim(nampind,2)+' FILES'

          ; Reference file
          gdref = where(ampfiles eq (refimbase+thisimager.separator+amps[k]),ngdref)

          ; We have a reference image for this amp
          ;if (ngdref gt 0 and nampind gt 1) then begin
          if (ngdref gt 0) then begin

            ; Making the input list
            inlist = ampfiles
            if nampind gt 1 then begin
              REMOVE,gdref[0],inlist
              si = sort(inlist)  ; Make sure they are in order!!!
              inlist = [ampfiles[gdref[0]],inlist[si]]
            endif
            inlist = inlist+'.als'
          
            ; Running daomatch
            ;DAOMATCH,inlist,logfile=logfile,error=daoerror,maxshift=maxshift,fake=fake,/verbose
            cmd1 = 'daomatch,["'+strjoin(inlist,'","')+'"],/verbose'
            if keyword_set(logfile) then cmd1+=',logfile="'+logfile+'"'
            if keyword_set(maxshift) then cmd1+=',maxshift='+strtrim(maxshift,2)
            if keyword_set(fake) then cmd1+=',/fake'

            mstr1 = {field:thisfield,dir:dirs[i],refbase:refimbase,chip:long(amps[k]),nfiles:long(n_elements(inlist)),cmd:cmd1,inlist:strjoin(inlist,','),mchbase:ampfiles[gdref[0]]}
            PUSH,mstr,mstr1

          ;; No reference image, or only 1 image
          endif else begin
            PUSH,failurelist,dirs[i]+'/'+ampfiles+'.als'
            ;if nampind eq 1 then printlog,logfile,'ONLY 1 FILE.  NEED AT LEAST 2'
            if ngdref eq 0 then printlog,logfile,'NO REFERENCE IMAGE FOR AMP=',amps[k]
          endelse

          BOMB1:
        Endfor ; amp loop


      ;;-------------------
      ;; Single-amp imagers
      ;;-------------------
      endif else begin

        ; Check that we have ALL files for this GROUP
        ;--------------------------------------------
        ; Some previous successes to check
        nsuccess = lists.nsuccesslines
        if (nsuccess gt 0) then begin
          printlog,logfile,''
          printlog,logfile,'Some previous successes.  Making sure we have all files for this field'

          successbase = FILE_BASENAME(lists.successlines,'.als')
          alstest = FILE_TEST(lists.successlines)
          matchind = where(stregex(successbase,'^'+thisfield+'-',/boolean) eq 1 and $
                           alstest eq 1,nmatchind)

          ; Found some matches
          if (nmatchind gt 0) then begin

            ; Check if these are already in the INLIST
            undefine,ind1,ind2,num_alreadyinlist
            MATCH,successbase[matchind],base,ind1,ind2,count=num_alreadyinlist
            num_notinlist = nmatchind - num_alreadyinlist

            ; Some not in INLIST yet
            if (num_notinlist gt 0) then begin      
              printlog,logfile,'Found '+strtrim(num_notinlist,2)+' previously successful file(s) for this group NOT YET in the '+$
                               'INLIST.  Adding.'
              indtoadd = matchind
              if num_alreadyinlist gt 0 then REMOVE,ind1,indtoadd
              PUSH,base,successbase[indtoadd]

              ; Setting REDO=1 so the files in the success list will be redone.
              if not keyword_set(redo) then begin
                printlog,logfile,'Setting REDO=1'
                redo = 1
              endif

            endif  ; some not in inlist yet
          endif  ; some files from this group in success file

          ; Make sure they are unique
          ui = uniq(base,sort(base))
          ui = ui[sort(ui)]
          base = base[ui]
          nbase = n_elements(base)

          printlog,logfile,''
        endif ; some successes 


        ; Getting the REFERENCE Image
        ;----------------------------------
        ;  use photred_pickreferenceframe.pro
        if not keyword_set(fake) then begin
          ; Get filters 
          ; Get fits/fits.fz file name
          basefits = strarr(nbase)
          for l=0,nbase-1 do $
            if FILE_TEST(base[l]+'.fits') eq 1 then basefits[l]=base[l]+'.fits' else basefits[l]=base[l]+'.fits.fz'
          filters = PHOTRED_GETFILTER(basefits)
          exptime = PHOTRED_GETEXPTIME(basefits)
          rexptime = round(exptime*10)/10.  ; rounded to nearest 0.1s
          utdate = PHOTRED_GETDATE(basefits)
          uttime = PHOTRED_GETUTTIME(basefits)
          dateobs = utdate+'T'+uttime
          jd = dblarr(nind1)
          for l=0,nind1-1 do jd[l]=DATE2JD(dateobs[l])

          ; Find matches to the reference filter in priority order
          ngdref=0 & refind=-1
          repeat begin
            refind++
            gdref = where(filters eq filtref[refind],ngdref)
          endrep until (ngdref gt 0) or (refind eq nfiltref-1)
          if ngdref gt 0 then usefiltref=filtref[refind]
          ;gdref = where(filters eq filtref,ngdref)
          ; No reference filters
          if ngdref eq 0 then begin
            printlog,logfile,'NO IMAGES IN REFERENCE FILTER - '+filtref
            printlog,logfile,'MODIFY photred.setup file parameter FILTREF'
            printlog,logfile,'FILTERS AVAILABLE: '+filters[uniq(filters,sort(filters))]
            printlog,logfile,'Failing field '+thisfield+' and going to the next'
            PUSH,failurelist,dirs[i]+'/'+base+'.als'
            goto,BOMB
          endif

          ; More than one exposure in reference filter
          ; Use image with LONGEST exptime
          if ngdref gt 1 then begin
            ; Getting image with longest exptime
            refbase = base[gdref]
            maxind = maxloc(rexptime[gdref])
            if n_elements(maxind) gt 1 then begin  ; pick chronological first image
              si = sort(jd[gdref[maxind]])
              maxind = maxind[si[0]]  
            endif 
            ;exptime2 = exptime[gdref]
            ;maxind = first_el(maxloc(exptime2))
            refimbase = refbase[maxind]
            refexptime = exptime[gdref[maxind]]

            printlog,logfile,'Two images in reference filter.'
            printlog,logfile,'Picking the image with the largest exposure time'

          ; Single frame
          endif else begin
            refimbase = base[gdref[0]]
            refexptime = exptime[gdref[0]]
          endelse

        ; FAKE, pick reference image of existing MCH file
        ;  This ensures that we use exactly the same reference frame.
        endif else begin
          ; Get fits/fits.fz file name
          basefits = strarr(nbase)
          for l=0,nbase-1 do $
            if FILE_TEST(base[l]+'.fits') eq 1 then basefits[l]=base[l]+'.fits' else basefits[l]=base[l]+'.fits.fz'
          filters = PHOTRED_GETFILTER(basefits)
          exptime = PHOTRED_GETEXPTIME(basefits)
          gdref = where(file_test(base+'.mch') eq 1,ngdref)
          if ngdref eq 0 then begin
            printlog,logfile,'/FAKE, no existing MCH file.'
            goto,BOMB
          endif
          if ngdref gt 1 then begin
            printlog,logfile,'/FAKE, '+strtrim(ngdref,2)+' MCH files. Too many!'
            goto,BOMB
          endif
          refimbase = base[gdref]
          usefiltref = filters[gdref]
          refexptime = exptime[gdref]
        endelse
        
        ; Reference image information
        printlog,logfile,'REFERENCE IMAGE = '+refimbase+' Filter='+usefiltref+' Exptime=',strtrim(refexptime,2)

        printlog,logfile,strtrim(nbase,2)+' FILES'

        ; Only 1 file
        if (nbase eq 1) then begin
          PUSH,failurelist,dirs[i]+'/'+base+'.als'
          printlog,logfile,'ONLY 1 FILE.  NEED AT LEAST 2.'
          goto,BOMB
        endif

        ; Making the input list
        inlist = base
        REMOVE,gdref[0],inlist
        si = sort(inlist)  ; Make sure they are in order!!!
        inlist = [base[gdref[0]],inlist]+'.als'

        ; Running daomatch
        ;DAOMATCH,inlist,logfile=logfile,error=daoerror,maxshift=maxshift,fake=fake,/verbose
        cmd1 = 'daomatch,["'+strjoin(inlist,'","')+'"],/verbose'
        if keyword_set(logfile) then cmd1+=',logfile="'+logfile+'"'
        if keyword_set(maxshift) then cmd1+=',maxshift='+strtrim(maxshift,2)
        if keyword_set(fake) then cmd1+=',/fake'

        mstr1 = {field:thisfield,dir:dirs[i],refbase:refimbase,chip:1L,nfiles:long(n_elements(inlist)),cmd:cmd1,inlist:strjoin(inlist,','),mchbase:base[gdref[0]]}
        PUSH,mstr,mstr1
      Endelse  ; single-amp imagers
    BOMB:
    CD,dirs[i]
  ENDFOR  ; field loop
  ; Go back to original directory
  CD,curdir
ENDFOR  ; directoryloop

;; Run PBS_DAEMON
print,strtrim(n_elements(mstr),2)+' matches to run'

cmd = "cd,'"+mstr.dir+"' & "+mstr.cmd  ; go to the directory
; Submit the jobs to the daemon
PBS_DAEMON,cmd,mstr.dir,jobs=jobs,nmulti=nmulti,prefix='match',hyperthread=hyperthread,/idle,waittime=1,scriptsdir=scriptsdir


;; Were we successful?
For i=0,n_elements(mstr)-1 do begin
  mchfile = mstr[i].dir+'/'+mstr[i].mchbase+'.mch'
  mchtest = FILE_TEST(mchfile)
  if mchtest eq 1 then mchlines=FILE_LINES(mchfile) else mchlines=0
  rawfile = mstr[i].dir+'/'+mstr[i].mchbase+'.raw'
  rawtest = FILE_TEST(rawfile)
  if rawtest eq 1 then rawlines=FILE_LINES(rawfile) else rawlines=0

  ;; Successful
  if ((mchlines eq mstr[i].nfiles) and (rawlines gt 3)) then begin
    PUSH,successlist,mstr[i].dir+'/'+strsplit(mstr[i].inlist,',',/extract)
    PUSH,outlist,mchfile
        
    ;; Getting total number of stars
    nrecords = FILE_LINES(rawfile)-3
        
    ;; Printing the results
    printlog,logfile,'NSTARS = ',strtrim(nrecords,2)
    printlog,logfile,'MCH file = ',mchfile
    printlog,logfile,'RAW file = ',rawfile
        
  ;; Failure
  endif else begin
    PUSH,failurelist,mstr[i].dir+'/'+strsplit(mstr[i].inlist,',',/extract)
    ;; failure information
    if mchtest eq 0 then printlog,logfile,mchfile+' NOT FOUND'
    if rawtest eq 0 then printlog,logfile,rawfile+' NOT FOUND'
    if rawlines le 3 then printlog,logfile,'NO SOURCES IN '+rawfile
  endelse

  printlog,logfile,''
Endfor



;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_MATCH Finished  ',systime(0)

if keyword_set(stp) then stop

end
