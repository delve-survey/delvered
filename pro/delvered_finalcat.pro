pro delvered_finalcat

;; Create the FINAL DELVE-MC catalog using the brick catalogs
;;  making sure to only use the objects from the brick's unique area

t0 = systime(1)
radeg = 180.0d0 / !dpi

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
finaldbfile = brickdir+'db/delvered_final.db'
tempdir = '/tmp/'

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1)

;; Brick directories
bdirs = file_search(brickdir+'????/????????',/test_directory,count=nbdirs)
bricks = file_basename(bdirs)
print,strtrim(nbdirs,2),' brick directories'
;; Match them to the list
MATCH,brickstr.brickname,bricks,ind1,ind2,/sort,count=nmatch
if nmatch eq 0 then begin
  print,'No brick directories that match brick list'
  return
endif
bricks = bricks[ind2]
bdirs = bdirs[ind2]
nbricks = n_elements(bricks)
bstr = brickstr[ind1]
if nmatch lt nbdirs then print,strtrim(nbricks,2),' bricks retained'

;; Final schema
fschema = {brick:'',id:'',ring256:0L,ra:0.0d0,dec:0.0d0,ndet:0L,gmag:99.99,gerr:9.99,gscatter:99.99,ndetg:0L,$
           rmag:99.99,rerr:9.99,rscatter:99.99,ndetr:0L,imag:99.99,ierr:9.99,iscatter:99.99,ndeti:0L,$
           zmag:99.99,zerr:9.99,zscatter:99.99,ndetz:0L,chi:0.0,sharp:0.0,prob:0.0,sfd_ebv:0.0,$
           mag_auto:999.0,magerr_auto:999.0,asemi:999.0,bsemi:999.0,theta:999.0,ellipticity:999.0,fwhm:999.0}
ftags = tag_names(fschema)

;; Brick loop
for i=0,nbdirs-1 do begin
  idir = bdirs[i]+'/'
  brick = bricks[i]
  print,strtrim(i+1,2),' ',brick,'   ',stringize(bstr[i].ra,ndec=4),' ',stringize(bstr[i].dec,ndec=4)
  catfile = idir+brick+'.fits.gz'
  metafile = idir+brick+'_meta.fits'
  test = file_test([catfile,metafile])
  if min(test) eq 0 then begin
    if test[0] eq 0 then print,catfile,' NOT FOUND'
    if test[1] eq 0 then print,metafile,' NOT FOUND'
    goto,BOMB
  endif
  cat = MRDFITS(catfile,1,/silent)
  ncat = n_elements(cat)
  ctags = tag_names(cat)
  meta = MRDFITS(metafile,1,/silent)
  nmeta = n_elements(meta)
  print,'  ',strtrim(ncat,2),' objects'

  ginside = where(cat.brickuniq eq 1,ninside)  
  print,'  Keeping '+strtrim(ninside,2)+' objects inside the unique brick area'
  if ninside eq 0 then begin
    print,'No objects left to save'
    return
  endif
  cat = cat[ginside]
  ncat = n_elements(cat)

  ;; Convert to final format
  new = replicate(fschema,ncat)
  struct_assign,cat,new,/nozero
  new.brick = brick
  new.id = brick+'.'+strtrim(cat.id,2)
  ANG2PIX_RING,256L,(90-cat.dec)/radeg,cat.ra/radeg,hpix
  new.ring256 = hpix
  new.sfd_ebv = cat.ebv

  ;; Concatenate

  ;; Write to a temporary file
  tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
  tmpfile += '.fits'
  MWRFITS,new,tmpfile,/create

  ;; maybe load the database every ~5th brick or so
  print,'  Writing to the database'
  pylines = 'python -c "from delvered import delvered_db as delvedb;'+$
            'from astropy.io import fits;'+$
            "cat = fits.getdata('"+tmpfile+"',1);"+$
            "delvedb.writecat2db(cat,'"+finaldbfile+"',table='object')"+'"'
  SPAWN,pylines,out,errout
  file_delete,tmpfile,/allow

  ;; KEEP TRACK OF INDIVIDUAL MEASUREMENTS? 

  BOMB:
endfor

;; Index database
print,'Indexing the table'
pylines = 'python -c "from delvered import delvered_db as delvedb;'+$
          "delvedb.createindexdb('"+finaldbfile+"',col='ra',table='object',unique=False);"+$
          "delvedb.createindexdb('"+finaldbfile+"',col='dec',table='object',unique=False)"+'"'
SPAWN,pylines,out,errout

;stop

end
