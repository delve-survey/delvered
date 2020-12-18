pro convert_mwmodels,files

;; Convert MW models from ASCII to FITS files

radeg = 180.0d0 / !dpi

delvereddir = '/home/dnidever/projects/delvered/'
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1,/silent)
brickstr.brickname = strtrim(brickstr.brickname,2)

nfiles = n_elements(files)
print,strtrim(nfiles,2),' file(s) to process'

for i=0,nfiles-1 do begin
  brick = (strsplit(file_basename(files[i]),'_',/extract))[1]
  print,strtrim(i+1,2),' ',files[i]
  info = file_info(files[i])
  if info.exists eq 0 or info.size eq 0 then begin
    print,files[i],' is empty'
    continue
  endif
  arr = importascii(files[i],/header,/silent)
  narr = n_elements(arr)
  print,'  ',strtrim(narr,2),' rows'
  ;; New schema
  schema = {brick:'',ra:0.0d0,dec:0.0d0,ring256:0L,u:0.0,g:0.0,r:0.0,i:0.0,z:0.0,feh:0.0,age:0.0,popid:0.0}
  new = replicate(schema,narr)
  struct_assign,arr,new,/nozero
  new.feh = arr._fe_h_
  new.brick = brick
  MATCH,brickstr.brickname,brick,ind1,ind2,/sort,count=nmatch
  if nmatch eq 0 then begin
    print,'No match for this brick'
    continue
  endif
  new.ra = brickstr[ind1[0]].ra
  new.dec = brickstr[ind1[0]].dec
  ANG2PIX_RING,256L,(90-brickstr[ind1[0]].dec)/radeg,brickstr[ind1[0]].ra/radeg,hpix
  new.ring256 = hpix[0]
  outfile = files[i]+'.fits'
  print,'  Writing to ',outfile
  MWRFITS,new,outfile,/create
  SPAWN,['gzip','-f',outfile],/noshell
endfor

;stop

end
