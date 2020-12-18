pro make_exposures_list_deep,version

;; Make list of DECAM MC exposures of DELVE+community data
;; using the NSC master decam list

if n_elements(version) eq 0 then version = 'v3'
nscdir = '/net/dl2/users/dnidever/nsc/instcal/'+version+'/lists/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;str = mrdfits(nscdir+'decam_instcal_list.fits.gz',1)
str1 = mrdfits('/net/dl2/dnidever/delve/deep/lists/decam_sexb_20200908.fits.gz',1)
str2 = mrdfits('/net/dl2/dnidever/delve/deep/lists/decam_ngc55_20200908.fits.gz',1)
str = [str1,str2]
str.prop_id = strtrim(str.prop_id,2)
nstr = n_elements(str)

;; Make sure they are public or are DELVE
; APPLY RELEASE DATE CUTS
release_date = strtrim(str.release_date,2)
release_year = long(strmid(release_date,0,4))
release_month = long(strmid(release_date,5,2))
release_day = long(strmid(release_date,8,2))
release_mjd = JULDAY(release_month,release_day,release_year)-2400000.5d0
;release_cutoff = [2019,10,17]    ; v3 - Oct 17, 2019
caldat,systime(/julian,/utc),month,day,year
release_cutoff = [year, month, day]
release_cutoff_mjd = JULDAY(release_cutoff[1],release_cutoff[2],release_cutoff[0])-2400000.5d0
;; public or DELVE or MagLiteS
gd = where(release_mjd le release_cutoff_mjd or str.prop_id eq '2019A-0305' or str.prop_id eq '2018A-0242',ngd)
str = str[gd]
nstr = ngd
print,strtrim(nstr,2),' total public DECam or DELVE/MagLiteS exposures'

jd = systime(/julian)
caldat,jd,month,day,year,hour
outfile = delvereddir+'data/decam_deep_'+string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')+'.fits'
print,'Writing to ',outfile
MWRFITS,str,outfile,/create
if file_test(outfile+'.gz') eq 1 then file_delete,outfile+'.gz'
spawn,['gzip',outfile],/noshell

;stop

end
