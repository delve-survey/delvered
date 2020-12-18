pro make_exposures_list,version

;; Make list of DECAM MC exposures of DELVE+community data
;; using the NSC master decam list

if n_elements(version) eq 0 then version = 'v3'
nscdir = '/net/dl2/users/dnidever/nsc/instcal/'+version+'/lists/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;str = mrdfits(nscdir+'decam_instcal_list.fits.gz',1)
;str = mrdfits('/net/dl2/dnidever/delve/lists/delvemc_info_20121022_20200216.fits.gz',1)
;str = mrdfits('/net/dl2/dnidever/delve/lists/delvemc_info_20121022_20200301.fits.gz',1)
str = mrdfits('/net/dl2/dnidever/delve/lists/delvemc_info_20121022_20200318.fits.gz',1)
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

;; Cutting on coordinates
glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
cel2lmc,str.ra,str.dec,lmcpa,lmcrad
cel2smc,str.ra,str.dec,smcpa,smcrad
gd = where(abs(glat) gt 10 and (lmcrad lt 25 or smcrad lt 15),ngd)
mc = str[gd]
print,strtrim(ngd,2),' MC exposures'


jd = systime(/julian)
caldat,jd,month,day,year,hour
outfile = delvereddir+'data/decam_mcs_'+string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')+'.fits'
print,'Writing to ',outfile
MWRFITS,mc,outfile,/create
if file_test(outfile+'.gz') eq 1 then file_delete,outfile+'.gz'
spawn,['gzip',outfile],/noshell

;stop

end
