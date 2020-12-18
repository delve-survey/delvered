pro dlvexpsumplot

;; DELVE exposures processing summary figure

delvereddir = '/home/dnidever/projects/delvered/'

;expfile = delvereddir+'data/decam_mcs_20191017.fits.gz'
;allexpstr = MRDFITS(expfile,1,/silent)
;nallexp = n_elements(allexpstr)
;allexpstr.fluxfile = strtrim(allexpstr.fluxfile,2)
;allexpstr.maskfile = strtrim(allexpstr.maskfile,2)
;allexpstr.wtfile = strtrim(allexpstr.wtfile,2)
;allexpstr.expnum = strtrim(allexpstr.expnum,2)
;;; Fix names for hulk/thing
;;;  /net/mss1 -> /mss1
;allexpstr.fluxfile = strmid(allexpstr.fluxfile,4)
;allexpstr.maskfile = strmid(allexpstr.maskfile,4)
;allexpstr.wtfile = strmid(allexpstr.wtfile,4)
;filt = strmid(allexpstr.filter,0,1)
;;gdexp = where((filt eq 'g' or filt eq 'r' or filt eq 'i') and
;;allexpstr.exposure ge 90,ngdexp)
;gdexp = where((filt eq 'u' or filt eq 'g' or filt eq 'r' or filt eq 'i' or filt eq 'z' or filt eq 'Y') and allexpstr.exposure ge 90,ngdexp)
;print,strtrim(ngdexp,2),' EXPOSURES pass the cuts'
;expstr = allexpstr[gdexp]
;add_tag,expstr,'mlon',0.0,expstr
;add_tag,expstr,'mlat',0.0,expstr
;glactc,expstr.ra,expstr.dec,2000.0,glon,glat,1,/deg
;gal2mag,glon,glat,mlon,mlat
;expstr.mlon = mlon
;expstr.mlat = mlat
;
;23277 exposures

delvedir = '/net/dl1/users/dnidever/delve/'
dirs = file_search(delvedir+'exposures/201?????',/test_directory,count=ndirs)
dirs = file_basename(dirs)

;; Wrap up the nightly exposures files
print,'Wrapping up '+strtrim(ndirs,2)+' exposure files'
expfiles = delvedir+'exposures/'+dirs+'/'+dirs+'_exposures.fits'
for i=0,ndirs-1 do begin
  if i mod 50 eq 0 then print,i
  if file_test(expfiles[i]) eq 1 then PUSH,expstr,mrdfits(expfiles[i],1,/silent)
endfor
add_tag,expstr,'mlon',0.0,expstr
add_tag,expstr,'mlat',0.0,expstr
glactc,expstr.ra,expstr.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
expstr.mlon = mlon
expstr.mlat = mlat
nexpstr = n_elements(expstr)
print,strtrim(nexpstr,2),' total exposures'

;; Wrap up all the nightly summary files
nightsumfiles = delvedir+'exposures/'+dirs+'/'+dirs+'_summary.fits'
test = file_test(nightsumfiles)
g = where(test eq 1,ng)
print,'Wrapping up '+strtrim(ng,2)+' summary files'
nights = strarr(ng)
nexp = lonarr(ng)
undefine,sumstr
for i=0,ng-1 do begin
  print,nightsumfiles[g[i]]
  inight = file_basename(nightsumfiles[g[i]],'_summary.fits')
  nights[i] = inight
  str1 = mrdfits(nightsumfiles[g[i]],1)
  nstr1 = n_elements(str1)
  if nstr1 gt 0 and size(str1,/type) eq 8 then begin
    nexp[i] = nstr1
    add_tag,str1,'night','',str1
    str1.night = inight
    push,sumstr,str1
  endif
endfor
nsumstr = n_elements(sumstr)
add_tag,sumstr,'mlon',0.0,sumstr
add_tag,sumstr,'mlat',0.0,sumstr
glactc,sumstr.ra,sumstr.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
sumstr.mlon = mlon
sumstr.mlat = mlat
print,strtrim(nsumstr,2),' finished exposures'

;; Make the figure
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
time = string(year,month,day,hour,minute,second,format='(I4,I02,I02,I02,I02,I02)')
setdisp
!p.font = 0
psfile = delvedir+'exposures/summary/plots/dlvexpsumplot_'+time
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=12.0,ysize=10.5
plotc,expstr.mlon,expstr.mlat,/nodata,ps=3,/xflip,xtit='L!dMS!n',ytit='B!dMS!n',charsize=1.8,tit='DELVERED Processing',position=[0.09,0.09,0.97,0.95]
filters = ['u','g','r','i','z','Y']
nfilters = n_elements(filters)
colors = [100,130,160,180,200,250]
filt = strmid(expstr.filter,0,1)
ntotexp = lonarr(nfilters)
for i=0,nfilters-1 do begin
  ind = where(filt eq filters[i],nind)
  ntotexp[i] = nind
  if nind gt 0 then oplot,[expstr[ind].mlon],[expstr[ind].mlat],ps=1,sym=0.2,co=colors[i],thick=2
endfor
nfinexp = lonarr(nfilters)
for i=0,nfilters-1 do begin
  ind = where(sumstr.filter eq filters[i],nind)
  nfinexp[i] = nind
  if nind gt 0 then oplot,[sumstr[ind].mlon],[sumstr[ind].mlat],ps=6,sym=0.5,co=colors[i],thick=2
endfor
legend_old,['All ('+strtrim(nexpstr,2)+')','Done ('+strtrim(nsumstr,2)+')'],textcolor=[0,250],psym=[1,6],colors=[0,250],symsize=[0.4,0.7],/top,/left,box=0,charsize=1.9
legend_old,filters+' ('+strtrim(nfinexp,2)+'/'+strtrim(ntotexp,2)+')',textcolor=colors,/top,/left,box=0,charsize=2.0,pos=[-22,28]
ps_close
ps2png,psfile+'.eps',/eps
spawn,['epstopdf',psfile+'.eps'],/noshell
print,'Summary figure saved to ',psfile+'.eps/png/pdf'

;stop

end
