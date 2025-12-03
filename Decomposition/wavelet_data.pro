pro wavelet_data,name,dimen=dimenn,kolom=kol,trange=trange,j=j,pdens=pdens,mother_wav=mother_wav,param=param
;get_data,name,t,y
;name=tname
;if dimen2(y) gt 1 then begin
;   if n_elements(dimenn) eq 0 then message,"Need dimen!"
;   y = y[*,dimenn]
;   name = string(name)+'('+strtrim(dimenn,2)+')'
;endif
help,scale
t=name.x
y=name.y
if n_elements(trange) eq 2 then begin
   w = where((t le trange[1]) and (t ge trange[0]),nw)
   if nw eq 0 then begin
      message,/info,'No data in time range'
      return
   endif
   t=t[w]
   y=y[w]
endif
nt = n_elements(t)
dt = (t[nt-1]-t[0])/(nt-1)
interp_gap,t,y,/verbose
n = n_elements(y)
if keyword_set(mother_wav) then mother=mother_wav else mother='MORLET'
wave = wavelet(y,dt,pad=2,period=period,coi=coi,signif=signif,scale=scale,/verb,j=j,mother=mother,param=param)
;print,signif
if keyword_set(pdens) then pwave=abs(wave*conj(wave)) else pwave=wave
;p=average(pwave,1)
freq=1/period
delta_freq=freq-shift(freq,-1)
delta_freq(n_elements(freq)-1)=delta_freq(n_elements(freq)-2)
;print,'freq',freq
;print,'delta_freq',delta_freq
pkg = .0003 * freq^(-5./3.)
if keyword_set(kol) then pwave = pwave/(replicate(1.,n) # pkg)
;if keyword_set(pdens) then pwave = pwave/(replicate(1.,n) # (delta_freq))/1000.
;if keyword_set(pdens) then pwave = pwave/(replicate(1.,n) # (1./(freq*1000.)))
;if keyword_set(pdens) then pwave = pwave*1000.

;if keyword_set(pdens) then pwave = pwave/(sqrt(2.)/(dt*!DPI));(replicate(1.,n) # sqrt(period/(2.*!DPI*scale))*1.0/dt*2.0/sqrt(!DPI))
;corrected normalisation based on the paper by Torrence and Compo 'A practical guide to wavelet analysis' Equation 6.
;The normalisation is approximately 1/sqrt(2.*!DPI)*1.0/dt*2.0/sqrt(!DPI) where dt is the sampling period of the data
;and the 2.0/sqrt(!DPI) is a correction close to 1 taken from 'Analysis methods for multi-sp. data' Eq 1.43  

;if keyword_set(pdens) then pwave = pwave/dt*sqrt(!DPI)/2.0;this seems to be correct normalisation the sqrt(!DPI)/2.0 is close to 1.0 and comes from 'Analysis methods for multi-sp. data' Eq 1.43 and is intended for Morlet wavelets
if keyword_set(pdens) then pwave = pwave*dt*2.0*sqrt(!DPI)/2.0

;print,'period',period
;print,'scale',scale
;print,'dt',dt
;specplot,t,freq,pwave,lim={ylog:1},/no_inter
if keyword_set(pdens) then store_data,'pow'+'_wvlt',data={x:t,y:pwave,v:period},dlim={spec:1,ylog:1,zlog:1,ystyle:1,no_interp:1}$
else store_data,'wvlt',data={x:t,y:pwave,v:period},dlim={spec:1,ylog:1,zlog:1,ystyle:1,no_interp:1} 

;store_data,name+'_coi' , data={x:t,y:coi}
;store_data,name+'_wvs',data=name+'_wvlt '+name+'_coi',dlim={yrange:minmax(period),ystyle:1}
end

