pro clu_pow_spec_fixed,tseries,av_spec,nptsfft,slide,$
store_struc=store_struc,changesf=changesf,wavelet=wavelet,median_val=median_val,mother_wav=mother_wav,name=name,param=param

;get uniform time steps and put e and b on same time scale

;n_points=n_elements(tseries.x)
;start_time=tseries.x(0)
;end_time=tseries.x(n_points-1)
;if keyword_set(changesf) then n_points=floor(n_points/changesf)
;deltat=(end_time-start_time)/n_points
;time_scale=start_time+findgen(n_points)*deltat

if not keyword_set(name) then name='pow_wvlt'


keep=where(finite(tseries.y) EQ 1)
out=interpol(tseries.y(keep),tseries.x(keep),tseries.x)
data_fix={x:tseries.x,y:out}
;tseries=data_fix
if not keyword_set(wavelet) and not keyword_set(mother_wav) then begin
	
	;sp=spec1(tseries,nptsfft,slide,averages=1,store_struc=1)
	out=ave_power_spec_dynamic(data_fix,npts=nptsfft,overlap=slide)
	av_spec=out
	if keyword_set(name) then begin
		get_data,'dy_sp',data=sp
		store_data,name,data=sp
	endif else begin
		name='dy_sp'
	endelse

endif else begin
	wavelet_data,data_fix,/pdens,mother_wav=mother_wav
	get_data,'pow_wvlt',data=sp
	sp.v=reverse(1./sp.v)
	npts=n_elements(sp.y(*,0))
	for j=0L,npts-1 do sp.y(j,*)=reverse(sp.y(j,*),2)
	store_data,name,data=sp
		av_spec_array=make_array(n_elements(sp.v),/double)
		for j=0, n_elements(sp.v)-1 do begin
			if not keyword_set(median_val) then temp=total(sp.y(*,j))/n_elements(sp.x) else temp=median(sp.y(*,j))
			if finite(temp) then av_spec_array(j)=temp else av_spec_array(j)=!values.f_nan
			endfor
	av_spec={x:sp.v,y:av_spec_array}

endelse

options,name,'spec',1
options,name,'ylog',1
options,name,'zlog',1
options,name,'yrange',[av_spec.x(1),av_spec.x(n_elements(av_spec.x)-1)]

return
end
