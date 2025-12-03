function ave_power_spec_dynamic, data, $
                NPTS = npts, $
               	ave_spec, $
               	ave_phase, $
               	fft_keep, $
                OVERLAP = overlap
                

if not keyword_set(npts) then npts=1024
fft_out = fltarr(npts)

W = npts*total((hanning(npts))^2)


sample=(data.x(1)-data.x(0))
    

; If overlap keyword is set, slide half of a segment. If not, slide a
; complete segment.


if not keyword_set(overlap) then overlap=0.5

n_loop=n_elements(data.x)/(npts-npts*overlap)-2.

seg_step=npts*overlap
;if keyword_set(overlap) then begin
;    n_loop=2*n_ave-1
;    seg_step=npts/2
; else begin
;    seg_step=npts
;    n_loop=n_ave
;endelse
ave_spec = fltarr(npts/2+1)
ave_phase=fltarr(npts/2+1)
dyn_spec = fltarr(n_loop,npts/2+1)
fft_keep=make_array(n_loop,npts,/dcomplex)
time=make_array(n_loop,/double)
; Initialize array indexing
start=long(0)
stop=long(npts-1)
print,n_loop

; Define frequency array and delta frequency between each point.
freq = findgen(npts/2+1)/(npts*sample)
delta_freq = freq(1)-freq(0)
; Begin loop for averaging segments through the time series.

for i = long(0),long((n_loop-1)) do begin
    ; Window the data using a Hanning envelope.
    ;wdata(start:stop) = data.y(start:stop)*hanning(npts)
    ; Take FFT of the windowed segment.
    time(i)=data.x(start)+(data.x(stop)-data.x(start))/2.0
    fft_out = fft((data.y(start:stop)*hanning(npts)),+1)
    dyn_spec(i,0) = (1./W)*fft_out(0)*conj(fft_out(0))
    dyn_spec(i,1:(npts/2-1)) = (1./W)*(fft_out(1:(npts/2-1))*conj(fft_out(1:(npts/2-1))) + $
      reverse(fft_out(npts/2+1:npts-1)) * $
      conj(reverse(fft_out(npts/2+1:npts-1))))
    dyn_spec(i,npts/2) = (1./W)*(fft_out(npts/2)*conj(fft_out(npts/2)))
    fft_keep(i,*)=fft_out
    ;plot,freq(1:(npts/2-1)),dyn_spec(i,1:(npts/2-1)),/xlog,/ylog
    ;print,time_to_str(time(i),/ms)
    stop=stop+seg_step
    start=start+seg_step
endfor



; Convert spectral power (amplitude^2) to (amplitude^2/Hz)
dyn_spec(*,1:npts/2-1) = dyn_spec(*,1:npts/2-1)/delta_freq
; The zeroeth and NPTS/2th PSD bins are only half as wide as the rest
; of the frequency bins -- add factor of 2 for (amplitude^2/Hz)
dyn_spec(*,0) = 2*dyn_spec(*,0)/delta_freq
dyn_spec(*,npts/2) = 2*dyn_spec(*,npts/2)/delta_freq

for j=0L,n_elements(freq)-1 do begin
	ave_phase(j)=total(atan(fft_keep(*,j),/phase)/!dtor)/n_loop
	ave_spec(j)=total(dyn_spec(*,j))/n_loop
endfor

ave_spec={x:freq,y:ave_spec}
ave_phase={x:freq,y:ave_phase}

data_spec={x:time,y:dyn_spec,v:freq,spec:1}
store_data,'dy_sp',data=data_spec
options,'dy_sp','ylog',1
options,'dy_sp','zlog',1
options,'dy_sp','yrange',[freq(1),freq(n_elements(freq)-1)]
options,'dy_sp','ystyle',1
options,'dy_sp','y_no_interp',1
options,'dy_sp','x_no_interp',1
options,'dy_sp','z_no_interp',1
return,ave_spec;data_spec
end
