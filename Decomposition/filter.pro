function filter, data,freq,name,type,notch=notch,band=band
;low pass filter data
totpoints=n_elements(data.x)
sampfreq=1/(data.x(1)-data.x(0))
print,sampfreq
if not keyword_set(notch) and not keyword_set(band) then begin
	cutoff=floor(freq*(totpoints/sampfreq))
	dataspec=fft(data.y)
	print,cutoff,totpoints/2.0
	if type EQ 'l' then begin
		if cutoff LT totpoints/2. then dataspec(cutoff:totpoints-1-cutoff)=0
	endif
	if type EQ 'h' then begin
		;plot,dataspec
		print,cutoff,n_elements(dataspec),freq,sampfreq
		if cutoff GE totpoints then cutoff=totpoints-1
		dataspec(0:cutoff)=0
		dataspec(totpoints-1-cutoff:totpoints-1)=0
	endif
endif
if keyword_set(notch) then begin
	dataspec=fft(data.y)
	cutoff_low=floor(notch(0)*(totpoints/sampfreq))
	cutoff_high=floor(notch(1)*(totpoints/sampfreq))
	dataspec(cutoff_low:cutoff_high)=0.0
	dataspec(totpoints-1-cutoff_high:totpoints-1-cutoff_low)=0.0
endif

if keyword_set(band) then begin
	dataspec=fft(data.y)
	cutoff_low=floor(band(0)*(totpoints/sampfreq))
	cutoff_high=floor(band(1)*(totpoints/sampfreq))
	dataspec(0:cutoff_low)=0.0
	dataspec(cutoff_high:totpoints-1-cutoff_high)=0.0
	dataspec(totpoints-1-cutoff_low:totpoints-1)=0.0
endif


;ncolumns=totpoints-2*cutoff+1
;window=1.0;hanning(ncolumns)
;dataspec(cutoff:totpoints-cutoff-2)=window*dataspec(cutoff:totpoints-cutoff-2)
;plot,dataspec

new=fft(dataspec,1)
out=data
out.y=new
store_data,name,data={x:out.x,y:out.y}
return,out
end
