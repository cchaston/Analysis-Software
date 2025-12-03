pro mhd_mode_decomp_te_minimization_angle,mode_in,bx,by,bz,$
		    vx,vy,vz,$
		    dens,te,ti,$
		    waveangle,frequency,wave_vec_sp,$
	  	    spec_av=spec_av,$
		    theta_minvar=theta_minvar,$
		    start_vector=start_vector
		    
		    
		

;this pro takes magnetic field and velocity measurements, determines the direction of k from min variance analysis and then uses the technique described
;in the ISSI book to separate out the mhd modes. The background magnetic field and density is taken to be frequency dependent consequently the magnitude
;of va and beta are functions of frequency. Default is to use wavelet transform however the same analysis can be done using a single fourier transform
;of the whole interval by setting the keyword fourier. Magnetic field, velocity and density data should be all on of the same cadence. 

;Unlike mhd_mode_decomp mhd_mode_decomp_te,pro Te is a variable with the same sampling as the b-field data 


;the direction of k is derived from the minimum variance direction. The Pi ambiguity in the direction for k given by this technique is not important
;and will just change the sense of propagation of the wave - i.e forward or backward propagating (see below) but not the mode. 
;this approach relies on elliptical polarisation of the magnetic field. The mhd model for mode recognition makes no assumptions with regard to circular or linear polarisation and so this will not effect the mode id.

;can run tests for synthetic data with the prpoerties of the slow fast and shear alfven wave by seeting the corresponding keywords.
;This test can be done for both wavelets (default) or fourier approaches using the fourier keyword. 
;In this case the min var analysis is skipped and a prescribed k direction defined in the code for each mode is used.

;can also set the k direction and again avoid min var analysis for use with data by setting the keyword wavector=[kx,ky,kz] in field-aligned coordinates

;INPUT- any orthogonal 3-D coordinate system
;input b is in nT
;input v is in km/s
;input dens is in cm-3
;input te is in eV

;OUTPUT
;output is in the variable called modes containing the fractional power in each mode with format: (number of times, number of frequencies,mode). Mode is from 0-5 with: 0,1 forward and backward propagating shear Alfven waves ;2,3 forward and backward propagating fast mode Alfven waves; 4,5 forward and backward propagating slow mode Alfven waves.  

;need to compile these first

;.run /home/ccc/idl_prog/arctangent.pro
;.run /home/ccc/idl_prog/fieldrot_fast.pro
;.run /home/ccc/idl_prog/interp_gap.pro
;.run /home/ccc/idl_prog/wavelet.pro
;.run /home/ccc/idl_prog/wavelet_data.pro
;.run /home/ccc/idl_prog/filter.pro
;.run /home/ccc/idl_prog/ave_power_spec_dynamic.pro
;.run /home/ccc/idl_prog/clu_pow_spec_fixed.pro
common fac,bx_fac_spec,by_fac_spec,bz_fac_spec,vx_fac_spec,vy_fac_spec,vz_fac_spec,rho_zero,dens_spec
common wave_mode,mode
common sv,state_vec
common wavenumber,k,kangle
common plasma_params,Bo,T_i,T_e,n_i
common minvar,wn_x,wn_y,wn_z
common coverage,mag_spec_frac
common wave_vec_sp,wave_vec
common in_vec_sp,in_vec
mi=1.67e-27;ion mass in kgs
mo=4.0*!DPI*1.0e-7

;*******************************************create a temperature array if not input*****************************************************
help,te
if (data_type(Te) EQ 2 or data_type(Te) EQ 4) then  Te_r=make_array(n_elements(bx.x),value=Te) else te_r=te.y  
if (data_type(Ti) EQ 2 or data_type(Ti) EQ 4) then  Ti_r=make_array(n_elements(bx.x),value=Ti) else ti_r=ti.y  
;***************************************NOW DO THE DECOMPOSITION*************************************************************************

;get the power spec in each quantity
clu_pow_spec_fixed,bx,name='bx_spec',/wavelet
clu_pow_spec_fixed,by,name='by_spec',/wavelet
clu_pow_spec_fixed,bz,name='bz_spec',/wavelet
clu_pow_spec_fixed,vx,name='vx_spec',/wavelet
clu_pow_spec_fixed,vy,name='vy_spec',/wavelet
clu_pow_spec_fixed,vz,name='vz_spec',/wavelet
clu_pow_spec_fixed,dens,name='dens_spec',/wavelet

;get the background magnetic field for each frequency


totpoints = n_elements(bx.x)
sampfreq=1.0/(bx.x(1)-bx.x(0))
mode=mode_in



;take the wavelet transform of each quantity

	;magnetic field quantities
	
	wavelet_data,bx,mother_wav=mother_wav
	get_data,'wvlt',data=bx_wvlt
	bx_wvlt.v=1.0/bx_wvlt.v
	frequency=bx_wvlt.v
	wavelet_data,by,mother_wav=mother_wav
	get_data,'wvlt',data=by_wvlt
	by_wvlt.v=1.0/by_wvlt.v
	wavelet_data,bz,mother_wav=mother_wav
	get_data,'wvlt',data=bz_wvlt
	bz_wvlt.v=1.0/bz_wvlt.v

	;velocity quantitites

	wavelet_data,vx,mother_wav=mother_wav
	get_data,'wvlt',data=vx_wvlt
	vx_wvlt.v=1.0/vx_wvlt.v
	wavelet_data,vy,mother_wav=mother_wav
	get_data,'wvlt',data=vy_wvlt
	vy_wvlt.v=1.0/vy_wvlt.v
	wavelet_data,vz,mother_wav=mother_wav
	get_data,'wvlt',data=vz_wvlt
	vz_wvlt.v=1.0/vz_wvlt.v

	;density

	wavelet_data,dens,mother_wav=mother_wav
	get_data,'wvlt',data=dens_wvlt
	dens_wvlt.v=1.0/dens_wvlt.v


	n_freq=n_elements(bx_wvlt.v)


	;get the  background magnetic field for each frequency

	bo_x=make_array(totpoints,n_freq,/double)
	bo_y=make_array(totpoints,n_freq,/double)
	bo_z=make_array(totpoints,n_freq,/double)
	denso=make_array(totpoints,n_freq,/double)
	rho_o=make_array(totpoints,n_freq,/double)
	ratio=10.;this sets the amount by which the frequency of the reference field is different from the frequency of the selected frequency bin


	for kk=0,n_freq-1 do begin
		
		freq_filt=bx_wvlt.v(kk)/ratio
		
		if freq_filt LT 1.0/((totpoints-1)/sampfreq) then freq_filt=1.0/((totpoints-1)/sampfreq)
		if n_elements(bx.y) LT 200000. then begin
		bx_f_l=filter(bx,freq_filt,'bx_f_l','l')
		by_f_l=filter(by,freq_filt,'by_f_l','l')
		bz_f_l=filter(bz,freq_filt,'bz_f_l','l')
		;print,'floor(sampfreq/freq_filt)',floor(sampfreq/freq_filt),sampfreq,freq_filt
		denso_f_l={x:dens.x,y:smooth(dens.y,floor(sampfreq/freq_filt),/edge_mirror)};filter(dens,freq_filt,'dens_f_l','l')	
			endif else begin
			bx_f_l=time_domain_filter(bx,0.,freq_filt)
			by_f_l=time_domain_filter(by,0.,freq_filt)
			bz_f_l=time_domain_filter(bz,0.,freq_filt)
			denso_f_l={x:dens.x,y:smooth(dens.y,floor(sampfreq/freq_filt),/edge_mirror)};time_domain_filter(dens,0.,freq_filt)
		endelse
		print,'filt_freq_bin,freq',kk,freq_filt		
		bo_x(*,kk)=bx_f_l.y
		bo_y(*,kk)=by_f_l.y
		bo_z(*,kk)=bz_f_l.y

		denso(*,kk)=denso_f_l.y
		rho_o(*,kk)=mi*denso(*,kk)*1.0e6;in kg/m^3
		;plot,denso(*,kk)
		;read,'tt',test
	endfor

	;now do the rotation to get spectra in FAC coords
	
	bx_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)
	by_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)
	bz_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)

	vx_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)
	vy_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)
	vz_fac_wvlt=make_array(totpoints,n_freq,/dcomplex)
	start_vector_array=make_array(totpoints,n_freq,3,/double)
	if keyword_set(start_vector) then start_vector_old_coords=[start_vector(0),start_vector(1),start_vector(2)] else start_vector_old_coords=[1,1,1]
	if 10 NE 10 then begin;keyword_set(alf_test) or keyword_set(fast_test) or keyword_set(slow_test) then begin
		;can use the direct allocations here by removing the 10 NE 10 part and using the rest of the line
		;for the mode test quantitites but will still work if leave as is for the full case
		
		bx_fac_wvlt=bx_wvlt.y*1.0e-9;in T
		by_fac_wvlt=by_wvlt.y*1.0e-9
		bz_fac_wvlt=bz_wvlt.y*1.0e-9

		vx_fac_wvlt=vx_wvlt.y*1000.;in m/s
		vy_fac_wvlt=vy_wvlt.y*1000.
		vz_fac_wvlt=vz_wvlt.y*1000.


	endif else begin

		b_rot_d=make_array(n_freq,3,/double)
		b_rot_i=make_array(n_freq,3,/double)

		v_rot_d=make_array(n_freq,3,/double)
		v_rot_i=make_array(n_freq,3,/double)

		for j=0L,totpoints-1 do begin  
			for kk=0L,n_freq-1 do begin
								
				;print,'bo',bo_x(j,kk),bo_y(j,kk),bo_z(j,kk)
				;print,'freq',bx_wvlt.v(kk)
				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),double(bx_wvlt.y(j,kk)),double(by_wvlt.y(j,kk)),double(bz_wvlt.y(j,kk)))
				b_rot_d(kk,0)=temp(0)
				b_rot_d(kk,1)=temp(1)
				b_rot_d(kk,2)=temp(2)

				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),imaginary(bx_wvlt.y(j,kk)),imaginary(by_wvlt.y(j,kk)),imaginary(bz_wvlt.y(j,kk)))
				b_rot_i(kk,0)=temp(0)
				b_rot_i(kk,1)=temp(1)
				b_rot_i(kk,2)=temp(2)

				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),double(vx_wvlt.y(j,kk)),double(vy_wvlt.y(j,kk)),double(vz_wvlt.y(j,kk)))
				v_rot_d(kk,0)=temp(0)
				v_rot_d(kk,1)=temp(1)
				v_rot_d(kk,2)=temp(2)

				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),imaginary(vx_wvlt.y(j,kk)),imaginary(vy_wvlt.y(j,kk)),imaginary(vz_wvlt.y(j,kk)))
				v_rot_i(kk,0)=temp(0)
				v_rot_i(kk,1)=temp(1)
				v_rot_i(kk,2)=temp(2)
				
				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),start_vector_old_coords(0),start_vector_old_coords(1),start_vector_old_coords(2))
				start_vector_array(j,kk,0)=temp(0)
				start_vector_array(j,kk,1)=temp(1)
				start_vector_array(j,kk,2)=temp(2)
				
			endfor

			;reassemble the rotated spectral data - 

			bx_fac_wvlt(j,*)=dcomplex(b_rot_d(*,0),b_rot_i(*,0))*1.0e-9;in T
			by_fac_wvlt(j,*)=dcomplex(b_rot_d(*,1),b_rot_i(*,1))*1.0e-9
			bz_fac_wvlt(j,*)=dcomplex(b_rot_d(*,2),b_rot_i(*,2))*1.0e-9

			vx_fac_wvlt(j,*)=dcomplex(v_rot_d(*,0),v_rot_i(*,0))*1000.;in m/s
			vy_fac_wvlt(j,*)=dcomplex(v_rot_d(*,1),v_rot_i(*,1))*1000.
			vz_fac_wvlt(j,*)=dcomplex(v_rot_d(*,2),v_rot_i(*,2))*1000.

		endfor

	endelse

;store spectra in FAC coords

store_data,'Bx_fac_wvlt',data={x:bx_wvlt.x,y:bx_fac_wvlt,v:bx_wvlt.v,spec:1}
store_data,'By_fac_wvlt',data={x:bx_wvlt.x,y:by_fac_wvlt,v:bx_wvlt.v,spec:1}
store_data,'Bz_fac_wvlt',data={x:bx_wvlt.x,y:bz_fac_wvlt,v:bx_wvlt.v,spec:1}

store_data,'Vx_fac_wvlt',data={x:bx_wvlt.x,y:vx_fac_wvlt,v:bx_wvlt.v,spec:1}
store_data,'Vy_fac_wvlt',data={x:bx_wvlt.x,y:vy_fac_wvlt,v:bx_wvlt.v,spec:1}
store_data,'Vz_fac_wvlt',data={x:bx_wvlt.x,y:vz_fac_wvlt,v:bx_wvlt.v,spec:1}
time=bx_wvlt.x

if keyword_set(spec_av) then begin

red_points=floor(totpoints/spec_av)
tmpbx=make_array(red_points,n_freq,/dcomplex)
tmpby=make_array(red_points,n_freq,/dcomplex)
tmpbz=make_array(red_points,n_freq,/dcomplex)
tmpvx=make_array(red_points,n_freq,/dcomplex)
tmpvy=make_array(red_points,n_freq,/dcomplex)
tmpvz=make_array(red_points,n_freq,/dcomplex)
tmpd=make_array(red_points,n_freq,/dcomplex)
tmpdo=make_array(red_points,n_freq,/double)
tmpbox=make_array(red_points,n_freq,/double)
tmpboy=make_array(red_points,n_freq,/double)
tmpboz=make_array(red_points,n_freq,/double)
time=make_array(red_points,/double)
Te_tmp=make_array(red_points,/double)
Ti_tmp=make_array(red_points,/double)
count=0
for jjj=0,totpoints-spec_av,spec_av do begin
	time(count)=total(bx_wvlt.x(jjj:jjj+spec_av-1))/spec_av
	Te_tmp(count)=total(Te_r(jjj:jjj+spec_av-1))/spec_av
	Ti_tmp(count)=total(Te_r(jjj:jjj+spec_av-1))/spec_av
	for kkk=0,n_freq-1 do begin
		
    		tmpbx(count,kkk)=total(bx_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpby(count,kkk)=total(by_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpbz(count,kkk)=total(bz_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av	

		tmpvx(count,kkk)=total(vx_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpvy(count,kkk)=total(vy_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpvz(count,kkk)=total(vz_fac_wvlt(jjj:jjj+spec_av-1,kkk))/spec_av
		
		tmpd(count,kkk)=total(dens_wvlt.y(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpdo(count,kkk)=total(denso(jjj:jjj+spec_av-1,kkk))/spec_av
		
		tmpbox(count,kkk)=total(bo_x(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpboy(count,kkk)=total(bo_y(jjj:jjj+spec_av-1,kkk))/spec_av
		tmpboz(count,kkk)=total(bo_z(jjj:jjj+spec_av-1,kkk))/spec_av
		
	endfor
	count=count+1
endfor
store_data,'Bx_fac_wvlt',data={x:time,y:tmpbx,v:bx_wvlt.v,spec:1}
store_data,'By_fac_wvlt',data={x:time,y:tmpby,v:bx_wvlt.v,spec:1}
store_data,'Bz_fac_wvlt',data={x:time,y:tmpbz,v:bx_wvlt.v,spec:1}
store_data,'Vx_fac_wvlt',data={x:time,y:tmpvx,v:bx_wvlt.v,spec:1}
store_data,'Vy_fac_wvlt',data={x:time,y:tmpvy,v:bx_wvlt.v,spec:1}
store_data,'Vz_fac_wvlt',data={x:time,y:tmpvz,v:bx_wvlt.v,spec:1}
store_data,'dens_wvlt',data={x:time,y:tmpd,v:bx_wvlt.v,spec:1}
denso=tmpdo
bo_x=tmpbox
bo_y=tmpboy
bo_z=tmpboz
Te_r=Te_tmp
Ti_r=Ti_tmp
totpoints=red_points

endif
;MAIN LOOP

matspec=make_array(totpoints,n_freq,3,3,/dcomplex)
waveangle=make_array(totpoints,n_freq,/double)
wnx=make_array(totpoints,n_freq,/double)
wny=make_array(totpoints,n_freq,/double)
wnz=make_array(totpoints,n_freq,/double)

bo_x_w=make_array(totpoints,n_freq,/double)
bo_y_w=make_array(totpoints,n_freq,/double)
bo_z_w=make_array(totpoints,n_freq,/double)


bx_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)
by_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)
bz_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)

vx_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)
vy_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)
vz_wnv_wvlt=make_array(totpoints,n_freq,/dcomplex)

modes=make_array(totpoints,n_freq,6,/double)
out=make_array(totpoints,n_freq,/double)

wave_vec_sp=make_array(totpoints,n_freq,6,/dcomplex)
in_vec_sp=make_array(totpoints,n_freq,6,/dcomplex)

;k_location=make_array(totpoints,n_freq,2,/double)
k_angle=make_array(totpoints,n_freq,/double)
;total_fraction=make_array(totpoints,n_freq,/double)
for j=0L,totpoints-1 do begin  

	;CALCULATION OF THE MINIMUM VARIANCE DIRECTION in B AND Wavenormal_Angle
      	
	if keyword_set(alf_test) or keyword_set(fast_test) or keyword_set(slow_test) or keyword_set(wavevector) then begin
		if not keyword_set(wavevector) then begin
			for kk=0,n_freq-1 do begin
				
				temp=fieldrot_fast(bo_x(j,kk),bo_y(j,kk),bo_z(j,kk),kx/sqrt(kx^2+ky^2+kz^2),ky/sqrt(kx^2+ky^2+kz^2),kz/sqrt(kx^2+ky^2+kz^2))
				wnx(j,kk)=temp(0)
				wny(j,kk)=temp(1)
				wnz(j,kk)=temp(2)
			endfor
			waveangle(j,*)=theta_rad 
		endif else begin
			wnx(j,*)=kx/sqrt(kx^2+ky^2+kz^2)
     			wny(j,*)=ky/sqrt(kx^2+ky^2+kz^2)
     			wnz(j,*)=kz/sqrt(kx^2+ky^2+kz^2)
			waveangle(j,*)=ATAN(Sqrt(wnx(j,*)^2+wny(j,*)^2),abs(wnz(j,*))) 
		endelse
		;uncomment the next three lines and comment out the for loop above if taking the mode test short cut from the previous if statement.
		
		;wnx(j,*)=kx/sqrt(kx^2+ky^2+kz^2)
     		;wny(j,*)=ky/sqrt(kx^2+ky^2+kz^2)
     		;wnz(j,*)=kz/sqrt(kx^2+ky^2+kz^2)
		;waveangle(j,*)=theta_rad 
	endif else begin
		;CALCULATION OF THE WAVELET SPECTRAL MATRIX IN B
 
    		matspec(j,*,0,0)=bx_fac_wvlt(j,*)*conj(bx_fac_wvlt(j,*))
    		matspec(j,*,1,0)=bx_fac_wvlt(j,*)*conj(by_fac_wvlt(j,*))
    		matspec(j,*,2,0)=bx_fac_wvlt(j,*)*conj(bz_fac_wvlt(j,*))
    		matspec(j,*,0,1)=by_fac_wvlt(j,*)*conj(bx_fac_wvlt(j,*))
    		matspec(j,*,1,1)=by_fac_wvlt(j,*)*conj(by_fac_wvlt(j,*))
    		matspec(j,*,2,1)=by_fac_wvlt(j,*)*conj(bz_fac_wvlt(j,*))
    		matspec(j,*,0,2)=bz_fac_wvlt(j,*)*conj(bx_fac_wvlt(j,*))
    		matspec(j,*,1,2)=bz_fac_wvlt(j,*)*conj(by_fac_wvlt(j,*))
    		matspec(j,*,2,2)=bz_fac_wvlt(j,*)*conj(bz_fac_wvlt(j,*))

		aaa2=SQRT(IMAGINARY(matspec(j,*,0,1))^2+IMAGINARY(matspec(j,*,0,2))^2+IMAGINARY(matspec(j,*,1,2))^2)
     		wnx(j,*)=IMAGINARY(matspec(j,*,1,2))/aaa2;ABS(IMAGINARY(matspec(j,*,1,2))/aaa2)
     		wny(j,*)=-IMAGINARY(matspec(j,*,0,2))/aaa2;-ABS(IMAGINARY(matspec(j,*,0,2))/aaa2)
     		wnz(j,*)=IMAGINARY(matspec(j,*,0,1))/aaa2
     		
		;setup so min var direction always points in the direction of bo - remember it is degenerate
		test=where(wnz(j,*) LT 0.)		
		wnx(j,test)=-1.0*wnx(j,test)		
		wny(j,test)=-1.0*wny(j,test)
		wnz(j,test)=-1.0*wnz(j,test)
     		
		
		WAVEANGLE(j,*)=ATAN(Sqrt(wnx(j,*)^2+wny(j,*)^2),(wnz(j,*)));because min var set to always be along bo this angle is always between 0 and 90 degrees. 
		for kk=0,n_freq-1 do begin
				;now rotate to get the dc magfile for each bin in the minvar coordiate system
				temp=fieldrot_fast(wnx(j,kk),wny(j,kk),wnz(j,kk),bo_x(j,kk),bo_y(j,kk),bo_z(j,kk))
				bo_x_w(j,kk)=temp(2);now rotate clockwise about y by 90 degrees so that the z component becomes the x component
				bo_y_w(j,kk)=temp(1)
				bo_z_w(j,kk)=-1.0*temp(0)
				
			
			endfor
	
	
	
	
	
	
	endelse

	
	
	;NOW DO THE MODE RECOGNITION
	
	power=make_array(n_freq,/double)

	;define the common block plasma parameters for this time step
	;help,te_r
	;print,j
	scale_factor=1.0;1.0/1.5;should be 1., use different values for referee to show no anisotropies	
	T_i=Ti_r(j)*scale_factor
	T_e=Te_r(j)*scale_factor

	spec_dens_scale_factor=1.;change this factor to change the amount of compression. For example to include a temperature variation a factor of 2 in here would 					 ;correspond to a normalised temperature variation equivalent to the normalised density variation.  

	if keyword_set(theta_minvar) then begin
		for kk=0,n_freq-1 do begin
			;define the common block background magnetic field for this time and frequency. 
			n_i=denso(j,kk)		
			Bo=1.0e-9*sqrt(bo_x(j,kk)^2+bo_y(j,kk)^2+bo_z(j,kk)^2)
			va=Bo/sqrt(mo*mi*n_i*1.0e6)
		        cyc_freq=1.6e-19*bo/1.67e-27
			cs=sqrt((T_e+T_i)*1.6e-19/mi)
			;print,'parameters',va,Bo,T_e(j),T_i(j),n_i
			;define the spectral compnents in fac coords
			bx_fac_spec=bx_fac_wvlt(j,kk)
			by_fac_spec=by_fac_wvlt(j,kk)
			bz_fac_spec=bz_fac_wvlt(j,kk)
			vx_fac_spec=vx_fac_wvlt(j,kk)
			vy_fac_spec=vy_fac_wvlt(j,kk)
			vz_fac_spec=vz_fac_wvlt(j,kk)
			dens_spec=dens_wvlt.y(j,kk)*spec_dens_scale_factor
			rho_zero=rho_o(j,kk)
			wn_x=wnx(j,kk)
			wn_y=wny(j,kk)
			wn_z=wnz(j,kk)
			;state_vec=[1.0/va*vy_wnv_wvlt(j,kk),1.0/va*by_wnv_wvlt(j,kk)/sqrt(mo*rho_o(j,kk)),1.0/va*vx_wnv_wvlt(j,kk),$
					;1.0/va*vz_wnv_wvlt(j,kk),1.0/va*bz_wnv_wvlt(j,kk)/sqrt(mo*rho_o(j,kk)),1.0/va*dens_wvlt.y(j,kk)/n_i*cs]
			if not keyword_set(start_vector) then begin
				kangle_in=waveangle(j,kk)
                        endif else begin
				wn_x=start_vector_array(j,kk,0)
				wn_y=start_vector_array(j,kk,1)
				wn_z=start_vector_array(j,kk,2)
		                kangle_in=ATAN(Sqrt(wn_x^2+wn_y^2),abs(wn_z))
				;print,'kangle',kangle
				;print,'freq',2.0*!DPI*bx_wvlt.v(kk)/cyc_freq
			endelse
	

			;print,state_vec
			;return	
			
			;k_in=1.0e-6		
			;guess=[k_in*cos(kangle_in),k_in*sin(kangle_in)]
			;if mode EQ 'alfven_backward' or mode EQ 'fast_backward' or mode EQ 'slow_backward' then sign=-1. else sign=1.


			;guess=[k_in*cos(kangle_in),sign*abs(k_in*sin(kangle_in))]			
			guess=kangle_in
			;print,'vy',vy_wnv_wvlt(0,0)
			;print,'kangle_in',kangle_in
			
			;power(kk)=total(reform(state_vec*conj(state_vec)))
			;print,'param+m',Bo,T_i,T_e,n_i
			fmin=fluid_kinetic_mode_spec_diff_MHD_angle(guess)
			
			;return

			k_angle_out=guess
        		;print,'frequencies',n_freq,' ind',kk,' power',power(kk)
  			;help,k_out
			;total_fraction(j,kk)=mag_spec_frac
  			k_angle(j,kk)=k_angle_out
  			out(j,kk)=fmin
  			;get the decompoed wave vector amplitudes for each mode
  			wave_vec_sp(j,kk,*)=wave_vec
  			in_vec_sp(j,kk,*)=in_vec
  			;print,'******k and kangle',k_location(kk,0),k_location(kk,1)/!DPI*180.
  			;print,'******kpar and kper',k_out
  			;print,'******fmin',fmin
  			
  			;get the svd estimate for k for each mode
  			if keyword_set(svd) then begin
  			
  			endif
  			
	endfor
	endif else begin
		for kk=0,n_freq-1 do begin
			;define the common block background magnetic field for this time and frequency. 
			n_i=denso(j,kk)		
			Bo=1.0e-9*sqrt(bo_x(j,kk)^2+bo_y(j,kk)^2+bo_z(j,kk)^2)
			va=Bo/sqrt(mo*mi*n_i*1.0e6)
		
			cs=sqrt((T_e+T_i)*1.6e-19/mi)
			


			;define the spectral compnents in fac coords
			;help,bx_fac_wvlt(j,kk)			
			bx_fac_spec=bx_fac_wvlt(j,kk)
			by_fac_spec=by_fac_wvlt(j,kk)
			bz_fac_spec=bz_fac_wvlt(j,kk)
			vx_fac_spec=vx_fac_wvlt(j,kk)
			vy_fac_spec=vy_fac_wvlt(j,kk)
			vz_fac_spec=vz_fac_wvlt(j,kk)
			rho_zero=rho_o(j,kk)
			dens_spec=dens_wvlt.y(j,kk)
			if not keyword_set(start_vector) then begin
				wn_x=wnx(j,kk)
				wn_y=wny(j,kk)
				wn_z=wnz(j,kk)
                        endif else begin
				wn_x=start_vector(0)
				wn_y=start_vector(1)
				wn_z=start_vector(2)
			endelse
			;state_vec=[1.0/va*vy_wnv_wvlt(j,kk),1.0/va*by_wnv_wvlt(j,kk)/sqrt(mo*rho_o(j,kk)),1.0/va*vx_wnv_wvlt(j,kk),$
					;1.0/va*vz_wnv_wvlt(j,kk),1.0/va*bz_wnv_wvlt(j,kk)/sqrt(mo*rho_o(j,kk)),1.0/va*dens_wvlt.y(j,kk)/n_i*cs]
			kangle_in=waveangle(j,kk)
			;k_in=1.0e-6		
			;power(kk)=total(reform(state_vec*conj(state_vec)))
  			;xi=transpose([[k_in*cos(kangle_in)/100.0,0.],[0.,k_in*sin(kangle_in)/100.0]])
  			;
			;xi=(kangle_in)/100.0
			;if mode EQ 'alfven_backward' or mode EQ 'fast_backward' or mode EQ 'slow_backward' then sign=-1. else sign=1.


			;guess=[(k_in*cos(kangle_in)),sign*abs(k_in*sin(kangle_in))]
			guess=kangle_in
			scale=[!DPI/2.]
			;powell,guess,xi,1.0,fmin,'fluid_kinetic_mode_spec_diff_MHD_angle',/double,itmax=1000.
			k_amoeba=amoeba(1.0e-5,scale=scale,P0=guess,function_name='fluid_kinetic_mode_spec_diff_MHD_angle',FUNCTION_VALUE=fval)			
			k_angle_out=k_amoeba
			k_angle(j,kk)=abs(atan(tan(k_angle_out)));get in range 0-Pi/2
  			out(j,kk)=fval[0]
  			;total_fraction(j,kk)=mag_spec_frac
			;print,'angle',k_amoeba
			;print,'min_value',fval[0]
			
			;k_out=guess
        		;print,'frequencies',n_freq,' ind',kk,' power',power(kk)
  			;help,k_out
  			;k_location(j,kk,*)=[sqrt(k_out(0)^2+k_out(1)^2),atan((k_out(1)),(k_out(0)))]
  			;out(j,kk)=fmin
  			;print,'******k and kangle',k_location(kk,0),k_location(kk,1)/!DPI*180.
  			;print,'******kpar and kper',k_out
  			;print,'******fmin',fmin
		endfor
	endelse
	print,'time',time_string(time(j)),' mode ',mode
	;!p.multi=[0,1,3]
	;!p.charsize=2.0
	;plot,frequency,power,/xlog,/ylog,xrange=[min(frequency),max(frequency)]
	;plot,frequency,out(j,*),/xlog,/ylog,xrange=[min(frequency),max(frequency)],yrange=[0.01,1]
	;oplot,frequency,out(j,*),psym=5.0
	
	;plot,frequency,k_angle(j,*),/xlog,xrange=[min(frequency),max(frequency)],thick=2
	;oplot,frequency,waveangle(j,*),color=240,thick=2
	
	;print,'f',frequency
	;print,'p',power

endfor
;create tplot spectral quantities
kill=where(out EQ 0.0)
out(kill)=!values.f_nan

theta_array=k_angle
;theta_array(kill)=!values.f_nan
store_data,mode +'_theta_wave',data={x:time,y:theta_array,v:bx_wvlt.v,spec:1}
store_data,mode +'_comp',data={x:time,y:(1.0-out),v:bx_wvlt.v,spec:1}
;store_data,mode +'_comp_tot_fraction',data={x:time,y:sqrt(1.0-total_fraction),v:bx_wvlt.v,spec:1}

options,mode +'_theta_wave','ylog',1
options,mode +'_comp','ylog',1


options,mode +'_theta_wave','yrange',[min(bx_wvlt.v),max(bx_wvlt.v)]
options,mode +'_comp','yrange',[min(bx_wvlt.v),max(bx_wvlt.v)]


options,mode +'_theta_wave','ystyle',1
options,mode +'_comp','ystyle',1



options,mode +'_theta_wave','ylog',1
options,mode +'_comp','zlog',0
options,mode +'_comp','zrange',[0,1]
tplot,[mode +'_theta_wave',mode +'_comp']

store_data,'out_in',data={x:bx_wvlt.x,y:(bo_z)/abs(bo_z),v:bx_wvlt.v,spec:1}
options,'out_in','zrange',[-1.0,1.0]
options,'out_in','ylog',1

;bandwidth=abs(bx_wvlt.v-shift(bx_wvlt.v,1))
;bandwidth(0)=bandwidth(1)
;bandwidth_sp=make_array(n_elements(time),n_elements(bandwidth),/double)
;for mmm=0,n_elements(time)-1 do bandwidth_sp(mmm,*)=bandwidth
bandwidth_sp=make_array(n_elements(time),n_elements(bx_wvlt.v),value=1.0/(1./sampfreq*2.0*sqrt(!DPI)/2.0))
help,wave_vec_sp
help,wave_vec_sp(*,*,0)*conj(wave_vec_sp(*,*,0))


;spectrum in div b coordinates

store_data,'divb_vy_spec',data={x:time,y:double(in_vec_sp(*,*,0)*conj(in_vec_sp(*,*,0)))/bandwidth_sp,v:bx_wvlt.v,spec:1};m/s^2/Hz
store_data,'divb_by_spec',data={x:time,y:double(in_vec_sp(*,*,1)*conj(in_vec_sp(*,*,1)))/bandwidth_sp,v:bx_wvlt.v,spec:1};T^2/Hz
store_data,'divb_vx_spec',data={x:time,y:double(in_vec_sp(*,*,2)*conj(in_vec_sp(*,*,2)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,'divb_vz_spec',data={x:time,y:double(in_vec_sp(*,*,3)*conj(in_vec_sp(*,*,3)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,'divb_bz_spec',data={x:time,y:double(in_vec_sp(*,*,4)*conj(in_vec_sp(*,*,4)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,'divb_n_spec',data={x:time,y:double(in_vec_sp(*,*,5)*conj(in_vec_sp(*,*,5)))/bandwidth_sp,v:bx_wvlt.v,spec:1};(cm-3)^2/Hz

options,'divb_vy_spec','ylog',1
options,'divb_vy_spec','zlog',1
options,'divb_by_spec','ylog',1
options,'divb_by_spec','zlog',1
options,'divb_vx_spec','ylog',1
options,'divb_vx_spec','zlog',1
options,'divb_vz_spec','ylog',1
options,'divb_vz_spec','zlog',1
options,'divb_bz_spec','ylog',1
options,'divb_bz_spec','zlog',1
options,'divb_n_spec','ylog',1
options,'divb_n_spec','zlog',1


;spectrum of decomposed mode
store_data,mode +'_vy_spec',data={x:time,y:double(wave_vec_sp(*,*,0)*conj(wave_vec_sp(*,*,0)))/bandwidth_sp,v:bx_wvlt.v,spec:1};m/s^2/Hz
store_data,mode +'_by_spec',data={x:time,y:double(wave_vec_sp(*,*,1)*conj(wave_vec_sp(*,*,1)))/bandwidth_sp,v:bx_wvlt.v,spec:1};T^2/Hz
store_data,mode +'_vx_spec',data={x:time,y:double(wave_vec_sp(*,*,2)*conj(wave_vec_sp(*,*,2)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,mode +'_vz_spec',data={x:time,y:double(wave_vec_sp(*,*,3)*conj(wave_vec_sp(*,*,3)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,mode +'_bz_spec',data={x:time,y:double(wave_vec_sp(*,*,4)*conj(wave_vec_sp(*,*,4)))/bandwidth_sp,v:bx_wvlt.v,spec:1}
store_data,mode +'_n_spec',data={x:time,y:double(wave_vec_sp(*,*,5)*conj(wave_vec_sp(*,*,5)))/bandwidth_sp,v:bx_wvlt.v,spec:1};(cm-3)^2/Hz

options,mode +'_vy_spec','ylog',1
options,mode +'_vy_spec','zlog',1
options,mode +'_by_spec','ylog',1
options,mode +'_by_spec','zlog',1
options,mode +'_vx_spec','ylog',1
options,mode +'_vx_spec','zlog',1
options,mode +'_vz_spec','ylog',1
options,mode +'_vz_spec','zlog',1
options,mode +'_bz_spec','ylog',1
options,mode +'_bz_spec','zlog',1
options,mode +'_n_spec','ylog',1
options,mode +'_n_spec','zlog',1

store_data,'bo_x_w',data={x:time,y:bo_x_w,v:bx_wvlt.v,spec:1}
store_data,'bo_y_w',data={x:time,y:bo_y_w,v:bx_wvlt.v,spec:1}
store_data,'bo_z_w',data={x:time,y:bo_z_w,v:bx_wvlt.v,spec:1}


store_data,'wnx',data={x:time,y:wnx,v:bx_wvlt.v,spec:1}
store_data,'wny',data={x:time,y:wny,v:bx_wvlt.v,spec:1}
store_data,'wnz',data={x:time,y:wnz,v:bx_wvlt.v,spec:1}


return
end
 

