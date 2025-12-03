function fluid_kinetic_mode_spec_diff_MHD_angle,k_angle_in,$
                                      alfven_forward=alfven_forward,$
                                      alfven_backward=alfven_backward,$
                                      fast_foward=fast_forward,$
                                      fast_backward=fast_backward,$
                                      slow_forward=slow_forward,$
                                      slow_backward=slow_backward,$
                                      v_b_n=v_b_n,$
                                      v_p_b=v_p_b
                                      
common fac,bx_fac_spec,by_fac_spec,bz_fac_spec,vx_fac_spec,vy_fac_spec,vz_fac_spec,rho_o,dens_spec 
common sv,state_vec
common wave_mode,mode
common wavenumber,k,kangle
common plasma_params,Bo,T_i,T_e,n_i
common minvar,wn_x,wn_y,wn_z
common coverage,mag_spec_frac
common wave_vec_sp,wave_vec
common in_vec_sp,in_vec
mo=4.0*!DPI*1.0e-7
;this pro dtermines the difference between the observed spectral energy density at a spacecraft frame frequency and that given by
;the filtering for either the shear, fast or slow modes with a specified value of k and wave angle
;correct for error in issi book normalization of 1/h_fast and 1/h_slow
;k=sqrt(kpar_kper(0)^2+kpar_kper(1)^2)
;if mode EQ 'alfven_backward' or mode EQ 'fast_backward' or mode EQ 'slow_backward' then sign=-1. else sign=1.
;kpar_kper(0)=sign*abs(kpar_kper(0))
;kpar_kper(0)=(kpar_kper(0))
;kpar_kper(1)=sign*abs(kpar_kper(1))
;kangle=atan(kpar_kper(1),kpar_kper(0))
kangle=abs(atan(tan(k_angle_in[0])));get as angle between 0 and pi/2.;maybe not the best way since this will always create a peak or a trough at 0,pi/2 when doing the minimization -let this angle var over entire 2pi and see where it goes - tried this and ended up with the same composition for forward and backward modes

;print,'k_angle_in',k_angle_in
;now get the spectral matrix from the data derived state vector
;get the fields and velocities into coord system with x along the estimate for the k direction
	b_rot_d=make_array(3,/double)
	b_rot_i=make_array(3,/double)

	v_rot_d=make_array(3,/double)
	v_rot_i=make_array(3,/double)

	;Need to define wnx,wny,wnz. Use the the minimum variance direction to define the relationship between wnx and wny (i.e the orientation of k in the plane perp to Bo) - this fixes the direction of k in the direction azimuthal to bo. Then use kangle to set the relationship between jkperp and kpara.  
		wnz=1.0
		wn_perp=tan(kangle)*wnz
		wnx=wn_perp*cos(atan(wn_y,wn_x))
		wny=wn_perp*sin(atan(wn_y,wn_x))
		
		;normalize
		wn=sqrt(wnx^2+wny^2+wnz^2)
		
		wnx=wnx/wn
		wny=wny/wn
		wnz=wnz/wn
		;print,wnx,wny,wnz
		;print,wn_x,wn_y,wn_z
		;print,'wnx,wn_x',wnx,wn_x		
		;print,'wny,wn_y',wny,wn_y
		;print,'wnz,wn_z',wnz,wn_z
		;return,1
		temp=fieldrot_fast(wnx,wny,wnz,double(bx_fac_spec),double(by_fac_spec),double(bz_fac_spec))
		b_rot_d(0)=temp(0)
		b_rot_d(1)=temp(1)
		b_rot_d(2)=temp(2)

		temp=fieldrot_fast(wnx,wny,wnz,imaginary(bx_fac_spec),imaginary(by_fac_spec),imaginary(bz_fac_spec))
		b_rot_i(0)=temp(0)
		b_rot_i(1)=temp(1)
		b_rot_i(2)=temp(2)

		;now rotate clockwise about y by 90 degrees so that the z component becomes the x component
		;i.e. z->x, x->-z, y->y
		b_rot_d(*)=[b_rot_d(2),b_rot_d(1),-1.0*b_rot_d(0)]
		b_rot_i(*)=[b_rot_i(2),b_rot_i(1),-1.0*b_rot_i(0)]


		temp=fieldrot_fast(wnx,wny,wnz,double(vx_fac_spec),double(vy_fac_spec),double(vz_fac_spec))
		v_rot_d(0)=temp(0)
		v_rot_d(1)=temp(1)
		v_rot_d(2)=temp(2)

		temp=fieldrot_fast(wnx,wny,wnz,imaginary(vx_fac_spec),imaginary(vy_fac_spec),imaginary(vz_fac_spec))
		v_rot_i(0)=temp(0)
		v_rot_i(1)=temp(1)
		v_rot_i(2)=temp(2)

		;now rotate clockwise about y by 90 degrees so that the z component becomes the x component
		;i.e. z->x, x->-z, y->y
		v_rot_d(*)=[v_rot_d(2),v_rot_d(1),-1.0*v_rot_d(0)]
		v_rot_i(*)=[v_rot_i(2),v_rot_i(1),-1.0*v_rot_i(0)]
	

		;reassemble the rotated spectral data - this is in a coordinate system where X points along the min var direction or k vector and Bo lies in the Y-Z plane 

		bx_wnv_spec=dcomplex(b_rot_d(0),b_rot_i(0))
		by_wnv_spec=dcomplex(b_rot_d(1),b_rot_i(1))
		bz_wnv_spec=dcomplex(b_rot_d(2),b_rot_i(2))

		vx_wnv_spec=dcomplex(v_rot_d(0),v_rot_i(0))
		vy_wnv_spec=dcomplex(v_rot_d(1),v_rot_i(1))
		vz_wnv_spec=dcomplex(v_rot_d(2),v_rot_i(2))
		;fraction of mag spec energy density covered
		;mag_spec_frac=abs(bx_wnv_spec)^2/(abs(bx_wnv_spec)^2+abs(by_wnv_spec)^2+abs(bz_wnv_spec)^2)
		;print,'mag_spec_frac',mag_spec_frac

	
;print,'kangle',kangle
;print,'norm_spec',norm_spec
;print,'state_vec',state_vec
;return,1
;now get the eigenvectors for theh wave modes
m_i=1.67e-27;ion mass in kgs
m_e=9.1e-31;e mass i kg
mo=4.0*!DPI*1.0e-7
qi=1.6e-19;prtn charge

rho_o=m_i*n_i*1.0e6;in kg/m^3
va=Bo/sqrt(mo*rho_o)
cyc_freq=1.6e-19*Bo/m_i
c_s=sqrt((T_e+T_i)*1.6e-19/m_i)
beta=(c_s/va)^2
theta_rad=kangle
;print,'theta_rad',theta_rad
;help,theta_rad
v_fast=(va*sqrt((1.+beta)/2.+0.5*sqrt(1.+beta^2-2.*beta*cos(2.*theta_rad))))
v_slow=(va*sqrt((1.+beta)/2.-0.5*sqrt(1.+beta^2-2.*beta*cos(2.*theta_rad))))

;print,'va','v_fast','v_slow','Bo','rho_o',va,v_fast,v_slow,Bo,rho_o
state_vec=[1.0/va*vy_wnv_spec,1.0/va*by_wnv_spec/sqrt(mo*rho_o),1.0/va*vx_wnv_spec,$
					1.0/va*vz_wnv_spec,1.0/va*bz_wnv_spec/sqrt(mo*rho_o),1.0/va*dens_spec/n_i*c_s]
to_phys_units=[va,va*sqrt(mo*rho_o),va,va,va*sqrt(mo*rho_o),va*n_i/c_s]


in_vec=state_vec*to_phys_units
;norm_state=sqrt(transpose(reform(state_vec))#conj(reform(state_vec)))
;print,'state_vec',[1.0/va*vy_wnv_spec,1.0/va*by_wnv_spec/sqrt(mo*rho_o),1.0/va*vx_wnv_spec,$
;					1.0/va*vz_wnv_spec,1.0/va*bz_wnv_spec/sqrt(mo*rho_o),1.0/va*dens_spec/n_i*c_s]/norm_state(0)
;print,'state_vec',state_vec
spectral_matrix=make_array(6,6,/dcomplex)
spectral_matrix=state_vec#conj(state_vec)

norm_spec=trace(spectral_matrix)

;h_fast=(sqrt(1.+va^2*beta/v_fast^2+sin(theta_rad)^2*((v_fast^2+va^2*cos(theta_rad)^2)/(v_fast^2-va^2*cos(theta_rad)^2))));error in nirmalization from ISSI book
h_fast=(sqrt(1.+va^2*beta/v_fast^2+va^2*sin(theta_rad)^2*((v_fast^2+va^2*cos(theta_rad)^2)/(v_fast^2-va^2*cos(theta_rad)^2)^2)))
;h_slow=(sqrt(1.+va^2*beta/v_slow^2+sin(theta_rad)^2*((v_slow^2+va^2*cos(theta_rad)^2)/(v_slow^2-va^2*cos(theta_rad)^2))));error in nirmalization from ISSI book
h_slow=(sqrt(1.+va^2*beta/v_slow^2+va^2*sin(theta_rad)^2*((v_slow^2+va^2*cos(theta_rad)^2)/(v_slow^2-va^2*cos(theta_rad)^2)^2)))
;print,'h_fast','h_slow',h_fast,h_slow

alf_eigen_forward=             [-1./sqrt(2.),  1./sqrt(2.), 0.0,  0.0,                                                                 0.0,                                                       0.0                     ]
alf_eigen_backward=            [ 1./sqrt(2.),  1./sqrt(2.), 0.0,  0.0,                                                                 0.0,                                                       0.0                     ]
fast_eigen_forward= 1.0/h_fast*[ 0.0,          0.0,         1.0, -va^2*sin(theta_rad)*cos(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2), va*v_fast*sin(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2), va*sqrt(beta)/v_fast    ]
fast_eigen_backward=1.0/h_fast*[ 0.0,          0.0,         1.0, -va^2*sin(theta_rad)*cos(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2),-va*v_fast*sin(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2),-va*sqrt(beta)/v_fast    ]   
slow_eigen_forward= 1.0/h_slow*[ 0.0,          0.0,         1.0, -va^2*sin(theta_rad)*cos(theta_rad)/(v_slow^2-va^2*cos(theta_rad)^2), va*v_slow*sin(theta_rad)/(v_slow^2-va^2*cos(theta_rad)^2), va*sqrt(beta)/v_slow    ]  
slow_eigen_backward=1.0/h_slow*[ 0.0,          0.0,         1.0, -va^2*sin(theta_rad)*cos(theta_rad)/(v_slow^2-va^2*cos(theta_rad)^2),-va*v_slow*sin(theta_rad)/(v_slow^2-va^2*cos(theta_rad)^2),-va*sqrt(beta)/v_slow    ]

;orthogonality check
;print,'ttt',total(fast_eigen_forward*slow_eigen_forward)
;print,'ttt',total(fast_eigen_forward*slow_eigen_backward)
;print,'ttt',total(fast_eigen_backward*slow_eigen_forward)
;print,'ttt',total(fast_eigen_backward*slow_eigen_backward)



;print,'alf_eigen_backward',alf_eigen_backward
;normalize

alf_eigen_forward= alf_eigen_forward;/norm(alf_eigen_forward);no need to normalize since these already have a norm of 1             
alf_eigen_backward= alf_eigen_backward;/norm(alf_eigen_backward)
fast_eigen_forward= fast_eigen_forward;/norm(fast_eigen_forward) 
fast_eigen_backward=fast_eigen_backward;/norm(fast_eigen_backward) 
slow_eigen_forward=slow_eigen_forward;/norm(slow_eigen_forward)  
slow_eigen_backward=slow_eigen_backward;/norm(slow_eigen_backward) 
;print,norm([ 0.0,          0.0,         1.0, -va^2*sin(theta_rad)*cos(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2), va*v_fast*sin(theta_rad)/(v_fast^2-va^2*cos(theta_rad)^2), va*sqrt;(beta)/v_fast    ])
;print,1.0/h_fast^2
;print,norm(alf_eigen_forward),norm(alf_eigen_backward),norm(fast_eigen_forward),norm(fast_eigen_backward),norm(slow_eigen_forward),norm(slow_eigen_backward)
;return,1

;print,mode
;print,theta_rad
;print,v_fast,v_slow,beta,va,t_e,t_i,Bo
;print,alf_eigen_forward            
;print,'alf_eigen_forward',alf_eigen_forward
;print,'fast_eigen_forward',fast_eigen_forward
;return,1
;print,fast_eigen_backward
;print,slow_eigen_forward
;print,slow_eigen_backward

;print,'eg',fast_eigen_forward
;print,'eg',fast_eigen_backward
if mode EQ 'alfven_forward' then begin
  ;now do the filtering for the dispersion solution
   gi_alf_forward=transpose((alf_eigen_forward))#spectral_matrix#conj(reform(alf_eigen_forward))
   diff=(norm_spec-gi_alf_forward)/norm_spec
   hi_alf_forward=total(alf_eigen_forward*state_vec);dot product state_vc and eigen_vec
   wave_vec=alf_eigen_forward*to_phys_units*hi_alf_forward;dot product of eigen_vec and physical units conv. then multiplied by the magnitude of the strte vecor along the eigen vector
   ;print,'gi',gi_alf_forward
  ;print,'hi',hi_alf_forward*conj(hi_alf_forward)
   
   endif    
if mode EQ 'alfven_backward' then begin 
  ;now do the filtering for the dispersion solution
  gi_alf_backward=transpose((alf_eigen_backward))#spectral_matrix#conj(reform(alf_eigen_backward))
  diff=(norm_spec-gi_alf_backward)/norm_spec
  hi_alf_backward=total(alf_eigen_backward*state_vec);dot product state_vc and eigen_vec
  wave_vec=alf_eigen_backward*to_phys_units*hi_alf_backward
  
  ;print,'gi',gi_alf_backward
  ;print,'hi',hi_alf_backward*conj(hi_alf_backward)
endif
if mode EQ 'fast_forward' then begin
  ;now do the filtering for the dispersion solution
  gi_fast_forward=transpose((fast_eigen_forward))#spectral_matrix#conj(reform(fast_eigen_forward))
  diff=(norm_spec-gi_fast_forward)/norm_spec
  hi_fast_forward=total(fast_eigen_forward*state_vec)
  wave_vec=fast_eigen_forward*to_phys_units*hi_fast_forward
endif
if mode EQ 'fast_backward' then begin
  ;now do the filtering for the dispersion solution
  gi_fast_backward=transpose((fast_eigen_backward))#spectral_matrix#conj(reform(fast_eigen_backward))
  diff=(norm_spec-gi_fast_backward)/norm_spec
  hi_fast_backward=total(fast_eigen_backward*state_vec)
  wave_vec=fast_eigen_backward*to_phys_units*hi_fast_backward
endif
if mode EQ 'slow_forward' then begin
  ;now do the filtering for the dispersion solution
  gi_slow_forward=transpose((slow_eigen_forward))#spectral_matrix#conj(reform(slow_eigen_forward))
  diff=(norm_spec-gi_slow_forward)/norm_spec
  hi_slow_forward=total(slow_eigen_forward*state_vec)
  wave_vec=slow_eigen_forward*to_phys_units*hi_slow_forward
endif
if mode EQ 'slow_backward' then begin
  ;now do the filtering for the dispersion solution
  gi_slow_backward=transpose((slow_eigen_backward))#spectral_matrix#conj(reform(slow_eigen_backward))
  diff=(norm_spec-gi_slow_backward)/norm_spec
  hi_slow_backward=total(slow_eigen_backward*state_vec)
  wave_vec=slow_eigen_backward*to_phys_units*hi_slow_backward
endif
if finite(diff) EQ 0 then diff=dcomplex(0,0)

;print
;print,'kangle',kangle
;print,'diff',diff
;print


;print,'param',Bo,T_i,T_e,n_i
;print,'state_vec,',state_vec
return, abs(diff)
end
