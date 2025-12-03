function combine2, quanx,quany
;this pro combine the three component E an B field data for input to
;poynting.pro 

datax=create_struct('time',quanx.x ,'comp1',quanx.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 
datay=create_struct('time',quany.x ,'comp1',quany.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 
;quanz=create_struct('time',quanz.x ,'comp1',quanz.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 

;quan1=create_struct('time',quan1.x ,'comp1',quan1.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 
;quan2=create_struct('time',quan2.x ,'comp1',quan2.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 
;quan3=create_struct('time',quan3.x ,'comp1',quan3.y, 'NCOMP',1,'Data_Name','Mag','Units_Name','','calibrated','') 




thelot=datax

fa_fields_combine,thelot, datay,/add,/interp,delt_t=100.
help,thelot,/st
;fa_fields_combine,thelot, quanz,/svy,/add

;fa_fields_combine,thelot,quan1,/interp,delt_t=100.,/add

;fa_fields_combine,thelot,quan2,/interp,delt_t=100.,/add

;fa_fields_combine,thelot,quan3,/interp,delt_t=100.,/add


return,thelot
end
