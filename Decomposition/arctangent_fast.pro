function arctangent_fast, upper, lower

;If (upper GE 0.0) then begin
;                     theta=Atan(Upper,lower)
                   
;endif else begin
      ;   theta=!DPI+(!DPI+Atan(upper,lower))
;endelse

theta=2.0*!DPI*(upper LT 0)+Atan(upper,lower)

return, theta
end
     
