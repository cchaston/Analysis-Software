;+
;Procedure interp_gap,x,y,index=wb,count=c
;replaces NANs with interpolated values.
;-

pro interp_gap,x,y,verbose=verbose,index=wb,count=c

; first check to make sure the x are ok
w = where(finite(x) eq 0,c)
if c ne 0 then message,'invalid x points'
wb=-1

s = size(/struc,y)
if s.n_dimensions ne 1 then begin
   for i=0,s.dimensions[1]-1 do begin
     foo=y[*,i]
     interp_gap,x,foo,verbose=verbose
     y[*,i] = foo
   endfor
   return
end

ny = n_elements(y)
wb = where(finite(y) eq 0,c)
if c eq 0 then return             ; no bad points

if keyword_set(verbose) then print,'Found ',c,' out of ',ny,' bad points'

wbp = 0 > [wb-1,wb+1] < (ny-1)
wbp = wbp[uniq(wbp,sort(wbp))]
w = where(finite(y[wbp]),c)
if c lt 2 then message,'Not enough valid points"

wbp = wbp[w]

yp = interp(y[wbp],x[wbp],x[wb])
y[wb] = yp

end


