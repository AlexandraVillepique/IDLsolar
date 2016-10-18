FUNCTION Com, data
;+
; NAME:
;       COM
; PURPOSE:
;       Find CenterOfMass of a one or two dimensional array
; CATEGORY:
;       
; CALLING SEQUENCE:
;       RESULT=COM(DATA)
; INPUTS:
;       DATA : 1 or 2-dim Array of data
; OUTPUTS:
;       1 or 2 element vector with index of the COM
; PROCEDURE:
;       
; MODIFICATION HISTORY:
;       14-Jul-1992  P.Suetterlin, KIS
;-

on_error, 2

s = size(data)

n = s(0)
sx = s(1)
sy = s(2)
sz = s(3)
rx = findgen(sx)

td = total(data)

IF n EQ 1 THEN return, total(data*rx)/td
ry = transpose(findgen(sy)) 
IF n EQ 2 THEN BEGIN
    x1 = total(data*rebin(rx, sx, sy))/td
    y1 = total(data*rebin(ry, sx, sy))/td
    return, [x1, y1]
ENDIF
IF n EQ 3 THEN BEGIN
    rz = fltarr(1, 1, sz) & rz(*) = findgen(sz)
    x1 = total(data*rebin(rx, sx, sy, sz))/td
    y1 = total(data*rebin(ry, sx, sy, sz))/td
    z1 = total(data*rebin(rz, sx, sy, sz))/td
    return, [x1, y1, z1]
ENDIF
END

