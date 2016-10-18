FUNCTION apodisate, r, l
;+
; NAME:
;       APODISATE
; PURPOSE:
;       Generate a square apodisation filter
; CALLING SEQUENCE:
;       Filter = Apodisate( RIM, LENGTH )
; INPUTS:
;       RIM   : Rim width of the filter. This is the width within wich
;               the filter rises from zero to one using a hanning shape
;       LENGTH: Length of the filter (will be length x length square)
;
; OUTPUTS:
;       Filter: square array (float) of dim length^2
;
; MODIFICATION HISTORY:
;       18-Jun-1996  P.Suetterlin, KIS
;-

h = hanning(2*r + 1)
h1 = replicate(1., l)
h1(0) = h(0:r)
h1(l-r-1:*) = h(r:*)
return, h1 # h1

END
