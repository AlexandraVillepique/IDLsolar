FUNCTION DIFFER,X
;+
; function differ, call is D=DIFFER(X)
; returns differences between elements in the first dimension only
; 2 and 3-D arrays are handled as series of 1-D's
;-
nd = SIZE(x)
nx = nd(1)
nd = nd(0)
CASE 1 OF
	nd EQ 0: RETURN,x
	nd EQ 1: dx = x(1:*) - x(0:(nx-2))
	nd EQ 2: dx = x(1:*,*) - x(0:(nx-2),*)
	nd EQ 3: dx = x(1:*,*,*) - x(0:(nx-2),*,*)
	ELSE: MESSAGE,'DIFFER is not intended for dimesnions >3'
ENDCASE
RETURN,dx
END		
