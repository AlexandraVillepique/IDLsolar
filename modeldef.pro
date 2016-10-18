;==================================================================================================
; ModelDef Function - Called by Powell routine in main program and is used to setup the model as
; according to the dictates of certain parameters given in the control file. (For IDL, functions's
; must be defined before the main program that calls them.)
;==================================================================================================

Function ModelDef, A

;define common variables used in main routine and ModefDef function
Common Share, nu_1, nu_2, Pow2, j, N1, N2, NNlow1, NNlow2, T, $
       PowWin_Short, sp, spv, asy, k, lm, mh, ich, sb, cv


;--------------------------------------------------------------------------------------------------
;for l=0/2(4/5) modes
;--------------------------------------------------------------------------------------------------

if (j le 32) then begin


;Transfer from log domain to absolute values in order to fit parameters
A[0]=exp(A[0])
A[3-asy]=exp(A[3-asy])
A[4-asy]=exp(A[4-asy])
A[6+sp]=exp(A[6+sp])
A[7+sp]=exp(A[7+sp])
if (lm ge 1) then begin
A[9+(2*sp)+asy]=exp(A[9+(2*sp)+asy])
A[10+(2*sp)+asy]=exp(A[10+(2*sp)+asy])
endif
if (lm eq 2) then begin
A[12+(3*sp)+(2*asy)]=exp(A[12+(3*sp)+(2*asy)])
A[13+(3*sp)+(2*asy)]=exp(A[13+(3*sp)+(2*asy)])
endif


;l=0

nu0a=A[2-asy]

Z0a=2.0*(nu_1-nu0a)/abs(A[3-asy])

if (asy eq 0) then begin
Q0a=(1.0+A[1]*Z0a)^2+A[1]^2
endif
if (asy eq 1) then begin
Q0a=(1.0+A[4]*Z0a)^2+A[4]^2
endif

L0a=((1.00*abs(A[4-asy]))/(1.0+Z0a^2))*Q0a

if (sb eq 1) then begin

nu0a_SBdn = A[2-asy]-11.574
nu0a_SBup = A[2-asy]+11.574

Z0a_SBdn =2.0*(nu_1-nu0a_SBdn)/abs(A[3-asy])
Z0a_SBup =2.0*(nu_1-nu0a_SBdn)/abs(A[3-asy])

if (asy eq 0) then begin
Q0a_SBdn = (1.0+A[1]*Z0a_SBdn)^2+A[1]^2
Q0a_SBup = (1.0+A[1]*Z0a_SBup)^2+A[1]^2
endif
if (asy eq 1) then begin
Q0a_SBdn = (1.0+A[4]*Z0a_SBdn)^2+A[4]^2
Q0a_SBup = (1.0+A[4]*Z0a_SBup)^2+A[4]^2
endif

L0a_SBdn = ((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[4-asy]))/(1.0+Z0a_SBdn^2))*Q0a_SBdn
L0a_SBup = ((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[4-asy]))/(1.0+Z0a_SBup^2))*Q0a_SBup

endif


;l=2

if (sp eq 0) then begin
nu2a=A[5]-2.0*spv
nu2b=A[5]
nu2c=A[5]+2.0*spv
endif
if (sp eq 1) then begin
nu2a=A[5]-2.0*A[6]
nu2b=A[5]
nu2c=A[5]+2.0*A[6]
endif

Z2a=2.0*(nu_1[*]-nu2a)/abs(A[6+sp])
Z2b=2.0*(nu_1[*]-nu2b)/abs(A[6+sp])
Z2c=2.0*(nu_1[*]-nu2c)/abs(A[6+sp])

if (asy eq 0) then begin
Q2a=(1.0+A[1]*Z2a)^2+A[1]^2
Q2b=(1.0+A[1]*Z2b)^2+A[1]^2
Q2c=(1.0+A[1]*Z2c)^2+A[1]^2
endif
if (asy eq 1) then begin
Q2a=(1.0+A[8+sp]*Z2a)^2+A[8+sp]^2
Q2b=(1.0+A[8+sp]*Z2b)^2+A[8+sp]^2
Q2c=(1.0+A[8+sp]*Z2c)^2+A[8+sp]^2
endif

L2a=((1.00*abs(A[7+sp]))/(1.0+Z2a^2))*Q2a
if (mh eq 0) then begin
L2b=((ich*abs(A[7+sp]))/(1.0+Z2b^2))*Q2b
endif
if (mh eq 1) then begin
L2b=((A[8+sp+asy+lm*(3+asy+sp)]*abs(A[7+sp]))/(1.0+Z2b^2))*Q2b
endif
L2c=((1.00*abs(A[7+sp]))/(1.0+Z2c^2))*Q2c


if (sb eq 1) then begin

if (sp eq 0) then begin
nu2a_SBdn = A[5]-2.0*spv-11.574
nu2a_SBup = A[5]-2.0*spv+11.574
nu2b_SBdn = A[5]-11.574
nu2b_SBup = A[5]+11.574
nu2c_SBdn = A[5]+2.0*spv-11.574
nu2c_SBup = A[5]+2.0*spv+11.574
endif

if (sp eq 1) then begin
nu2a_SBdn = A[5]-2.0*A[6]-11.574
nu2a_SBup = A[5]-2.0*A[6]+11.574
nu2b_SBdn = A[5]-11.574
nu2b_SBup = A[5]+11.574
nu2c_SBdn = A[5]+2.0*A[6]-11.574
nu2c_SBup = A[5]+2.0*A[6]+11.574
endif

Z2a_SBdn = 2.0*(nu_1[*]-nu2a_SBdn)/abs(A[6+sp])
Z2a_SBup = 2.0*(nu_1[*]-nu2a_SBup)/abs(A[6+sp])
Z2b_SBdn = 2.0*(nu_1[*]-nu2b_SBdn)/abs(A[6+sp])
Z2b_SBup = 2.0*(nu_1[*]-nu2b_SBup)/abs(A[6+sp])
Z2c_SBdn = 2.0*(nu_1[*]-nu2c_SBdn)/abs(A[6+sp])
Z2c_SBup = 2.0*(nu_1[*]-nu2c_SBup)/abs(A[6+sp])

if (asy eq 0) then begin
Q2a_SBdn = (1.0+A[1]*Z2a_SBdn)^2+A[1]^2
Q2a_Sbup = (1.0+A[1]*Z2a_SBup)^2+A[1]^2
Q2b_SBdn = (1.0+A[1]*Z2b_SBdn)^2+A[1]^2
Q2b_SBup = (1.0+A[1]*Z2b_SBup)^2+A[1]^2
Q2c_SBdn = (1.0+A[1]*Z2c_SBdn)^2+A[1]^2
Q2c_SBup = (1.0+A[1]*Z2c_SBup)^2+A[1]^2
endif
if (asy eq 1) then begin
Q2a_SBdn = (1.0+A[8+sp]*Z2a_SBdn)^2+A[8+sp]^2
Q2a_Sbup = (1.0+A[8+sp]*Z2a_SBup)^2+A[8+sp]^2
Q2b_SBdn = (1.0+A[8+sp]*Z2b_SBdn)^2+A[8+sp]^2
Q2b_SBup = (1.0+A[8+sp]*Z2b_SBup)^2+A[8+sp]^2
Q2c_SBdn = (1.0+A[8+sp]*Z2c_SBdn)^2+A[8+sp]^2
Q2c_SBup = (1.0+A[8+sp]*Z2c_SBup)^2+A[8+sp]^2
endif

L2a_SBdn =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[7+sp]))/(1.0+Z2a_SBdn^2))*Q2a_SBdn
L2a_SBup =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[7+sp]))/(1.0+Z2a_SBup^2))*Q2a_SBup
if (mh eq 0) then begin
L2b_SBdn =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*ich*abs(A[7+sp]))/(1.0+Z2b_SBdn^2))*Q2b_SBdn
L2b_SBup =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*ich*abs(A[7+sp]))/(1.0+Z2b_SBup^2))*Q2b_SBup
endif
if (mh eq 1) then begin
L2b_SBdn =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*A[8+sp+asy+lm*(3+asy+sp)]*abs(A[7+sp]))/(1.0+Z2b_SBdn^2))*Q2b_SBdn
L2b_SBup =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*A[8+sp+asy+lm*(3+asy+sp)]*abs(A[7+sp]))/(1.0+Z2b_SBup^2))*Q2b_SBup
endif
L2c_SBdn =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[7+sp]))/(1.0+Z2c_SBdn^2))*Q2c_SBdn
L2c_SBup =((abs(A[8+sp+asy+lm*(3+asy+sp)+mh])*1.00*abs(A[7+sp]))/(1.0+Z2c_SBup^2))*Q2c_SBup

endif


;l=4 (outer components only)

if (lm ge 1) then begin

if (sp eq 0) then begin
nu4a=A[8+sp+asy]-4.0*spv
nu4b=A[8+sp+asy]+4.0*spv
endif
if (sp eq 1) then begin
nu4a=A[8+sp+asy]-4.0*A[9+sp+asy]
nu4b=A[8+sp+asy]+4.0*A[9+sp+asy]
endif

Z4a=2.0*(nu_1-nu4a)/abs(A[9+(2*sp)+asy])
Z4b=2.0*(nu_1-nu4b)/abs(A[9+(2*sp)+asy])

if (asy eq 0) then begin
Q4a=(1.0+A[1]*Z4a)^2+A[1]^2
Q4b=(1.0+A[1]*Z4b)^2+A[1]^2
endif
if (asy eq 1) then begin
Q4a=(1.0+A[11+(2*sp)+asy]*Z4a)^2+A[11+(2*sp)+asy]^2
Q4b=(1.0+A[11+(2*sp)+asy]*Z4b)^2+A[11+(2*sp)+asy]^2
endif

L4a=((1.00*abs(A[10+(2*sp)+asy]))/(1.0+Z4a^2))*Q4a
L4b=((1.00*abs(A[10+(2*sp)+asy]))/(1.0+Z4b^2))*Q4b



;l=5 (outer components only)

if (lm eq 2) then begin

if (sp eq 0) then begin
nu5a=A[11+(2*sp)+(2*asy)]-5.0*spv
nu5b=A[11+(2*sp)+(2*asy)]+5.0*spv
endif
if (sp eq 1) then begin
nu5a=A[11+(2*sp)+(2*asy)]-5.0*A[12+(2*sp)+(2*asy)]
nu5b=A[11+(2*sp)+(2*asy)]+5.0*A[12+(2*sp)+(2*asy)]
endif

Z5a=2.0*(nu_1[*]-nu5a)/abs(A[12+(3*sp)+(2*asy)])
Z5b=2.0*(nu_1[*]-nu5b)/abs(A[12+(3*sp)+(2*asy)])

if (asy eq 0) then begin
Q5a=(1.0+A[1]*Z5a)^2+A[1]^2
Q5b=(1.0+A[1]*Z5b)^2+A[1]^2
endif
if (asy eq 1) then begin
Q5a=(1.0+A[14+(3*sp)+(2*asy)]*Z5a)^2+A[14+(3*sp)+(2*asy)]^2
Q5b=(1.0+A[14+(3*sp)+(2*asy)]*Z5b)^2+A[14+(3*sp)+(2*asy)]^2
endif

L5a=((1.00*abs(A[13+(3*sp)+(2*asy)]))/(1.0+Z5a^2))*Q5a
L5b=((1.00*abs(A[13+(3*sp)+(2*asy)]))/(1.0+Z5b^2))*Q5b

endif
endif


;Sum components to give Power Model
if (lm eq 0) and (sb eq 0) then begin
PowModel = L0a + L2a + L2b + L2c + A[0]
endif
if (lm eq 1) and (sb eq 0) then begin
PowModel = L0a + L2a + L2b + L2c + L4a + L4b + A[0]
endif
if (lm eq 2) and (sb eq 0) then begin
PowModel = L0a + L2a + L2b + L2c + L4a + L4b + L5a + L5b + A[0]
endif

if (lm eq 0) and (sb eq 1) then begin
PowModel = L0a + L0a_SBdn + L0a_SBup + $
           L2a + L2a_SBdn + L2a_SBup + $
           L2b + L2b_SBdn + L2b_SBup + $
           L2c + L2c_SBdn + L2c_SBup + A[0]
endif
if (lm eq 1) and (sb eq 1) then begin
PowModel = L0a + L0a_SBdn + L0a_SBup + $
           L2a + L2a_SBdn + L2a_SBup + $
           L2b + L2b_SBdn + L2b_SBup + $
           L2c + L2c_SBdn + L2c_SBup + $
           L4a + L4b + A[0]
endif
if (lm eq 2) and (sb eq 1) then begin
PowModel = L0a + L0a_SBdn + L0a_SBup + $
           L2a + L2a_SBdn + L2a_SBup + $
           L2b + L2b_SBdn + L2b_SBup + $
           L2c + L2c_SBdn + L2c_SBup + $
           L4a + L4b + L5a + L5b + A[0]
endif

;Transfer back to log domain so the log of certain values (background, linewidth and heights) are
;varied rather than the absolute values.
A[0]=alog(A[0])
A[3-asy]=alog(A[3-asy])
A[4-asy]=alog(A[4-asy])
A[6+sp]=alog(A[6+sp])
A[7+sp]=alog(A[7+sp])
if (lm ge 1) then begin
A[9+(2*sp)+asy]=alog(A[9+(2*sp)+asy])
A[10+(2*sp)+asy]=alog(A[10+(2*sp)+asy])
endif
if (lm eq 2) then begin
A[12+(3*sp)+(2*asy)]=alog(A[12+(3*sp)+(2*asy)])
A[13+(3*sp)+(2*asy)]=alog(A[13+(3*sp)+(2*asy)])
endif

endif


;--------------------------------------------------------------------------------------------------
;for l=1/3 modes
;--------------------------------------------------------------------------------------------------

if (j ge 33) then begin

;Transfer from log domain to absolute values in order to fit parameters
A[0]=exp(A[0])
A[3+sp-asy]=exp(A[3+sp-asy])
A[4+sp-asy]=exp(A[4+sp-asy])
A[6+(2*sp)]=exp(A[6+(2*sp)])
A[7+(2*sp)]=exp(A[7+(2*sp)])

;l=1

if (sp eq 0) then begin
nu1a=A[2-asy]-spv
nu1b=A[2-asy]+spv
endif
if (sp eq 1) then begin
nu1a=A[2-asy]-A[3-asy]
nu1b=A[2-asy]+A[3-asy]
endif

Z1a=2.0*(nu_1-nu1a)/abs(A[3+sp-asy])
Z1b=2.0*(nu_1-nu1b)/abs(A[3+sp-asy])

if (asy eq 0) then begin
Q1a=(1.0+A[1]*Z1a)^2+A[1]^2
Q1b=(1.0+A[1]*Z1b)^2+A[1]^2
endif
if (asy eq 1) then begin
Q1a=(1.0+A[4+sp]*Z1a)^2+A[4+sp]^2
Q1b=(1.0+A[4+sp]*Z1b)^2+A[4+sp]^2
endif

L1a=((1.00*abs(A[4+sp-asy]))/(1.0+Z1a^2))*Q1a
L1b=((1.00*abs(A[4+sp-asy]))/(1.0+Z1b^2))*Q1b


if (sb eq 1) then begin

if (sp eq 0) then begin
nu1a_SBdn = A[2-asy]-spv-11.574
nu1a_SBup = A[2-asy]-spv+11.574
nu1b_SBdn = A[2-asy]+spv-11.574
nu1b_SBup = A[2-asy]+spv+11.574
endif
if (sp eq 1) then begin
nu1a_SBdn = A[2-asy]-A[3-asy]-11.574
nu1a_SBup = A[2-asy]-A[3-asy]+11.574
nu1b_SBdn = A[2-asy]+A[3-asy]-11.574
nu1b_SBup = A[2-asy]+A[3-asy]+11.574
endif

Z1a_SBdn = 2.0*(nu_1-nu1a_SBdn)/abs(A[3+sp-asy])
Z1a_SBup = 2.0*(nu_1-nu1a_SBup)/abs(A[3+sp-asy])
Z1b_SBdn = 2.0*(nu_1-nu1b_SBdn)/abs(A[3+sp-asy])
Z1b_SBup = 2.0*(nu_1-nu1b_SBup)/abs(A[3+sp-asy])

if (asy eq 0) then begin
Q1a_SBdn = (1.0+A[1]*Z1a_SBdn)^2+A[1]^2
Q1a_SBup = (1.0+A[1]*Z1a_SBup)^2+A[1]^2
Q1b_SBdn = (1.0+A[1]*Z1b_SBdn)^2+A[1]^2
Q1b_SBup = (1.0+A[1]*Z1b_SBup)^2+A[1]^2
endif
if (asy eq 1) then begin
Q1a_SBdn = (1.0+A[4+sp]*Z1a_SBdn)^2+A[5]^2
Q1a_SBup = (1.0+A[4+sp]*Z1a_SBup)^2+A[5]^2
Q1b_SBdn = (1.0+A[4+sp]*Z1b_SBdn)^2+A[5]^2
Q1b_SBup = (1.0+A[4+sp]*Z1b_SBup)^2+A[5]^2
endif

L1a_SBdn = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[4+sp-asy]))/(1.0+Z1a_SBdn^2))*Q1a_SBdn
L1a_SBup = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[4+sp-asy]))/(1.0+Z1a_SBup^2))*Q1a_SBup
L1b_SBdn = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[4+sp-asy]))/(1.0+Z1b_SBdn^2))*Q1b_SBdn
L1b_SBup = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[4+sp-asy]))/(1.0+Z1b_SBup^2))*Q1b_SBup

endif


;l=3

if (sp eq 0) then begin
nu3a=A[5+sp]-3.0*spv
nu3b=A[5+sp]-1.0*spv
nu3c=A[5+sp]+1.0*spv
nu3d=A[5+sp]+3.0*spv
endif
if (sp eq 1) then begin
nu3a=A[5+sp]-3.0*A[6+sp]
nu3b=A[5+sp]-1.0*A[6+sp]
nu3c=A[5+sp]+1.0*A[6+sp]
nu3d=A[5+sp]+3.0*A[6+sp]
endif

Z3a=2.0*(nu_1-nu3a)/abs(A[6+(2*sp)])
Z3b=2.0*(nu_1-nu3b)/abs(A[6+(2*sp)])
Z3c=2.0*(nu_1-nu3c)/abs(A[6+(2*sp)])
Z3d=2.0*(nu_1-nu3d)/abs(A[6+(2*sp)])

if (asy eq 0) then begin
Q3a=(1.0+A[1]*Z3a)^2+A[1]^2
Q3b=(1.0+A[1]*Z3b)^2+A[1]^2
Q3c=(1.0+A[1]*Z3c)^2+A[1]^2
Q3d=(1.0+A[1]*Z3d)^2+A[1]^2
endif
if (asy eq 1) then begin
Q3a=(1.0+A[8+(2*sp)]*Z3a)^2+A[8+(2*sp)]^2
Q3b=(1.0+A[8+(2*sp)]*Z3b)^2+A[8+(2*sp)]^2
Q3c=(1.0+A[8+(2*sp)]*Z3c)^2+A[8+(2*sp)]^2
Q3d=(1.0+A[8+(2*sp)]*Z3d)^2+A[8+(2*sp)]^2
endif

L3a=((1.00*abs(A[7+(2*sp)]))/(1.0+Z3a^2))*Q3a
if (mh eq 0) then begin
L3b=((ich*abs(A[7+(2*sp)]))/(1.0+Z3b^2))*Q3b
L3c=((ich*abs(A[7+(2*sp)]))/(1.0+Z3c^2))*Q3c
endif
if (mh eq 1) then begin
L3b=((A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3b^2))*Q3b
L3c=((A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3c^2))*Q3c
endif
L3d=((1.00*abs(A[7+(2*sp)]))/(1.0+Z3d^2))*Q3d

if (sb eq 1) then begin

if (sp eq 0) then begin
nu3a_SBdn = A[5+sp]-3.0*spv-11.574
nu3a_SBup = A[5+sp]-3.0*spv+11.574
nu3b_SBdn = A[5+sp]-1.0*spv-11.574
nu3b_SBup = A[5+sp]-1.0*spv+11.574
nu3c_SBdn = A[5+sp]+1.0*spv-11.574
nu3c_SBup = A[5+sp]+1.0*spv+11.574
nu3d_SBdn = A[5+sp]+3.0*spv-11.574
nu3d_SBup = A[5+sp]+3.0*spv+11.574
endif
if (sp eq 1) then begin
nu3a_SBdn = A[5+sp]-3.0*A[6+sp]-11.574
nu3a_SBup = A[5+sp]-3.0*A[6+sp]+11.574
nu3b_SBdn = A[5+sp]-1.0*A[6+sp]-11.574
nu3b_SBup = A[5+sp]-1.0*A[6+sp]+11.574
nu3c_SBdn = A[5+sp]+1.0*A[6+sp]-11.574
nu3c_SBup = A[5+sp]+1.0*A[6+sp]+11.574
nu3d_SBdn = A[5+sp]+3.0*A[6+sp]-11.574
nu3d_SBup = A[5+sp]+3.0*A[6+sp]+11.574
endif

Z3a_SBdn = 2.0*(nu_1-nu3a_SBdn)/abs(A[6+(2*sp)])
Z3a_SBup = 2.0*(nu_1-nu3a_SBup)/abs(A[6+(2*sp)])
Z3b_SBdn = 2.0*(nu_1-nu3b_SBdn)/abs(A[6+(2*sp)])
Z3b_SBup = 2.0*(nu_1-nu3b_SBup)/abs(A[6+(2*sp)])
Z3c_SBdn = 2.0*(nu_1-nu3c_SBdn)/abs(A[6+(2*sp)])
Z3c_SBup = 2.0*(nu_1-nu3c_SBup)/abs(A[6+(2*sp)])
Z3d_SBdn = 2.0*(nu_1-nu3d_SBdn)/abs(A[6+(2*sp)])
Z3d_SBup = 2.0*(nu_1-nu3d_SBup)/abs(A[6+(2*sp)])

if (asy eq 0) then begin
Q3a_SBdn = (1.0+A[1]*Z3a_SBdn)^2+A[1]^2
Q3a_SBup = (1.0+A[1]*Z3a_SBup)^2+A[1]^2
Q3b_SBdn = (1.0+A[1]*Z3b_SBdn)^2+A[1]^2
Q3b_SBup = (1.0+A[1]*Z3b_SBup)^2+A[1]^2
Q3c_SBdn = (1.0+A[1]*Z3c_SBdn)^2+A[1]^2
Q3c_SBup = (1.0+A[1]*Z3c_SBup)^2+A[1]^2
Q3d_SBdn = (1.0+A[1]*Z3d_SBdn)^2+A[1]^2
Q3d_SBup = (1.0+A[1]*Z3d_SBup)^2+A[1]^2
endif
if (asy eq 1) then begin
Q3a_SBdn = (1.0+A[8+(2*sp)]*Z3a_SBdn)^2+A[8+(2*sp)]^2
Q3a_SBup = (1.0+A[8+(2*sp)]*Z3a_SBup)^2+A[8+(2*sp)]^2
Q3b_SBdn = (1.0+A[8+(2*sp)]*Z3b_SBdn)^2+A[8+(2*sp)]^2
Q3b_SBup = (1.0+A[8+(2*sp)]*Z3b_SBup)^2+A[8+(2*sp)]^2
Q3c_SBdn = (1.0+A[8+(2*sp)]*Z3c_SBdn)^2+A[8+(2*sp)]^2
Q3c_SBup = (1.0+A[8+(2*sp)]*Z3c_SBup)^2+A[8+(2*sp)]^2
Q3d_SBdn = (1.0+A[8+(2*sp)]*Z3d_SBdn)^2+A[8+(2*sp)]^2
Q3d_SBup = (1.0+A[8+(2*sp)]*Z3d_SBup)^2+A[8+(2*sp)]^2
endif

L3a_SBdn = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[7+(2*sp)]))/(1.0+Z3a_SBdn^2))*Q3a_SBdn
L3a_SBup = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[7+(2*sp)]))/(1.0+Z3a_SBup^2))*Q3a_SBup
if (mh eq 0) then begin
L3b_SBdn = ((abs(A[8+(2*sp)+asy+mh])*ich*abs(A[7+(2*sp)]))/(1.0+Z3b_SBdn^2))*Q3b_SBdn
L3b_SBup = ((abs(A[8+(2*sp)+asy+mh])*ich*abs(A[7+(2*sp)]))/(1.0+Z3b_SBup^2))*Q3b_SBup
L3c_SBdn = ((abs(A[8+(2*sp)+asy+mh])*ich*abs(A[7+(2*sp)]))/(1.0+Z3c_SBdn^2))*Q3c_SBdn
L3c_SBup = ((abs(A[8+(2*sp)+asy+mh])*ich*abs(A[7+(2*sp)]))/(1.0+Z3c_SBup^2))*Q3c_SBup
endif
if (mh eq 1) then begin
L3b_SBdn = ((abs(A[8+(2*sp)+asy+mh])*A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3b_SBdn^2))*Q3b_SBdn
L3b_SBup = ((abs(A[8+(2*sp)+asy+mh])*A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3b_SBup^2))*Q3b_SBup
L3c_SBdn = ((abs(A[8+(2*sp)+asy+mh])*A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3c_SBdn^2))*Q3c_SBdn
L3c_SBup = ((abs(A[8+(2*sp)+asy+mh])*A[8+(2*sp)+asy]*abs(A[7+(2*sp)]))/(1.0+Z3c_SBup^2))*Q3c_SBup
endif
L3d_SBdn = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[7+(2*sp)]))/(1.0+Z3d_SBdn^2))*Q3d_SBdn
L3d_SBup = ((abs(A[8+(2*sp)+asy+mh])*1.00*abs(A[7+(2*sp)]))/(1.0+Z3d_SBup^2))*Q3d_SBup

endif


;Sum components to give Power Model
if (sb eq 0) then begin
PowModel = L1a + L1b + L3a + L3b + L3c + L3d + A[0]
endif

if (sb eq 1) then begin
PowModel = L1a + L1a_SBdn + L1a_SBup + $
           L1b + L1b_SBdn + L1b_SBup + $
           L3a + L3a_SBdn + L3a_SBup + $
           L3b + L3b_SBdn + L3b_SBup + $
           L3c + L3c_SBdn + L3c_SBup + $
           L3d + L3d_SBdn + L3d_SBup + A[0]
endif

;Transfer back to log domain so the log of certain values (background, linewidth and heights) are
;varied rather than the absolute values.
A[0]=alog(A[0])
A[3+sp-asy]=alog(A[3+sp-asy])
A[4+sp-asy]=alog(A[4+sp-asy])
A[6+(2*sp)]=alog(A[6+(2*sp)])
A[7+(2*sp)]=alog(A[7+(2*sp)])


endif

;--------------------------------------------------------------------------------------------------

if (cv eq 0) then begin
PowModel_Fit=PowModel[NNlow2-NNlow1:(NNlow2-NNlow1)+(N2-1)]
S=total(alog(PowModel_Fit)+Pow2/PowModel_Fit)
endif

if (cv eq 1) then begin
Temp=(fft(PowModel,/inverse)/(2.0*T))*(fft(PowWin_short,/inverse)/(2.0*T))
PowModel_Convol=fft(Temp)*2.0*T
PowModel_Convol_Fit=PowModel_Convol[NNlow2-NNlow1:(NNlow2-NNlow1)+(N2-1)]
S=total(alog(PowModel_Convol_Fit)+Pow2/PowModel_Convol_fit)
endif


if (A[0] lt 0.0) then S=1.0e6
if (k eq 0) then begin
if (lm ge 1) then begin
if (exp(A[10+(2*sp)+asy]) gt 3.0*exp(A[4-asy])*0.028) then S=1.0e6
if (lm eq 2) then begin
if (exp(A[13+(3*sp)+(2*asy)]) gt 3.0*exp(A[4-asy])*0.005) then S=1.0e6
endif
endif
endif


Return, S

end

;--------------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------------

