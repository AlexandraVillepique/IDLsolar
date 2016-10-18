
;05-08-2008 - A.Andic (Univeristy of Birmingham) - following changes introduced:

	;to start the code first compile it with .r filename.pro and then start with just STANDARD_FITTING  -code itself start where this name is written


;1. added the directory so that initialisation files can be readed in when code is started outside the directory where they are.


;2. path for reading the data  changed

;3. Added some interfaces to make control of the paramethers easier.

;==================================================================================================
; PowFit_Pair routine - A new 'peak-bagging' fitting code for low-l 'sun-as-a-star' data.
; Improvements over older fitting code include the option of directly fitting the l=4 and 5 and
; the option of using a spectral window convolved model to account for gaps in the data, as
; opposed to fitting sidebands.
;==================================================================================================


pro Standard_Fitting_v1


;Define common variables used in main routine and ModefDef function
Common Share, nu_1, nu_2, Pow2, j, N1, N2, NNlow1, NNlow2, T, $
       PowWin_Short, sp, spv, asy, k, lm, mh, ich, sb, cv

answer='y'
;read, answer, prompt='plot profile? (y/n) '

str1='/home/andic/idl/amb_fitting/' ;path where initialization files are

str2='/home/andic/podaci/' ;path for data

fit_sym='n'; -adds the noise in simulated data when value is different than 'n'
;--------------------------------------------------------------------------------------------------
;Program initialization - Read in control file parameters and filenames
;--------------------------------------------------------------------------------------------------

;Determine number of lines in control file. A new line must be entered every time a change is
;made to one of the control parameters.



Ctrl1 = read_ascii(str1+'Control.dat',count=nR)

Ctrl=fltarr(15,nR)

;Read in the control file
close, 1
openr, 1, str1+'Control.dat'
readf, 1, Ctrl
close, 1

;Determine the number of input and output files being used.
nF_in=total(ctrl[13,*])
nF_out=total(ctrl[14,*])


FilesIn=strarr(nF_in)

;Read in the input file names.
close, 2
openr, 2, str1+'FilesIn.dat'
readf, 2, FilesIn
close, 2

FilesOut=strarr(nF_out)

;Read in the output file names.
close, 3
openr, 3, str1+'FilesOut.dat'
readf, 3, FilesOut
close, 3

;Set file name counters to zero
Fi=0
Fo=0

;--------------------------------------------------------------------------------------------------
;Main program loop - used if fitting multiple spectra or if changing control parameters for
;different degree modes.
;--------------------------------------------------------------------------------------------------


for R=0,nR-1 do begin

print,''
print,'***********************************'
print,''
print, 'Processing spectra no: (R)', R
print,''
print,'***********************************'
print,''
;--------------------------------------------------------------------------------------------------
;Read control parameters into more easily tracked variable names
;--------------------------------------------------------------------------------------------------
print,''
print,'----------------------------------------------------------------------------------------------'
print,''
print,'Parameters for calculations:'
print,''

;Number of days in time series.
D=ctrl[0,R]
print,'Number of days in time series:',D

;Cadence of time series in seconds.
cad=ctrl[1,R]
print,'Cadence of time series:',cad

;Which degree (n) to start fitting at.
nlow=ctrl[2,R]

;Which degree (n) to end fitting at.
nhigh=ctrl[3,R]

print,'Fitting will start at the degree (n):',nlow, ' and end at:',nhigh

;Choice on whether to fit l=0/2 pairs (0) or l=1/3 pairs (1).
k=ctrl[4,R]

case k of 
0: print,'Fitting l=0/2 pairs'
1: print,'Fitting l=1/3 pairs'
endcase

;Choice on whether to fit rotational splitting (1) or used a fixed value (0).
sp=ctrl[5,R]

;Value of fixed rotational splitting, ignored if sp=1.
spv=ctrl[6,R]

case sp of 
0: print,'Using fixed value for rotational splitting of:',spv
1: print,'Fitting rotational splitting'
endcase

;Choice on whether to fit asymmetries individually (1) or in pairs (0).
asy=ctrl[7,R]

case asy of 
0: print,'Fitting asymmetries in pairs'
1: print,'Fitting asymetries individually'
endcase

;Choice on whether to also fit l=4 modes (1), l=4 and 5 modes (2) or neither (0)
;when fitting a section of spectrum centred on a l=0/2 pair, ignored if k=1.
lm=ctrl[8,R]

if k eq 0 then begin
	case lm of
	0: print,'Do not fit 4 and 5 modes'
	1: print,'Fitting 4 modes'
	2: print,'Fitting 4 and 5 modes'
	endcase

endif

;Choice on whether to fit inner component height ratios (1) or used fixed a value (0).
mh=ctrl[9,R]

;Value of fixed inner component height ratio (either l=2 or l=3 depending on k), ignored if mh=1
ich=ctrl[10,R]

case mh of 
0: print,'Using fixed value for height ratio of:',ich
1: print,'Fitting inner component height ratios'
endcase

;Choice on whether to fit first order sidebands (1) or not (0).
sb=ctrl[11,R]

case sb of 
0: print,'Do not fit first order sidebands'
1: print,'Fitting first order sidebands'
endcase

;Choice on whether to fit spectral window convolved model (1) or not (0).
cv=ctrl[12,R]

case cv of 
0: print,'Do not fit spectral window convolved model'
1: print,'Fitting spectral window convolved model'
endcase

;Set whether a new input file is to be opened (1) or not (0).
F_in=ctrl[13,R]

case F_in of 
0: print,'Do not open new input file'
1: print,'Open new input file'
endcase

;Set whether a new output file is to be opened (1) or not (0).
F_out=ctrl[14,R]

case F_out of 
0: print,'Do not open new output file'
1: print,'Open new output file'
endcase

print,''
print,'---------------------------------------------------------------------------------------------'
print,''
odgovor='y'
read, odgovor, prompt='Are parameters set as you wish? (y/n) '
print,''
print,''
if odgovor eq 'n' then begin
print,'Please change values in Control.dat so that parameters are as you wish and then compile the program again and restart it.'
print,''
stop
endif


;--------------------------------------------------------------------------------------------------
;Program response to certain control parameter combinations
;--------------------------------------------------------------------------------------------------

;Program ignores lm=1 or 2 if l=1/3 pair (k=1) being fitted
if (k eq 1) and (lm eq 1) then begin
print, ''
print, 'l=4 modes are not fitted when fitting l=1/3 pairs'
endif
if (k eq 1) and (lm eq 2) then begin
print, ''
print, 'l=4 or 5 modes are not fitted when fitting l=1/3 pairs'
endif

;Program ignores preset values height ratios is they are chosen to be fitted
if (mh eq 1) and (k eq 0) then begin
print, ''
print, 'Inner component height ratios are being fitted, preset values in control file ignored'
endif

;Prevent sidebands from being fitted if cv want to fit the spectral window convolved model
if (cv eq 1) and (sb eq 1) then begin
sb=0
print, ''
print, 'You were attempting to fit both a spectral window convolved model and first order'
print, 'sidebands. The option to fit sidebands has therefore been turned off. Please restart'
print, 'the program if you wish to fit sidebands rather than the spectral window convolved model.'
endif




;--------------------------------------------------------------------------------------------------
;Define some constants and arrays
;--------------------------------------------------------------------------------------------------

;Calculate number of points in time series
NN=D*(86400L/cad)

;Calculate resolution of power spectrum
res=1.0e06/(NN*cad)

;Make array of frequencies
nu=dindgen(NN)*res

;Determine length in seconds of time series
T=NN*cad

print,''
print,'*****************************************************'
print,''
print,'Code defined the resolution of the power spectrum as:',res,'[sec]'
print,''
print,'*****************************************************'
print,''
;Setup array for Powers to be read into
Pow=dblarr(NN)

;Setup array for FLAG.con values to be read into ready for use as initial guess values
Guess=dblarr(12,85)


;Counter for input file names
Fi=Fi+F_in

;Counter for output file names.
Fo=Fo+F_out

;--------------------------------------------------------------------------------------------------
;Read in time series (name taken from input filename file)
;--------------------------------------------------------------------------------------------------

; If F_in set to 1 a new input time series will be opened
;if (F_in eq 1) then begin

;Input real or simulated time series - enter path of time series fill here
;DataFull = read_binary(+FilesIn[Fi-1],data_type=4)
;Time series to take window function from - usually use same as input time series
;DataWin = read_binary(+FilesIn[Fi-1],data_type=4)

;endif

;Data=DataFull[0:NN-1]

Data=readfits(str2+FilesIn,header)

;Create window function from time series
Win = (Data ne 0.0)*1.0


if fit_sym eq 'n' then goto, jump
;--------------------------------------------------------------------------------------------------
;Add Background (1/f noise)
;Only needed when fitting simulated data
;--------------------------------------------------------------------------------------------------

NoiseFunc=complexarr(NN)
NoiseFunc[0:NN/2]=sqrt(1.0/nu[0:NN/2])
NoiseFunc[0]=NoiseFunc[1]
NoiseFunc[NN/2+1:NN-1]=reverse(NoiseFunc[1:NN/2-1])

;Scale noise
Fnorm=3000.0
nnorm=floor(Fnorm/res)
Norm=0.25/sqrt(2.0*NN)
NormFactor=NoiseFunc[nnorm]/Norm
NoiseFunc=NoiseFunc/NormFactor

seedR=Fi+1

FNoise=complex(NoiseFunc*randomn(seedR,NN),NoiseFunc*randomn(seedR,NN))

Tnoise=fft(FNoise,/inverse)

Tnoise=real_part(Tnoise)

Data=Data+Tnoise


;--------------------------------------------------------------------------------------------------
;Multiply by window function
;--------------------------------------------------------------------------------------------------

Data=Data*Win

jump:
;--------------------------------------------------------------------------------------------------
;Create power spectrum
;--------------------------------------------------------------------------------------------------

;Function is multiplied by T to giver power per Hz scaling and doubled to include negative side
;of the transform
Pow = (abs(fft(Data)))^2 * T * 2.0

;--------------------------------------------------------------------------------------------------
;Setup output files (names taken from output filename file)
;--------------------------------------------------------------------------------------------------

if (F_out eq 1) then begin

;open files to output data
openw, 4, 'Main_'+FilesOut[Fo-1]

if (mh eq 1) then begin
openw, 5, 'mhratio_'+FilesOut[Fo-1]
endif
if (sb eq 1) then begin
openw, 6, 'sbratio_'+FilesOut[Fo-1]
endif

endif


;--------------------------------------------------------------------------------------------------
;Read in main setup files FLAG.con
;--------------------------------------------------------------------------------------------------

;File is used to define initial guess values for certain parameters and to setup fitting windows
close, 7
openr, 7, str1+'FLAG.con'
readf, 7, Guess, format='(2I8,4G13.6)'
close, 7

;stop, guess[*,0]
;--------------------------------------------------------------------------------------------------
;Start of pair-by-pair fitting loop
;--------------------------------------------------------------------------------------------------

;k=0: l=0/2(4/5) fits
;k=1: l=1/3 fits
nbroj=0.
nraz=nhigh-nlow+1.
aa=dblarr(8+sp+asy+lm*(3+asy+sp)+mh+sb,nraz)
granice=dblarr(2,nraz)
for j=(k*33)+(nlow-7),(k*33)+(nhigh-7) do begin
print,''
print,'Calculating for n of value:',nlow+nbroj
print,''
nbroj=nbroj+1.
;--------------------------------------------------------------------------------------------------
;Input first guess values
;--------------------------------------------------------------------------------------------------

;Some initial guess values are taken from FLAG.con file while others are fixed at set values. Note
;that the natural log of the background, linewidths and heights are taken and varied during the
;fitting process rather than the absolute values.


;for l=0/2(4/5) fitting window
;-----------------------------

if (j le 32) then begin

if lm ge 1 and j le 3 then print, 'Trying to fit l=4/5 modes at too low frequencies. Will only fit l=0/2 pair.'
if lm ge 1 and j le 3 then lm=0ull

A=dblarr(8+sp+asy+lm*(3+asy+sp)+mh+sb)
B=dblarr(8+sp+asy+lm*(3+asy+sp)+mh+sb)
E=dblarr(8+sp+asy+lm*(3+asy+sp)+mh+sb)


;Background
A[0] = 4.5

;Asymmetry (if a single value is fitted)
if (asy eq 0) then begin
A[1] = 0.0
endif

;l=0 Frequency
A[2-asy] = Guess[8,j]

;l=0 Linewidth
A[3-asy] = alog( exp(Guess[10,j])*1000.0 )

;l=0 Height
A[4-asy] = alog( 1.000*(2.0*exp(Guess[11,j]))/(!pi*exp(Guess[10,j])*0.001+(2.0/T)) )

;l=0 Asymmetry
if (asy eq 1) then begin
A[4] = 0.0
endif

;l=2 Frequency
A[5] = Guess[2,j]

;l=2 Rotational Splitting
if (sp eq 1) then begin
A[6] = 0.4
endif

;l=2 Linewidth
A[6+sp] = alog( exp(Guess[ 4,j])*1000.0 )

;l=2 Height
A[7+sp] = alog( 0.504*(2.0*exp(Guess[5,j]))/(!pi*exp(Guess[4,j])*0.001+(2.0/T)) )

;l=2 Asymmetry
if (asy eq 1) then begin
A[8+sp] = 0.0
endif

if (lm ge 1) then begin

;l=4 Frequency
A[8+sp+asy] = Guess[8,j+62]

;l=4 Rotational Splitting
if (sp eq 1) then begin
A[9+sp+asy] = 0.4
endif

;l=4 Linewidth
A[9+(2*sp)+asy] = alog( exp(Guess[10,j+62])*1000.0 )

;l=4 Height
A[10+(2*sp)+asy] = alog( 0.028*(2.0*exp(Guess[11,j+62]))/(!pi*exp(Guess[10,j+62])*0.001+(2.0/T)) )

;l=4 Asymmetry
if (asy eq 1) then begin
A[11+(2*sp)+asy] = 0.0
endif

if (lm eq 2) then begin

;l=5 Frequency
A[11+(2*sp)+(2*asy)] = Guess[2,j+63]

;l=5 Rotational Splitting
A[12+(2*sp)+(2*asy)] = 0.4

;l=5 Linewidth
A[12+(3*sp)+(2*asy)] = alog( exp(Guess[ 4,j+63])*1000.0 )

;l=5 Height
A[13+(3*sp)+(2*asy)] = alog( 0.005*(2.0*exp(Guess[5,j+63]))/(!pi*exp(Guess[4,j+63])*0.001+(2.0/T)) )

;l=5 Asymmetry
if (asy eq 1) then begin
A[14+(3*sp)+(2*asy)] = 0.0
endif

endif
endif

;l=2 inner component fractional height (if fitted)
if (mh eq 1) then begin
A[8+sp+asy+lm*(3+asy+sp)] = 0.55
endif

;fractional sideband height (if fitted)
if (sb eq 1) then begin
A[8+sp+asy+lm*(3+asy+sp)+mh] = 0.05
endif

endif

;for l=1/3 fitting window
;-------------------------

if (j ge 33) then begin

A=dblarr(8+(2*sp)+asy+mh+sb)
B=dblarr(8+(2*sp)+asy+mh+sb)
E=dblarr(8+(2*sp)+asy+mh+sb)


;Backgroundm.mathioudakis@qub.ac.uk
A[0] = 4.5

;Asymmetry (if a single value is fitted)
if (asy eq 0) then begin
A[1] = 0.0
endif

;l=1 Frequency
A[2-asy] = Guess[8,j]

;l=1 Rotational Splitting
if (sp eq 1) then A[3-asy] = 0.4

;l=1 Linewidth
A[3+sp-asy] = alog( exp(Guess[10,j])*1000.0 )

;l=1 Height (maximum power spectral density)
A[4+sp-asy] = alog( 1.000*(2.0*exp(Guess[11,j]))/(!pi*exp(Guess[10,j])*0.001+(2.0/T)) )

;l=1 Asymmetry
if (asy eq 1) then A[4+sp] = 0.0

;l=3 Frequency
A[5+sp] = Guess[2,j]

;l=3 Rotational Splitting
if (sp eq 1) then A[6+sp] = 0.4

;l=3 Linewidth
A[6+(2*sp)] = alog( exp(Guess[ 4,j])*1000.0 )

;l=3 Height (maximum power spectral density)
A[7+(2*sp)] = alog( 0.198*(2.0*exp(Guess[5,j]))/(!pi*exp(Guess[4,j])*0.001+(2.0/T)) )

;l=3 Asymmetry
if (asy eq 1) then A[8+(2*sp)] = 0.0

;l=3 inner component fractional height (if fitted)
if (mh eq 1) then A[8+(2*sp)+asy] = 0.55

;fractional sideband height (if fitted)
if (sb eq 1) then A[8+(2*sp)+asy+mh] = 0.05

endif
;stop, 'A', A

;--------------------------------------------------------------------------------------------------
;Setup of fitting windows using frequencies from FLAG.con file
;--------------------------------------------------------------------------------------------------

;Fitting window sizes are chosen to give a total number of points that has a small sum of prime
;factors, in order to reduce computing time

sumpf=fltarr(100)

N1=floor(24576/(3456/D))

for f=0,99 do begin
N1f=N1-49+f
factor, N1f, ppp, nnn, /quiet
sumpf(f)=total(ppp*nnn)
endfor

min_sumpf=min(sumpf,f)
N1=N1-49+f


if (j le 33) then begin

if (lm eq 0) then begin
N2=floor(13200/(3456/D))
flow1=((Guess[8,j]+Guess[2,j])/2.0)-(N1*res)/2
flow2=((Guess[8,j]+Guess[2,j])/2.0)-(N2*res)/2
endif

if (lm ge 1) then begin
N2=floor(19440/(3456/D))
flow1=((Guess[8,j+62]+Guess[2,j+63])/2.0)-(N1*res)/2
flow2=((Guess[8,j+62]+Guess[2,j+63])/2.0)-(N2*res)/2
endif

endif

if (j ge 33) then begin
N2=floor(16000/(3456/D))
flow1=((Guess[8,j]+Guess[2,j])/2.0)-(N1*res)/2
flow2=((Guess[8,j]+Guess[2,j])/2.0)-(N2*res)/2
endif

;Define start and end of convolution window in points
NNlow1=floor((flow1*T)/1.0e6)

;Define start and end of fitting window in points
NNlow2=floor((flow2*T)/1.0e6)

;Create convolution window cropped power spectrum
Pow1=Pow[NNlow1:NNlow1+(N1-1)]
nu_1=(lindgen(N1)+NNlow1)*res

;Create fitting window cropped power spectrum
Pow2=Pow[NNlow2:NNlow2+(N2-1)]
nu_2=(lindgen(N2)+NNlow2)*res

granice(0,nbroj-1)=NNlow2
granice(1,nbroj-1)=N2

;Make power spectrum of window function
PowWin=(abs(fft(Win)))^2 * T * 2.0
PowWin_Short=complexarr(N1)
PowWin_Short[0:(N1/2)-1]  = PowWin[0:(N1/2)-1]

if ( n_elements(PowWin_Short[(N1/2):N1-1]) eq n_elements(PowWin[NN-(N1/2):NN-1]) ) then begin
PowWin_Short[(N1/2):N1-1] = PowWin[NN-(N1/2):NN-1]
endif
if ( n_elements(PowWin_Short[(N1/2):N1-1]) ne n_elements(PowWin[NN-(N1/2):NN-1]) ) then begin
PowWin_Short[(N1/2):N1-2] = PowWin[NN-(N1/2):NN-1]
endif

;tolerance used by powell rotuine to define end criteria, should be set at at least 1.0e-08 but can
;be reduced to speed up computing time when testing program
ftol=1.0e-08

;set variable to track time of fitting and error calculation
T1 = systime(1)

;initialization matrix for IDL Powell minimisation routine
xi=Identity(n_elements(A))
xi=xi*0.1

;perform fitting
resolve_routine, 'modeldef', /is_function
powell, A, xi, ftol, fmin, 'ModelDef', iter=iter, itmax=500, /double

;--------------------------------------------------------------------------------------------------
;Calculate formal errors
;--------------------------------------------------------------------------------------------------

;Set 'E' column vector to final fitted values
E=A

;Hessian matrix used to estimate the errors on the various fitted parameters.

p=n_elements(A)
derr=dblarr(p)
hess=dblarr(p,p)
hessI=dblarr(p,p)
corr=dblarr(p,p)
sigma=dblarr(p)


;set small changes in parameter values

; l=0/2(4/5)
if (j le 32) then begin

                     derr[0]                         = 1.0e-03  ;Background
if (asy eq 0) then   derr[1]                         = 1.0e-05  ;Asymmetry

                     derr[2-asy]                     = 1.0e-02  ;l=0 Frequency
                     derr[3-asy]                     = 1.0e-03  ;l=0 Width
                     derr[4-asy]                     = 1.0e-03  ;l=0 Height
if (asy eq 1) then   derr[4]                         = 1.0e-05  ;l=0 Asymmetry

                     derr[5]                         = 1.0e-02  ;l=2 Frequency
if (sp eq 1)  then   derr[6]                         = 1.0e-05  ;l=2 Rotational Splitting
                     derr[6+sp]                      = 1.0e-03  ;l=2 Width
                     derr[7+sp]                      = 1.0e-03  ;l=2 Height
if (asy eq 1) then   derr[8+sp]                      = 1.0e-05  ;l=2 Asymmetry

if (lm ge 1)  then begin
                     derr[8+sp+asy]                  = 1.0e-02  ;l=4 Frequency
if (sp eq 1)  then   derr[9+sp+asy]                  = 1.0e-05  ;l=4 Rotational Splitting
                     derr[9+(2*sp)+asy]              = 1.0e-03  ;l=4 Width
                     derr[10+(2*sp)+asy]             = 1.0e-03  ;l=4 Height
if (asy eq 1) then   derr[11+(2*sp)+asy]             = 1.0e-05  ;l=4 Asymmetry

if (lm eq 2)  then begin
                     derr[11+(2*sp)+(2*asy)]         = 1.0e-02  ;l=5 Frequency
if (sp eq 1)  then   derr[12+(2*sp)+(2*asy)]         = 1.0e-05  ;l=5 Rotational Splitting
                     derr[12+(3*sp)+(2*asy)]         = 1.0e-03  ;l=5 Width
                     derr[13+(3*sp)+(2*asy)]         = 1.0e-03  ;l=5 Height
if (asy eq 1) then   derr[14+(3*sp)+(2*asy)]         = 1.0e-05  ;l=5 Asymmetry

endif
endif

if (mh eq 1)  then   derr[8+sp+asy+lm*(3+asy+sp)]    = 1.0e-04  ;l=2 inner component fractional height
if (sb eq 1)  then   derr[8+sp+asy+lm*(3+asy+sp)+mh] = 1.0e-04  ;Fractional sideband height

endif


; l=1/3
if (j ge 33)  then begin

                     derr[0]               = 1.0e-03   ;Background
if (asy eq 0) then   derr[1]               = 1.0e-05   ;Asymmetry
                     derr[2-asy]           = 1.0e-02   ;l=1 Frequency
if (sp eq 1)  then   derr[3-asy]           = 1.0e-05   ;l=1 Rotational Splitting
                     derr[3+sp-asy]        = 1.0e-03   ;l=1 Width
                     derr[4+sp-asy]        = 1.0e-03   ;l=1 Height
if (asy eq 1) then   derr[4+sp]            = 1.0e-05   ;l=1 Asymmetry
                     derr[5+sp]            = 1.0e-02   ;l=3 Frequency
if (sp eq 1)  then   derr[6+sp]            = 1.0e-05   ;l=3 Rotational Splitting
                     derr[6+(2*sp)]        = 1.0e-03   ;l=3 Width
                     derr[7+(2*sp)]        = 1.0e-03   ;l=3 Height
if (asy eq 1) then   derr[8+(2*sp)]        = 1.0e-05   ;l=3 Asymmetry
if (mh eq 1)  then   derr[8+(2*sp)+asy]    = 1.0e-04   ;l=3 inner component fractional height
if (sb eq 1)  then   derr[8+(2*sp)+asy+mh] = 1.0e-04   ;Fractional sideband height

endif

for ii=0,p-1 do begin
for jj=0,p-1 do begin

E[ii]=E[ii]+derr[ii]
E[jj]=E[jj]+derr[jj]
S1=ModelDef(E)

E[jj]=E[jj]-derr[jj]
S2=ModelDef(E)

E[ii]=E[ii]-derr[ii]
E[jj]=E[jj]+derr[jj]
S3=ModelDef(E)

E[jj]=E[jj]-derr[jj]
S4=ModelDef(E)

;construct Hessian matrix
hess[ii,jj]=((S1-S2)-(S3-S4))/(derr[ii]*derr[jj])

endfor
endfor

;for ii=0,p-1 do begin
;for jj=0,ii do begin;
;hess[ii,jj]=hess[ii,jj]
;endfor
;endfor

;invert Hessian matrix
hessI=invert(hess)

for ii=0,p-1 do begin
for jj=0,p-1 do begin
corr[ii,jj]=hessI[ii,jj]/(sqrt(hessI[ii,ii]*hessI[jj,jj]))
endfor
endfor

;for ii=0,p-1 do begin
;for jj=0,ii do begin
;corr[ii,jj]=corr[ii,jj]
;endfor
;endfor

;obtain uncertainties on fitted parameter
for ii=0,p-1 do begin
sigma[ii]=sqrt(abs(hessI[ii,ii]))
endfor

;--------------------------------------------------------------------------------------------------
;print results to files
;--------------------------------------------------------------------------------------------------

print, ''
;print, 'No. of Iterations', Iter
print, 'Time to perform fit and uncertainties:', systime(1) - T1, ' Seconds'
dimenzije=size(A)
krajj=dimenzije(1)-1
aa(0:krajj,nbroj-1.)=A
if answer eq 'n' then goto, jump2
print,''
print,'plots only part fitted so far'
print,''
nbrojac=nbroj-1.+nlow
plot_fitted_profiles, nu, A, nbrojac

;stop, ''
jump2:

if nbroj+nlow lt nhigh then begin
answer1=''
read, answer1, prompt='continue with the fitting of next n?  (y/n) '
if answer1 eq 'n' then goto, jump1
endif

;l=0/2(4/5)
if (j le 32) then begin

PF = dblarr(n_elements(A)*2)

for pp=0,2*(n_elements(A)-1),2 do begin
PF[pp]   = A[pp/2]
PF[pp+1] = sigma[pp/2]
endfor


if (sp eq 0) then begin

if (asy eq 0) then begin

if (lm eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:15], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 1) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[16:17], spv, 0.0, PF[18:21], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:15], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 2) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[16:17], spv, 0.0, PF[18:21], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:15], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j+63], Guess[1,j+63], PF[22:23], spv, 0.0, PF[24:27], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif

endif  ; end asy0

if (asy eq 1) then begin

if (lm eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:17], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 1) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[18:19], spv,0.0, PF[20:25], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:17], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 2) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[18:19], spv,0.0, PF[20:25], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv, 0.0, PF[12:17], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j+63], Guess[1,j+63], PF[26:27], spv,0.0, PF[28:33], PF[0:1], format='(2I,12G16.8)'
endif

endif  ; end asy1

endif  ; end sp0


if (sp eq 1) then begin

if (asy eq 0) then begin

if (lm eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[10:17], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 1) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[18:25], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:17], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 2) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[18:25], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:17], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], 0.0, 0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j+63], Guess[1,j+63], PF[26:33], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif

endif  ; end asy0

if (asy eq 1) then begin

if (lm eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[10:19], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 1) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[20:29], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:19], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
endif
if (lm eq 2) then begin
printf, 4, Guess[6,j+62], Guess[7,j+62], PF[20:29], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j], Guess[1,j], PF[10:19], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], 0.0, 0.0, PF[4:9], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[0,j+63], Guess[1,j+63], PF[30:39], PF[0:1], format='(2I,12G16.8)'
endif

endif  ; end asy

endif  ; end sp1


if (mh eq 1) then begin
printf, 5, Guess[0,j], Guess[1,j], abs(PF[16+2*(asy+sp+lm*(3+asy+sp)):17+2*(asy+sp+lm*(3+asy+sp))]), $
format='(2I,2G16.8)'
endif

if (sb eq 1) then begin
printf, 6, Guess[0,j], Guess[1,j], abs(PF[16+2*(asy+sp+lm*(3+asy+sp)+mh):17+2*(asy+sp+lm*(3+asy+sp)+mh)]), $
format='(2I,2G16.8)'
endif


endif   ; end j<33


;l=1/3
if (j ge 33) then begin

PF = dblarr(n_elements(A)*2)

for pp=0,2*(n_elements(A)-1),2 do begin
PF[pp]   = A[pp/2]
PF[pp+1] = sigma[pp/2]
endfor


if (sp eq 0) then begin

if (asy eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[10:11], spv,0.0, PF[12:15], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:5], spv,0.0, PF[6:9], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif

if (asy eq 1) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[8:9], spv,0.0, PF[10:13], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:3], spv,0.0, PF[4:7], PF[0:1], format='(2I,12G16.8)'
endif

endif ; end sp0


if (sp eq 1) then begin

if (asy eq 0) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[12:19], PF[2:3], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[4:11], PF[2:3], PF[0:1], format='(2I,12G16.8)'
endif

if (asy eq 1) then begin
printf, 4, Guess[0,j], Guess[1,j], PF[12:21], PF[0:1], format='(2I,12G16.8)'
printf, 4, Guess[6,j], Guess[7,j], PF[2:11], PF[0:1], format='(2I,12G16.8)'
endif

endif ; end sp1

if (mh eq 1) then begin
printf, 5, Guess[0,j], Guess[1,j], abs(PF[16+2*((2*sp)+asy):17+2*((2*sp)+asy)]), $
format='(2I,2G16.8)'
endif

if (sb eq 1) then begin
printf, 6, Guess[6,j], Guess[7,j], abs(PF[16+2*((2*sp)+asy+mh):17+2*((2*sp)+asy+mh)]), $
format='(2I,2G16.8)'
endif

endif    ; end j>33

;reset lm value
lm=ctrl[8,R]

;end j loop
endfor

;Make sure output files are closed when instructed by the control file
if (R lt nR-1) then begin
if (ctrl[14,R+1] eq 1) then begin
close, 4
close, 5
close, 6
endif
endif

if answer eq 'n' then goto, jump3
print,''
print,'plots only part fitted so far'
print,''
pom=min(granice(0,*))
pom2=max(granice(0,*),pom3)+granice(1,pom3)

snaga=pow[pom:pom2]
frekvencije=nu[pom:pom2]
lbroj=fix(k)
save,granice,aa,snaga,frekvencije,filename=str2+'granice'+strtrim(lbroj,1)+'.dat'

jump3:

answer1=''
read, answer1, prompt='continue with the fitting of next pair of l?  (y/n) '
if answer1 eq 'n' then goto, jump1
;end R loop
endfor

plot_fitted_profiles2,nlow,nhigh
print,'all calculations done!'
jump1: print, 'skipped to end'

close, 4
close, 5
close, 6


print, 'all done!'

end

;--------------------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------------------
