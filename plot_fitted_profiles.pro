pro plot_fitted_profiles, nu, A, nbrojac

Common Share, nu_1, nu_2, Pow2, j, N1, N2, NNlow1, NNlow2, T, $
       PowWin_Short, sp, spv, asy, k, lm, mh, ich, sb, cv
platforma=0 ;set 0 for linux platform and 1 for windows platform
lbrojac=k
case platforma of
	0: begin 
	loadct,3
	window,0,title='Result of fitting',xsize=900,ysize=500
	end
	1: window,0
endcase
plot, nu[NNlow2:NNlow2+(N2-1)], pow2, xtitle='Frequency [micro Hz]', ytitle='Power'
print,'---------------------------------------------------------------------'
print,''

if j le 32 then begin
	print, 'l=0/2 pair	has been fit - see plot for result, frequencies were', A[2-asy], '	and', A[5]
	if asy eq 0 then print, 'a single value of the asymmetry has been fitted for the pair'
	if asy eq 1 then print, 'the asymmetry has been fitted seperately for each mode'
	if sp eq 0 then print, 'the rotational splitting has been fixed at 0.4!4l!3Hz'
	if sp eq 1 then print, 'the rotational splitting has been fitted for each mode'
	if mh eq 0 then print, 'the height of the inner component of the l=2 mode was fixed with a fractional height of ', ich
	if mh eq 1 then print, 'the height of the inner component of the l=2 mode was fitted'
	if lm eq 0 then print, 'only l=0 and 2 modes were fitted'
	if lm eq 1 then print, 'the l=4 mode was also fitted at a frequency', A[8+sp+asy]
	if lm eq 2 then print, 'the l=4 and 5 modes were also fitted at frequencies', A[8+sp+asy], '	and', A[11+2*sp+2*asy]

	;l=0 mode
	part1=2*(nu[NNlow2:NNlow2+(N2-1)]-A[2-asy])/exp(A[3-asy])
	if asy eq 0 then model_l0=(exp(A[4-asy])/(1+part1^2))*((1+A[1]*part1)^2+A[1]^2)
	if asy eq 1 then model_l0=(exp(A[4-asy])/(1+part1^2))*((1+A[4]*part1)^2+A[4]^2)

	;l=2 mode
	if sp eq 0ull then begin 
		freq_neg2=A[5]-0.8
		freq_pos2=A[5]+0.8
	endif
	if sp eq 1ull then begin
		freq_neg2=A[5]-2*A[6]
		freq_pos2=A[5]+2*A[6]
	endif
	part_neg2=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_neg2)/exp(A[6+sp])
	part_pos2=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_pos2)/exp(A[6+sp])
	part_m0=2*(nu[NNlow2:NNlow2+(N2-1)]-A[5])/exp(A[6+sp])

	if asy eq 0 then begin 
		model_l2_neg2=(exp(A[7+sp])/(1+part_neg2^2))*((1+A[1]*part_neg2)^2+A[1]^2)
		model_l2_pos2=(exp(A[7+sp])/(1+part_pos2^2))*((1+A[1]*part_pos2)^2+A[1]^2)
		if mh eq 0 then model_l2_m0=(ich*exp(A[7+sp])/(1+part_m0^2))*((1+A[1]*part_m0)^2+A[1]^2)
		if mh eq 1 then model_l2_m0=(a[8+sp+asy+lm*(3+asy+sp)]*exp(A[7+sp])/(1+part_m0^2))*((1+A[1]*part_m0)^2+A[1]^2)
	endif
	if asy eq 1 then begin
		model_l2_neg2=(exp(A[7+sp])/(1+part_neg2^2))*((1+A[8+sp]*part_neg2)^2+A[8+sp]^2)
		model_l2_pos2=(exp(A[7+sp])/(1+part_pos2^2))*((1+A[8+sp]*part_pos2)^2+A[8+sp]^2)
		if mh eq 0 then model_l2_m0=(ich*str1+exp(A[7+sp])/(1+part_m0^2))*((1+A[8+sp]*part_m0)^2+A[8+sp]^2)
		if mh eq 1 then model_l2_m0=(a[8+sp+asy+lm*(3+asy+sp)]*exp(A[7+sp])/(1+part_m0^2))*((1+A[8+sp]*part_m0)^2+A[8+sp]^2)
	endif

	;sum l=0 and l=2 modes
	model_plot=model_l0+model_l2_neg2+model_l2_m0+model_l2_pos2+exp(A[0])

	;if fitted l=4 modes, note only the outer components are fitted as they are the only ones large enough to see
	if lm ge 1 then begin
	if sp eq 0ull then begin 
		freq_4neg4=A[8+sp+asy]-1.6
		freq_4pos4=A[8+sp+asy]+1.6
	endif
	if sp eq 1ull then begin 
		freq_4neg4=A[8+sp+asy]-4*A[9+sp+asy]
		freq_4pos4=A[8+sp+asy]+4*A[9+sp+asy]
	endif

	part_4neg4=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_4neg4)/exp(A[9+2*sp+asy])
	part_4pos4=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_4pos4)/exp(A[9+2*sp+asy])

	if asy eq 0 then begin
		model_l4_neg4=(exp(A[10+2*sp+asy])/(1+part_4neg4^2))*((1+A[1]*part_4neg4)^2+A[1]^2)
		model_l4_pos4=(exp(A[10+2*sp+asy])/(1+part_4pos4^2))*((1+A[1]*part_4pos4)^2+A[1]^2)
	endif
	if asy eq 1 then begin 
		model_l4_neg4=(exp(A[10+2*sp+asy])/(1+part_4neg4^2))*((1+A[11+2*sp+asy]*part_4neg4)^2+A[11+2*sp+asy]^2)
		model_l4_pos4=(exp(A[10+2*sp+asy])/(1+part_4pos4^2))*((1+A[11+2*sp+asy]*part_4pos4)^2+A[11+2*sp+asy]^2)
	endif

	;add l=4 mode to sum
	model_plot=model_plot+model_l4_neg4+model_l4_pos4

	;if fitted l=5 modes, note only the outer components are fitted as they are the only ones large enough to see
	if lm ge 2 then begin
		if sp eq 0ull then begin
			freq_5neg5=A[11+2*sp+2*asy]-2.0
			freq_5pos5=A[11+2*sp+2*asy]+2.0
		endif
		if sp eq 1ull then begin 
			freq_5neg5=A[11+2*sp+2*asy]-5*A[14+3*sp+2*asy]
			freq_5pos5=A[11+2*sp+2*asy]+5*A[14+3*sp+2*asy]
		endif
		part_5neg5=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_5neg5)/exp(A[12+3*sp+2*asy])
		part_5pos5=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_5pos5)/exp(A[12+3*sp+2*asy])
		if asy eq 0 then begin 
			model_l5_neg5=(exp(A[13+3*sp+2*asy])/(1+part_5neg5^2))*((1+A[1]*part_5neg5)^2+A[1]^2)
			model_l5_pos5=(exp(A[13+3*sp+2*asy])/(1+part_5pos5^2))*((1+A[1]*part_5pos5)^2+A[1]^2)
		endif
		if asy eq 1 then begin 
			model_l5_neg5=(exp(A[13+3*sp+2*asy])/(1+part_5neg5^2))*((1+A[14+3*sp+2*asy]*part_5neg5)^2+A[14+3*sp+2*asy]^2)
			 model_l5_pos5=(exp(A[13+3*sp+2*asy])/(1+part_5pos5^2))*((1+A[14+3*sp+2*asy]*part_5pos5)^2+A[14+3*sp+2*asy]^2)
		endif
		model_plot=model_plot+model_l5_neg5+model_l5_pos5
	endif

	endif		;end for plotting l=4 and 5 modes

endif		;end for l=2/0 pair

;if fitted l=1/3 pair
if j gt 32 then begin
	print, 'l=1/3 pair	has been fit - see plot for result, frequencies were', A[2-asy], '	and', A[5+sp]
	if asy eq 0 then print, 'a single value of the asymmetry has been fitted for the pair'
	if asy eq 1 then print, 'the asymmetry has been fitted seperately for each mode'
	if sp eq 0 then print, 'the rotational splitting has been fixed at 0.4!4l!3Hz'
	if sp eq 1 then print, 'the rotational splitting has been fitted for each mode'
	if mh eq 0 then print, 'the height of the inner components of the l=3 mode were fixed with a fractional height of', ich
	if mh eq 1 then print, 'the height of the inner components of the l=3 mode were fitted'

	;l=1 modes
	if sp eq 0ull then begin 
		freq_neg1=A[2-asy]-0.4
		freq_pos1=A[2-asy]+0.4
	endif
	if sp eq 1ull then begin 
		freq_neg1=A[2-asy]-A[3-asy]
		 freq_pos1=A[2-asy]+A[3-asy]
	endif
	part_neg1=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_neg1)/exp(A[3+sp-asy])
	part_pos1=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_pos1)/exp(A[3+sp-asy])
	if asy eq 0 then begin 
		model_l1_neg1=(exp(A[4+sp-asy])/(1+part_neg1^2))*((1+A[1]*part_neg1)^2+A[1]^2)
		model_l1_pos1=(exp(A[4+sp-asy])/(1+part_pos1^2))*((1+A[1]*part_pos1)^2+A[1]^2)
	endif
	if asy eq 1 then begin 
		model_l1_neg1=(exp(A[4+sp-asy])/(1+part_neg1^2))*((1+A[4+sp]*part_neg1)^2+A[4+sp]^2)
		 model_l1_pos1=(exp(A[4+sp-asy])/(1+part_pos1^2))*((1+A[4+sp]*part_pos1)^2+A[4+sp]^2)
	endif

	;l=3 modes
	if sp eq 0ull then begin 
		freq_3neg1=A[5+sp]-0.4
		freq_3pos1=A[5+sp]+0.4
		freq_3neg3=A[5+sp]-1.2
		freq_3pos3=A[5+sp]+1.2
	endif
	if sp eq 1ull then begin 
		freq_3neg1=A[5+sp]-A[6+sp]
		freq_3pos1=A[5+sp]+A[6+sp]
		freq_3neg3=A[5+sp]-3*A[6+sp]
		freq_3pos3=A[5+sp]+3*A[6+sp]
	endif
	part_3neg1=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_3neg1)/exp(A[6+2*sp])
	part_3pos1=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_3pos1)/exp(A[6+2*sp])
	part_3neg3=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_3neg3)/exp(A[6+2*sp])
	part_3pos3=2*(nu[NNlow2:NNlow2+(N2-1)]-freq_3pos3)/exp(A[6+2*sp])
	if asy eq 0 then begin 
		if  mh eq 0 then begin 
			model_l3_neg1=(ich*exp(A[7+2*sp])/(1+part_3neg1^2))*((1+A[1]*part_3neg1)^2+A[1]^2)
			model_l3_pos1=(ich*exp(A[7+2*sp])/(1+part_3pos1^2))*((1+A[1]*part_3pos1)^2+A[1]^2)
		endif
		if mh eq 1 then begin
			model_l3_neg1=(A[8+2*sp+asy]*exp(A[7+2*sp])/(1+part_3neg1^2))*((1+A[1]*part_3neg3)^2+A[1]^2)
			model_l3_pos1=(A[8+2*sp+asy]*exp(A[7+2*sp])/(1+part_3pos1^2))*((1+A[1]*part_3pos3)^2+A[1]^2)
		endif
		model_l3_neg3=(exp(A[7+2*sp])/(1+part_3neg3^2))*((1+A[1]*part_3neg3)^2+A[1]^2)
		model_l3_pos3=(exp(A[7+2*sp])/(1+part_3pos3^2))*((1+A[1]*part_3pos3)^2+A[1]^2)
	endif
	
	if asy eq 1 then begin 
		if mh eq 0 then begin
			model_l3_neg1=(ich*exp(A[7+2*sp])/(1+part_3neg1^2))*((1+A[8+2*sp]*part_3neg1)^2+A[8+2*sp]^2)
			model_l3_pos1=(ich*exp(A[7+2*sp])/(1+part_3pos1^2))*((1+A[8+2*sp]*part_3pos1)^2+A[8+2*sp]^2)
		endif
		if mh eq 1 then begin
			model_l3_neg1=(A[8+2*sp+asy]*exp(A[7+2*sp])/(1+part_3neg1^2))*((1+A[8+2*sp]*part_3neg1)^2+A[8+2*sp]^2)
			model_l3_pos1=(A[8+2*sp+asy]*exp(A[7+2*sp])/(1+part_3pos1^2))*((1+A[8+2*sp]*part_3pos1)^2+A[8+2*sp]^2)
		endif
		model_l3_neg3=(exp(A[7+2*sp])/(1+part_3neg3^2))*((1+A[8+2*sp]*part_3neg3)^2+A[8+2*sp]^2)
		model_l3_pos3=(exp(A[7+2*sp])/(1+part_3pos3^2))*((1+A[8+2*sp]*part_3pos3)^2+A[8+2*sp]^2)
	endif

;sum l=1 and l=3 components
model_plot=model_l1_neg1+model_l1_pos1+model_l3_neg1+model_l3_pos1+model_l3_neg3+model_l3_pos3+exp(A[0])

endif
print,''
print,'---------------------------------------------------------------------'

case platforma of
	0: oplot, nu[NNlow2:NNlow2+(N2-1)], model_plot, color=200
	1: oplot, nu[NNlow2:NNlow2+(N2-1)], model_plot, color=245
endcase
frekvencije=nu[NNlow2:NNlow2+(N2-1)]
snaga=pow2
profili=model_plot
lbroj=fix(lbrojac)
nbroj=fix(nbrojac)
save,frekvencije,snaga,profili,filename='/home/andic/podaci/profili_l'+strtrim(lbroj,1)+'n'+strtrim(nbroj,1)+'.dat'
print,'saved file:'+'/home/andic/podaci/profili_l'+strtrim(lbroj,1)+'n'+strtrim(nbroj,1)+'.dat'

end