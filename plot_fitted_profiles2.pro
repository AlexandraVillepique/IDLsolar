pro plot_fitted_profiles2,nlow,nhigh

Common Share, nu_1, nu_2, Pow2, j, N1, N2, NNlow1, NNlow2, T, $
       PowWin_Short, sp, spv, asy, k, lm, mh, ich, sb, cv

set_plot,'x'
put='/home/andic/podaci/'
ime3='profili_l0n'
prez='.dat'
ime2='granice0.dat'
restore,put+ime2
uksnaga=snaga
ukfrekvencije=frekvencije
dim=size(snaga)
pom=dim(1)
profil_model=dblarr(pom)
loadct,3
i=fix(nlow)
poc=fix(nlow)
kra=fix(nhigh)
j=0.
for i=poc,kra do begin
ime=ime3+strtrim(i,1)+prez
restore,put+ime
pocetak=granice(0,0)
print,pocetak
raz1=granice(1,j)
print,raz1
po=granice(0,j)-pocetak
kraj=raz1-1.+po
print,po,kraj
profil_model(po:kraj)=profili
j=j+1.
endfor

plot,ukfrekvencije,uksnaga
oplot,ukfrekvencije,profil_model,color=150

ime3='profili_l1n'
loadct,1
ime2='granice1.dat'
restore,put+ime2
test1=max(frekvencije)
test2=min(frekvencije)
test3=max(ukfrekvencije)
test4=min(ukfrekvencije)
if test2 lt test4 then begin
	print,'druge frekvencije manji minimum'
endif
if test2 gt test4 then begin
	print,'prve frekvencije manji minumum'
	pomo=where(ukfrekvencije lt test2)
	dimi1=size(pomo)
	if test1 lt test3 then begin
		print,'druge frekvencije veci maksimum'
	endif
	if test1 gt test3 then begin
		pomo2=where(frekvencije gt test3)
		dimi2=size(pomo2)
		dimi3=size(ukfrekvencije)
		dimi=dimi2(1)*1.+dimi3(1)*1.
		ukfrek=dblarr(dimi)
		uksnag=ukfrek
		ukfrek(0:dimi3(1)-1)=ukfrekvencije
		ukfrek(dimi3(1):dimi3(1)+dimi2(1)-1)=frekvencije(pomo2)
		uksnag(0:dimi3(1)-1)=uksnaga
		uksnag(dimi3(1):dimi3(1)+dimi2(1)-1)=snaga(pomo2)
		prof_mod=dblarr(dimi)
		prof_mod(0:pom-1)=profil_model
	endif
endif
j=0.
for i=poc,kra do begin
ime=ime3+strtrim(i,1)+prez
restore,put+ime
print,pocetak
raz1=granice(1,j)
print,raz1
po=granice(0,j)-pocetak
kraj=raz1-1.+po
print,po,kraj
prof_mod(po:kraj)=profili
j=j+1.
endfor
window,0,xsize=1500,ysize=1000
loadct,1
plot,ukfrek,uksnag,title='All fitted range',ytitle='power',xtitle='Frequencies, micro Hz'
oplot,ukfrek,prof_mod,color=150
xyouts,3900,10000,'l=1/3',charsize=2,color=150 
loadct,3
oplot,ukfrek,profil_model,color=150
xyouts,3900,10700,'l=0/2',charsize=2,color=150

odgovor=''
read, odgovor, prompt='Do you wish to print the example?  (y/n) '
if odgovor eq 'n' then goto, jump5

str='example_fitting.ps'
set_plot,'ps'
     device,filename=str,encapsulated=1,/color,bits_per_pixel=24,xsize=20,ysize=15
loadct,1
plot,ukfrek,uksnag,title='All fitted range',ytitle='power',xtitle='Frequencies, micro Hz'
oplot,ukfrek,prof_mod,color=150
xyouts,3900,10000,'l=1/3',charsize=2,color=150 
loadct,3
oplot,ukfrek,profil_model,color=150
xyouts,3900,10700,'l=0/2',charsize=2,color=150
 device,/close
set_plot,'x'
jump5:

end