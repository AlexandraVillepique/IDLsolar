function fshft_qub, line, sh, nst, nend, filt=filt, run=run, linear=linear
;+
;
;	function:  fshft
;
;	purpose:  shift array line by non-integer pixel shift sh by fourier
;		  or linear interpolation; uses wraparound for ends
;
;	notes:  1/95, lites@ncar (& rob@ncar) -
;			"extend"ed array to 256 points;
;			added frequency filtering option;
;			removed extraneous code and variables.
;		10/96, lites, put in keyword for possibility of linear
;			interpolation.
;
;==============================================================================
;
;	Check number of parameters.
;
if n_params() ne 4 then begin
	print
	print, "usage:  ret = fshft(line, sh, nst, nend)"
	print
	print, "	Shift array line by non-integer pixel shift sh by"
	print, "	linear interpolation.  Uses wraparound for ends."
	print, "	The parameter sh would be the negative of the result"
	print, "	of corshft to shift line2 back onto line1."
	print
	print, "	Arguments"
	print, "		line	- one-dimensional array (vector)"
	print, "		sh	- fractional pixel shift"
	print, "		nst	- starting index for shift"
	print, "		nend	- ending index for shift"
	print
	print, "	Keywords"
	print, "		filt	- if set, filter out fringes"
	print, "			  (hardwired range; def=no filtering)"
	print, "		run	- date of run for run-specific"
	print, "			  processing (def='may94')"
	print, "		linear  - if set, use linear instead of "
	print, "			  fourier interpolation   "
	print
	return, 0
endif
;-
;
;	Set general parameters.
;
nx = n_elements(line)
nx1 = nx - 1
if n_elements(run) eq 0 then run = 'may94'
;
;	Always interpolate to 1024 points.
;
nmask = 2048 - (nend - nst + 1)
if nmask le 14 then $
	message, 'no. pts ' + stringit(nmask) + ' is too close to 1024 limit.'
;
;	Extend the array so that there is smooth transition in wraparound.
;
linesh = fltarr(2048)
linesh(nst:nend) = line(nst:nend)
linesh = extend(linesh, nst, nend)
;
;	Shift the array by the *integer* pixel.
;
nsh = fix(sh)
del = sh
;
if abs(sh) ge 1.0 then begin
	linesh = shift(linesh, nsh)
	del = sh - nsh
endif
;
;	Shift the array by the *non-integer* remainder.
;
if (del ne 0.0) then begin
	if keyword_set(linear) then begin
;  optional linear interpolation
		if del gt 0. then dss = 1 else dss = -1
		temp1 = shift(linesh,dss)
		if del gt 0. then begin
			linesh = linesh*(1-del) + temp1*del
		endif else begin
			dds = abs(del)
			linesh = linesh*(1.-dds) + temp1*dds
		endelse

	endif else begin
;  fourier interpolation
		linesh = ffterpol(linesh, del, nst, nend)
	endelse
endif


;
;	Optionally filter out fringes.
;
if keyword_set(filt) then begin
	if comp_runs(run, 'lt', 'may94') then begin
		linesh = filtfreq(linesh, 81, 87)
	endif else begin
		linesh = filtfreq(linesh, 76, 77)
		linesh = filtfreq(linesh, 85, 86)
	endelse
endif

return, linesh(0:nx1)
end
