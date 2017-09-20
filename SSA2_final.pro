;IDL program for Radiative Processes project 1
;Author: Daniel Herman
;Start date: 29-08-17


function partitions, temp
;calculates partition function for Schadee element (E)
;input: temperature (K)
;output: values for partition function in 4 array elements u1, u2, u3, u4

u = fltarr(4)
chiion=[7,16,31,51]	;ionization energies for Schadee element
k=8.61734D-5		;Boltzmann constant in eV/K

for r=0,3 do $
	for s=0, chiion[r]-1 do u[r]=u[r] + exp(-s/(k*temp))

return, u

end

;----------------------------------------------------------------------

function boltz, temp, r, s

;computes the Boltzmann population for level r,s of E
;inputs: temperature (k)
;	 r (ionization stage)
;	 s (level, starting at s=1)
;output: normalized population  n_(r,s)/N_r

u=partitions(temp)
kev=8.61734D-5 ;Bolztmann constant in ev/deg
nmbrs = 1./u[r-1]*exp(-(s-1)/(kev*temp))
return, nmbrs

end

;----------------------------------------------------------------------

function saha,temp,elpress,ionstage

; computes saha population fraction N_r/N of E
;input: temperature, electron pressure, ion stage

;constants:

kerg=1.380658D-16	;Boltzmann constant in erg*K
kev=8.61734D-5		;Boltzmann constant in eV/K
h=6.62607D-27		;Planck constant in erg*s
me=9.109390D-28		;electron mass in grams

;other values:

kevt=kev*temp
kergt=kerg*temp
eldens=elpress/kergt

chiion=[7,16,31,51]	;ionization energies for E
u=partitions(temp)	;store partition function values
u=[u,2]
sahaconst=(2*!pi*me*kergt/(h*h))^1.5*2./eldens
nstage=dblarr(5)	;double precision fltarr
nstage[0]=1.	
for r=0,3 do $
	nstage[r+1] = nstage[r]*sahaconst*u[r+1]/u[r]*exp(-chiion[r]/kevt)
ntotal=total(nstage)	;sum all stages = element density
nstagerel=nstage/ntotal	;fractions of the element density
return,nstagerel[ionstage-1]
end


;--------------------------------------------------------------------------


function sahabolt,temp,elpress,ion,level

;compute Saha-Boltzmann population n_(r,s)/N for level r,s of E
;input: temperature, electron pressure, ionization stage, level nr

return, saha(temp,elpress,ion) * boltz(temp,ion,level)

end

;-------------------------------------------------------------------------


function sahabolt_H,temp,elpress,level

; computes the Saha-Boltzmann populations (n_(1,s)/N_H) for Hydrogen level

; constants 

kerg=1.380658D-16		;Boltzmann constant in erg*K
kev=8.61734D-5			;Boltzmann constant in eV/K
h=6.62607D-27			;Planck constant in erg*s
me=9.109390D-28			;electron mass in grams

;other values:

kevt=kev*temp
kergt=kerg*temp
eldens=elpress/kergt

; energy levels and weights for hydrogen
nrlevels = 100		
g = intarr(2,nrlevels)		; declaration of weights
chiexc = fltarr(2,nrlevels)	; declaration of excitation energies
for s=0, nrlevels-1 do begin
	g[0,s]=2*(s+1)^2
	chiexc[0,s]=13.598*(1-1./(s+1)^2)
endfor
g[1,0]=1
chiexc[1,0]=0.

; partition functions
u = fltarr(2)
u[0] = 0
for s=0, nrlevels-1 do u[0]=u[0] + g[0,s]*exp(-chiexc[0,s]/kevt)
u[1]=g[1,0]

; Saha values

sahaconst=(2*!pi*me*kergt/(h*h))^1.5*2./eldens
nstage=dblarr(5)		;double precision fltarr
nstage[0]=1.
nstage[1]=nstage[0]*sahaconst*u[1]/u[0]*exp(-13.598/kevt)
ntotal = total(nstage)		; sum the two stages = total hydrogen density

; Boltzmann	
nlevel = nstage[0]*g[0,level-1]/u[0]*exp(-chiexc[0,level-1]/kevt)
nlevelrel=nlevel/ntotal		; fraction of total hydrogen density

;STOP

return,nlevelrel
end

;---------------------- Recreation for Calcium ------------------------

function partitions_ca, temp
;calculates partition function for Calcium
;input: temperature (K)
;output: values for partition function in 4 array elements u1, u2, u3, u4

u = fltarr(4)
chiion=[6.113,11.871,50.91,67.15]	;ionization energies for Ca
k=8.61734D-5		;Boltzmann constant in eV/K

for r=0,3 do $
	for s=0, chiion[r]-1 do u[r]=u[r] + exp(-s/(k*temp))

return, u

end

;----------------------------------------------------------------------

function boltz_ca, temp, r, s

;computes the Boltzmann population for level r,s of Ca
;inputs: temperature (k)
;	 r (ionization stage)
;	 s (level, starting at s=1)
;output: normalized population  n_(r,s)/N_r

u=partitions_ca(temp)
kev=8.61734D-5 ;Bolztmann constant in ev/deg
nmbrs = 1./u[r-1]*exp(-(s-1)/(kev*temp))
return, nmbrs

end

;----------------------------------------------------------------------

function saha_ca,temp,elpress,ionstage

; computes saha population fraction N_r/N of Ca
;input: temperature, electron pressure, ion stage

;constants:

kerg=1.380658D-16	;Boltzmann constant in erg*K
kev=8.61734D-5		;Boltzmann constant in eV/K
h=6.62607D-27		;Planck constant in erg*s
me=9.109390D-28		;electron mass in grams

;other values:

kevt=kev*temp
kergt=kerg*temp
eldens=elpress/kergt

chiion=[6.113,11.871,50.91,67.15]	;ionization energies for Ca
u=partitions_ca(temp)	;store partition function values
u=[u,2]
sahaconst=(2*!pi*me*kergt/(h*h))^1.5*2./eldens
nstage=dblarr(5)	;double precision fltarr
nstage[0]=1.	
for r=0,3 do $
	nstage[r+1] = nstage[r]*sahaconst*u[r+1]/u[r]*exp(-chiion[r]/kevt)
ntotal=total(nstage)	;sum all stages = element density
nstagerel=nstage/ntotal	;fractions of the element density
return,nstagerel[ionstage-1]
end


;--------------------------------------------------------------------------


function sahabolt_ca,temp,elpress,ion,level

;compute Saha-Boltzmann population n_(r,s)/N for level r,s of E
;input: temperature, electron pressure, ionization stage, level nr

return, saha_ca(temp,elpress,ion) * boltz_ca(temp,ion,level)

end

;-------------------------------------------------------------------------


;temp=1000*indgen(31)
;print, temp
;pop=fltarr(5,31)
;for T=1,30 do $
;	for r=1,4 do pop[r,T]=sahabolt(temp[T],131.,r,1)
;plot,temp,pop[1,*],/ylog,yrange=[1E-3,1.1], $
;	xtitle='Temperature',ytitle='Population'
;oplot,temp,pop[2,*]
;oplot,temp,pop[3,*]
;oplot,temp,pop[4,*]

;STOP

;for T=1,30 do $
;        for r=1,4 do pop[r,T]=sahabolt(temp[T],131.,r,2)
;plot,temp,pop[1,*],/ylog,yrange=[1E-3,1.1], $
;	xtitle='Temperature',ytitle='Population'
;oplot,temp,pop[2,*]
;oplot,temp,pop[3,*]
;oplot,temp,pop[4,*]

;STOP

;for T=1,30 do $
;        for r=1,4 do pop[r,T]=sahabolt(temp[T],131.,r,4)
;plot,temp,pop[1,*],/ylog,yrange=[1E-3,1.1], $
;	xtitle='Temperature',ytitle='Population'
;oplot,temp,pop[2,*]
;oplot,temp,pop[3,*]
;oplot,temp,pop[4,*]


;end


;-------------------------------------------------------------------------

;temp=indgen(191)*100.+1000.       ; T = 1000-20000 in delta T = 100
;CaH = temp
;Caabund=2.E-6                     ; A_Ca = N_Ca / N_H
;for i=0,190 do begin
;NCa = sahabolt_ca(temp[i],1e2,2,1)
;NH = sahabolt_H(temp[i],1e2,2)
;CaH[i]=NCa*Caabund/NH
;endfor
;plot,temp,CaH,/ylog,$
;xtitle='temperature',ytitle="Ca II K / H alpha'

;end

;--------------------------------------------------------

;  graphing the temperature sensitivity for CaIIK and Halpha
;temp = indgen(101)*100.+2000.
;dNCadT = temp
;dNHdT = temp
;dT = 1.
;for i=0,100 do begin
;NCa = sahabolt_ca(temp[i],1e2,2,1)
;NCa2 = sahabolt_ca(temp[i]-dT,1e2,2,1)
;dNCadT[i] = (NCa - NCa2)/dT/NCa
;NH = sahabolt_H(temp[i],1e2,2)
;NH2 = sahabolt_H(temp[i]-dT,1e2,2)
;dNHdT[i] = (NH-NH2)/dT/NH
;endfor
;plot,temp,abs(dNHdT),/ylog,yrange=[1e-5,1],$
;xtitle='temperature',ytitle='abs d n(r,s) / n(r,s)'
;oplot,temp,abs(dNCadT),linestyle=2

;STOP

; recompute the arrays and overplot relative populations

;NCa = temp
;NH = temp
;for i=0,100 do begin
;NCa[i] = sahabolt_ca(temp[i],1e2,2,1)
;NH[i] = sahabolt_H(temp[i],1e2,2)
;endfor
;oplot,temp,NH/max(NH)
;oplot,temp,NCa/max(NCa),linestyle=2
;-------------------------------------------------------------------------

;temp = indgen(191)*100.+1000.
;nH=temp
;for i=0,190 do nH[i]=sahabolt_H(temp[i],1e2,1)
;plot,temp,nH,$
;xtitle='temperature',ytitle='neutral hydrogen fraction'


end

