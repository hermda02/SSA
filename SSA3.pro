;----------------------Chapter 3 Code-----------------------

function planck,temp,wav

;Computes the planck function as a function of wavelength

; Defining constants
c = 2.99792D10 ;speed of light in cm/s
h = 6.62607D-27 ;Planck constant in ergs
k = 1.38065D-16 ;boltzmann constant in ergK-1

planck=(2*h*c^2)/((wav^5)*(exp((h*c)/(wav*k*temp))-1))

return,planck
end

;-------------------------Planck Graphs------------------------
;wav = indgen(100)*200.+1000.  ;makes a wavefunction array 1000-20800
;print, wav
;b = wav
;for i=0,99 do b[i]=planck(8000,wav[i]*1E-8)
;plot,wav,b,/ylog,xtitle='wavelength (Angstrom)',ytitle='Planck function',$
;xmargin=[15,15],$
;charsize = 1.2
;for T=8000,5000,-200 do begin 
;for i =0,99 do b[i]=planck(T,wav[i]*1E-8)
;oplot,wav,b
;endfor

;-------------------------Attenuation--------------------------------

;B=2.
;tau=indgen(101)/10.+0.01        ; tau array from 0.01 - 10
;int=tau
;for I0=4,0,-1 do begin          ; step down from I0=4 ... I0=0
;for i=0,100 do int[i]=I0*exp(-tau[i])+B*(1-exp(-tau[i]))
;if (i0 eq 4) then plot,tau,int,/xlog,/ylog,$
;xtitle='tau',ytitle='Intensity',charsize=1.3
;if (i0 ne 4) then oplot,tau,int
;endfor

;-------------------------Voigt Function----------------------------

;u=indgen(201)/10.-10.                   ; u = -10 to 10 in 0.1 steps
;vau=u
;a=0.1
;for i=0,200 do vau[i]=voigt(a,abs(u[i]))
;plot,u,vau,yrange=[0,1]

;-----------------------Schuster-Schwarzschild----------------------
;Ts=5700                              ; solar surface temp
;T1=4200                              ; solar T-min temperature 
;a=0.1
;wav=10000.D-8                         ; wavelength in cm
;tau0=1                               ; reversing layer thickness at line center
;u=indgen(201)/10.-10.
;int=u
;for i=0,200 do begin
;tau=tau0*voigt(a,abs(u[i]))
;int[i]=planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
;endfor
;plot,u,int

;------------------------------------------------------------------

;tau0=[0.01,0.05,0.1,0.5,1,5,10,50,100]
;for itau=0,8 do begin
;for i=0,200 do begin
;tau=tau0[itau]*voigt(a,abs(u[i]))
;int[i]=planck(Ts,wav)*exp(-tau) + planck(T1,wav)*(1.-exp(-tau))
;endfor
;oplot,u,int
;endfor

;-----------------------------------------------------------------

;Ts=5700                              ; solar surface temp
;T1=4200                              ; solar T-min temperature 
;a=0.1                        
;tau0=[0.01,0.05,0.1,0.5,1,5,10,50,100]
;u=indgen(201)/10.-10.
;int=u
;for iwav=1,3 do begin
;wav=(iwav^2+1)*1.D-5
;for itau=0,8 do begin
;for i=0,200 do begin
;tau=tau0[itau]*voigt(a,abs(u[i]))
;int[i]=planck(Ts,wav)*exp(-tau)+planck(T1,wav)*(1.-exp(-tau))
;endfor
;int=int/int[0]
;if (iwav eq 1 and itau eq 0) then plot,u,int
;if (iwav eq 1 and itau gt 0) then oplot,u,int
;if (iwav eq 2) then oplot,u,int,linestyle=1
;if (iwav eq 3) then oplot,u,int,linestyle=4
;endfor
;endfor

;---------------------------Profile Function-------------------------

function profile,a,tau0,u
; return a Schuster-Schwarzschild profile
; input: a = damping parameter
; tau0 = SS layer thickness at line center
; u = wavelength array in Doppler units
; output: int = intensity array
Ts=5700
T1=4200
wav=5000.E-8
int=u
usize=SIZE(u)
for i=0,usize[1]-1 do begin
tau=tau0*voigt(a,abs(u[i]))
int[i]=planck(Ts,wav)*exp(-tau)+planck(T1,wav)*(1.-exp(-tau))
endfor
return,int

end

;-------------------------------------------------------------------

u=indgen(1001)/2.5-200.
a=0.1
tau0=1e2
int=profile(a,tau0,u)
plot,u,int

STOP

;-------------------------------------------------------------------

reldepth=(int[0]-int)/int[0]
plot,u,reldepth
eqw=total(reldepth)*0.4
print,eqw

STOP

tau0=10^(indgen(61)/10.-2.)  ;10^-2 to 10^4, 0.1 steps in the log
eqw=tau0                     ; same size array
for i=0,60 do begin
int=profile(a,tau0[i],u)
reldepth=(int[0]-int)/int[0]
eqw[i]=total(reldepth)*0.4
endfor
plot,tau0,abs(eqw),xtitle='tau0',ytitle='abs(equivalent width)',/xlog,/ylog

end

