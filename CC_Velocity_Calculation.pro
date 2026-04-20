; ---------------------------------------------------------------------------
; CROSS-CORRELATION VELOCITY CALCULATIONS
;
; This procedure estimates propagation speeds by measuring time delays
; between correlated signals along a time-distance map.
;
; INPUTS:
;   xx          - 2D time-distance data array [time, distance]
;   selected    - Indices of valid spatial positions for correlation
;   period      - Period of oscillation [minutes]
;   cad_s       - Cadence [seconds]
;   pixel_size  - Pixel size [Mm]
;
; OUTPUTS:
;   velkmps     - Computed velocity [km/s]
;   vel_error   - Uncertainty in velocity [km/s]
;   length      - Distances along slit [Mm]
;   times       - Time delays for max cross-correlation [s]
;   gradient    - Fitted distance-time curve [Mm]
;
; OPTIONAL KEYWORDS:
;   VISUALISE   - If numeric, displays the plot in new window.
;   SAVENAME    - If string, saves the plot as a JPEG.
;
; ---------------------------------------------------------------------------

pro cc_vel, xx, selected, period, cad_s, pixel_size, velkmps, $
  vel_error, length, times, gradient, max_height, FILTER=filter, VIEW_OG=view_og, $
  VISUALISE=visualise, SAVENAME=savename

  maxx=[]  ; Array to store the time delays of maximum cross-correlation
  sz=size(xx)
  xlen=round(max_height/pixel_size)
;  xlen=sz[2]   ; Number of spatial positions
  tlen=sz[1]   ; Number of time steps
  pixels=indgen(xlen)
  tt=indgen(tlen)
  frame_P=period*60/cad_s  ; Convert period from minutes to frame units
  
  xx_f=xx
  if keyword_set(filter) then begin
    for i =0, n_elements(xx[0, *]) -1 do begin
      to_filt=xx[*,i]
      ff=1./(period*60.)
      ff_width=ff/5.
      band=[ff-ff_width/2., ff+ff_width/2.]

      filtsig, tt*cad_s, to_filt, band, xf, sd, filt='ideal'
      xx_f[*,i]=xf
    endfor
  endif

  
  
  ; Loop through spatial positions and compute cross-correlation
  for i = 0, xlen-1 do begin  
    cc = c_correlate(xx[*,0], xx[*,i],tt)   ; Cross-correlate with reference column
    minv = min(cc[0:frame_P], fmin)         ; Find first minimum within a period window
    cccrop=cc[fmin:fmin+frame_P]            ; Crop correlation function around this min
    maxv = max(cccrop, ff)                  ; Find maximum cross-correlation
    fmax=ff+fmin                            ; Convert index to absolute time shift
    
    ; Handle cases where the detected max jumps backwards too far
    if i ne 0 then begin
      if fmax lt (maxx[i-1]-50/cad_s) then begin
        cccrop=cc[fmax:fmax+frame_P]
        minv=min(cccrop, ff)
        fmin=fmax+ff
        cccrop=cc[fmin:fmin+frame_P]
        maxv = max(cccrop, ff)
        fmax=ff+fmin
      endif
    endif
    
    maxx=[maxx, fmax]  ; Store the detected max shift
    
  endfor
  
  ; Select only the required indices
  maxx=maxx[selected]
  times=maxx*cad_s                          ; Convert to seconds
  
  pixels=pixels[selected]
  length=pixels*pixel_size                  ; Convert to distance in Mm
  
  ;Flip fit other way around
  fit_coeffs = poly_fit(times, length, 1, sigma=sigma)
;  print, fit_coeffs
  grad = fit_coeffs[1]
  grad_error = sigma[1]
  t_10min = indgen(10)
  t_20min = t_10min+10
  t_30min = t_10min+20
  t_40min = t_10min+30
  t_50min = t_10min+40
  t_60min = t_10min+50
  gradient = fit_coeffs[0]+t_10min*60*grad
  velkmps = grad*1000
  vel_error = grad_error*1000
  
  
;  ; Fit a linear model to extract the velocity
;  fit_coeffs = poly_fit(length, times, 1, sigma=sigma)
;  grad=fit_coeffs[1]                        ; Gradient of the linear fit
;  grad_error=sigma[1]                       ; Error in gradient
;  vel_errorMmps=grad_error/(grad^2)
;  vel_error=vel_errorMmps*1000*2            ; Convert error to km/s
;  gradient = poly(tt, fit_coeffs)           ; Generate fitted line for plotting
;  velMmps = 1./grad                         ; Velocity in Mm/s
;  velkmps = velMmps*1000                    ; Convert to km/s

    ; Plot results if requested
 if keyword_set(visualise) then begin
    if visualise eq 'in_cur_win' then begin
      ;Use current window
      if keyword_set(filter) then loadct, 39, /silent else aia_lct, wave=171, /load
      if keyword_set(view_og) then xx_f=xx 
      if keyword_set(view_og) then aia_lct, wave=171, /load
      plot_image, xx_f, title='Cross-Correlation Speed', $
        scale=[cad_s/60., pixel_size], xtitle='Time [min]', ytitle='Distance [Mm]', /nosquare  
      loadct, 39, /silent
      if keyword_set(filter) then begin
        oplot, times/60., length, color=255, psym=2, symsize=2, thick=1
        oplot, t_10min, gradient, color=0, linestyle=2, thick=2
        oplot, t_20min, gradient, color=0, linestyle=2, thick=2
        oplot, t_30min, gradient, color=0, linestyle=2, thick=2
        oplot, t_40min, gradient, color=0, linestyle=2, thick=2
        oplot, t_50min, gradient, color=0, linestyle=2, thick=2
        oplot, t_60min, gradient, color=0, linestyle=2, thick=2
      endif else begin
        oplot, times/60., length, color=100, psym=2, symsize=1, thick=1
        oplot, t_10min, gradient, color=100, linestyle=2, thick=2
        oplot, t_20min, gradient, color=100, linestyle=2, thick=2
        oplot, t_30min, gradient, color=100, linestyle=2, thick=2
        oplot, t_40min, gradient, color=100, linestyle=2, thick=2
        oplot, t_50min, gradient, color=100, linestyle=2, thick=2
        oplot, t_60min, gradient, color=100, linestyle=2, thick=2
      endelse
      
    endif else begin
      ; Open specified window number
      wdef, visualise, 700, 500
      !p.multi = 0
      aia_lct, wave=171, /load
      plot_image, xx[0:10*60/float(cad_s),*], title='Cross-Correlation Speed', $
        scale=[cad_s, pixel_size], xtitle='Time [min]', ytitle='Distance [Mm]', /nosquare
      loadct, 39, /silent
      oplot, times, length, color=100, psym=2, symsize=2, thick=2
      oplot, t_10min, gradient, color=100, linestyle=2, thick=2
    endelse
  endif

  ; Output result
;  print, 'Velocity = ' + strtrim(velkmps, 2) + ' ± ' + strtrim(vel_error, 2) + ' km/s'

  ; Save JPEG if specified
  if keyword_set(savename) then $
    write_jpeg, strtrim(savename,2)+'.jpg', tvrd(/true), /true, quality=100

end
