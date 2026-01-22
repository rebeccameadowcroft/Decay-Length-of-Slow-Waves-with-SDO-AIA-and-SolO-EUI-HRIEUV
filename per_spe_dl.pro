;Per_spe_dl procedure
;Creates a plot that has determined wave period, speed and decay length from TD map for AIA and EUI data.
;

function propagate_division, x, dx, y, dy
  ; Check for zero denominator
  if y eq 0 then begin
    message, 'Error: denominator y is zero.'
  endif
  ; Calculate z
  z = x / y
  ; Propagate error
  dz = abs(z) * sqrt( (dx/x)^2 + (dy/y)^2 )
  return, [z, dz]
end

function round_num, value, sigfigs
  v_string=strtrim(round(value*10^(sigfigs-1))/10.^(sigfigs-1),2)
  v_string=v_string.substring(0, sigfigs)
  return, v_string
end

pro per_spe_dl, aia_td, aia_index, aia_pixel_size, aia_cad, $
                eui_td, eui_index, eui_pixel_size, eui_cad, $
                winnum, SAVENAME=savename, crop_option=crop_option, dl_region=dl_region
  
;  wdef, winnum, 700, 700
;  !p.multi=[0, 2, 4]
;  !p.charsize=1.7
;  !p.background=255
;  !p.color=0
;  !p.charthick=1
;  !p.thick=1
;  
  aia_og=aia_td
  eui_og=eui_td
;  
;  aia_lct, wave=171, /load
;  plot_image, aia_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='AIA TD Map '+strtrim(savename, 2), /nosquare, $
;    scale=[aia_cad, aia_pixel_size];, max=750
;    
;  plot_image, eui_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='EUI TD Map '+strtrim(savename, 2), /nosquare, $
;    scale=[eui_cad, eui_pixel_size];, max=750
  
  ;Detrend td map
  
  aia_td=aia_td-smooth(aia_td, [250./aia_cad,1], /edge_mirror)
  eui_td=eui_td-smooth(eui_td, [250./eui_cad,1], /edge_mirror)
  
;  plot_image, aia_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='Detrended AIA TD Map', /nosquare, $
;    scale=[aia_cad, aia_pixel_size]
;
;  plot_image, eui_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='Detrended EUI TD Map', /nosquare, $
;    scale=[eui_cad, eui_pixel_size]
  
  
  ;interpolate in time
  nT_a = n_elements(aia_td[*,0])
  tA = findgen(nT_a) * aia_cad          ; seconds
  nT_e = n_elements(eui_td[*,0])
  tE = findgen(nT_e) * eui_cad
  tmax = min([max(tA), max(tE)])
  dt_new = 1.0                          ; seconds (choose 1 or 5 etc.)
  nTnew = floor(tmax/dt_new) + 1
  tNew = findgen(nTnew) * dt_new
  aia_inter = fltarr(nTnew, n_elements(aia_td[0,*]))
  eui_inter = fltarr(nTnew, n_elements(eui_td[0,*]))
  for h=0, n_elements(aia_td[0,*])-1 do aia_inter[*,h] = interpol(aia_td[*,h], tA, tNew)
  for h=0, n_elements(eui_td[0,*])-1 do eui_inter[*,h] = interpol(eui_td[*,h], tE, tNew)

  
  ;interpolate in height
  nHa = n_elements(aia_td[0,*])
  hA = findgen(nHa) * aia_pixel_size
  nHe = n_elements(eui_td[0,*])
  hE = findgen(nHe) * eui_pixel_size
  Hmax = min([max(hA), max(hE)])
  dh = 0.1
  nHnew = floor(Hmax/dh) + 1
  hNew = findgen(nHnew) * dh
  aia_hires = fltarr(nTnew, nHnew)
  eui_hires = fltarr(nTnew, nHnew)
  for t=0, nTnew-1 do begin
    aia_hires[t,*] = interpol(aia_inter[t,*], hA, hNew)
    eui_hires[t,*] = interpol(eui_inter[t,*], hE, hNew)
  endfor
  
;  print, 'aia_hires: '
;  print, size(aia_hires)
;  print, 'eui_hires: '
;  print, size(eui_hires)
  
  ;cross correlate
  dists = []
  a_xs =[]
  xs=[]
  lag_n = 100
  plus_lag=indgen(lag_n)
  minus_lag=-reverse(plus_lag[1:-1])
  lag=[minus_lag, plus_lag]

  aia_max=(n_elements(aia_td[0,*])*aia_pixel_size)
  eui_max=(n_elements(eui_td[0,*])*eui_pixel_size)
  ;  print, 'aia_max='+strtrim(aia_max)+'  &   eui_max='+strtrim(eui_max)
  max_height_Mm=round(min([aia_max,eui_max]))
  ;  print, 'height [Mm] = ' +strtrim(max_height_Mm,2)

  ;Varying height for AIA/EUI
  eui_xs=[]
  aia_xs=[]
  for i = 0, nHnew-1 do begin
    ;Changing aia height
    aa=aia_hires[*,i]
    bb=eui_hires[*,0]
    aia_cc = c_correlate(aa, bb, lag)
    aia_maxx = max(aia_cc, index)
    aia_lag = lag[index]
;    print, 'aia_lag = '+strtrim(aia_lag)
    aia_xs=[aia_xs, aia_maxx]
    ;Changing eui height
    aa=aia_hires[*,0]
    bb=eui_hires[*,i]
    eui_cc = c_correlate(aa, bb, lag)
    eui_maxx = max(eui_cc, index)
    eui_lag=lag[index]
;    print, 'eui_lag = '+strtrim(eui_lag)
    eui_xs=[eui_xs, eui_maxx]
    dists = [dists, hNew[i]]
  endfor

;  plot, dists, aia_xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='Fix EUI, AIA change height w/o AC', yrange=[0,1]
;  plot, dists, eui_xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='Fix AIA, EUI change w/o AC', yrange=[0, 1]
  
 
  
  ;crop to TD map to height were peak is, print which is being cropped OR option to choose either AIA, EUI or neither
  if keyword_set(crop_option) then crop_option=crop_option else read,'Which TD map should be cropped? (1 - AIA, 2 - EUI, 3 - neither)',crop_option
  crop_option=strtrim(crop_option,2)
  
  ; --- decide crop ---
  if crop_option eq 1 then begin
    max_value = max(aia_xs, max_index)
    aia_crop_idx = max_index[0]
    eui_crop_idx = 0
    print, aia_crop_idx
  endif else if crop_option eq 2 then begin
    max_value = max(eui_xs, max_index)
    aia_crop_idx = 0
    eui_crop_idx = max_index[0]
    print, eui_crop_idx
  endif else begin
    aia_crop_idx = 0
    eui_crop_idx = 0
  endelse
  
  aia_align = aia_hires[*, aia_crop_idx:-1]
  eui_align = eui_hires[*, eui_crop_idx:-1]
  
  aia_td = aia_td[*,aia_crop_idx*0.1/aia_pixel_size:-1]
  eui_td = eui_td[*,eui_crop_idx*0.1/eui_pixel_size:-1]
  
  ; --- match heights (same number of height bins) ---
  new_height = min([n_elements(aia_align[0,*]), n_elements(eui_align[0,*])])
  aia_match = aia_align[*, 0:new_height-1]
  eui_match = eui_align[*, 0:new_height-1]
  H_keep = (new_height-1) * dh   ; Mm max height in matched arrays
  aia_keep_pix = fix(H_keep / aia_pixel_size)
  eui_keep_pix = fix(H_keep / eui_pixel_size)
  aia_td = aia_td[*, 0:aia_keep_pix]
  eui_td = eui_td[*, 0:eui_keep_pix]

  ;UP TO HERE OK********************************************************************************

  ;autocorrelate new cropped array
  aia_auto=aia_match
  eui_auto=eui_match
  for i = 0, n_elements(aia_auto[0,*])-1 do begin
    aia_auto[*,i]=a_correlate(aia_match[*,i], indgen(n_elements(aia_match[*,0])))
  endfor
  for i = 0, n_elements(eui_auto[0,*])-1 do begin
    eui_auto[*,i]=a_correlate(eui_match[*,i], indgen(n_elements(eui_match[*,0])))
  endfor

  
  ;cross correlate on matched (cropped) arrays at the same height
  dists = []
  xs = []
  a_xs = []
  
  for i = 0, n_elements(aia_match[0,*]) - 1 do begin
    ; CC with detrended signals
    cc = aia_match[*, i]
    dd = eui_match[*, i]
    xx = c_correlate(cc, dd, lag)
    xs = [xs, max(xx, index)]
  
    ; CC with autocorrelated signals
    aa = aia_auto[*, i]
    bb = eui_auto[*, i]
    a_xx = c_correlate(aa, bb, lag)
    a_xs = [a_xs, max(a_xx, index)]
  
    dists = [dists, i*dh]
  endfor


;  plot, dists, xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='AIA & EUI CC same height without AC', yrange=[0, 1]
;  plot, dists, a_xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='AIA & EUI CC same height with AC', yrange=[0,1]

  
;  
;  wdef, winnum+1, 800, 800
;  !p.multi=[0,2,6]
;  !p.charsize=2
  
  ;Plotting OG td maps
  aia_lct, wave=171, /load
  plot_image, aia_og, xtitle='Time [min]', ytitle='Distance [Mm]', title='Original AIA TD Map', /nosquare, $
    scale=[aia_cad/60., aia_pixel_size]
  loadct, 0, /silent
  maxxx=n_elements(aia_og[*,0])*aia_cad/60.
  oplot, [0, maxxx, maxxx, 0, 0], [0,0, aia_crop_idx*0.1, aia_crop_idx*0.1, 0], color=180, thick=1.5
  polyfill, [0, maxxx, maxxx, 0, 0], [0,0, aia_crop_idx*0.1, aia_crop_idx*0.1, 0], color=180, $
    /line_fill, orientation=45, spacing = 0.15, thick=1.5
  polyfill, [0, maxxx, maxxx, 0, 0], [0,0, aia_crop_idx*0.1, aia_crop_idx*0.1, 0], color=180, $
    /line_fill, orientation=-45, spacing = 0.15, thick=1.5
  aia_lct, wave=171, /load
  plot_image, eui_og, xtitle='Time [min]', ytitle='Distance [Mm]', title='Original EUI TD Map', /nosquare, $
    scale=[eui_cad/60., eui_pixel_size]
  loadct, 0, /silent
  oplot, [0, maxxx, maxxx, 0, 0], [0,0, eui_crop_idx*0.1, eui_crop_idx*0.1, 0], color=180, thick=1.5
  polyfill, [0, maxxx, maxxx, 0, 0], [0,0, eui_crop_idx*0.1, eui_crop_idx*0.1, 0], color=180, $
    /line_fill, orientation=45, spacing = 0.15, thick=1.5
  polyfill, [0, maxxx, maxxx, 0, 0], [0,0, eui_crop_idx*0.1, eui_crop_idx*0.1, 0], color=180, $
    /line_fill, orientation=-45, spacing = 0.15, thick=1.5
  
  ;Plotting cropped and detrended td maps
  aia_lct, wave=171, /load
  plot_image, aia_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='Cropped, Detrended AIA TD Map', /nosquare, $
    scale=[aia_cad/60., aia_pixel_size]
  plot_image, eui_td, xtitle='Time [s]', ytitle='Distance [Mm]', title='Cropped, Detrended EUI TD Map', /nosquare, $
    scale=[eui_cad/60., eui_pixel_size]
  
  
  plot, dists, xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='AIA & EUI CC same height without AC', yrange=[0, 1]
  plot, dists, a_xs, color=0, xtitle='Height [Mm]', ytitle='R', psym=2, title='AIA & EUI CC same height with AC', yrange=[0,1]
  
  
  ;Determining wave period
  fourier_td_cum, aia_td, aia_pixel_size, aia_cad, aia_period, aia_perr_p, aia_perr_m, 171, $
    AC=AC, spectrum, VISUALISE='in_cur_win', half=half

  
  fourier_td_cum, eui_td, eui_pixel_size, eui_cad, eui_period, eui_perr_p, eui_perr_m, 171, $
    AC=AC, spectrum, VISUALISE='in_cur_win', half=half
  
  
;  cc_vel, aia_td, [0:15], aia_period, aia_cad, aia_pixel_size, aia_speed, $
;    aia_speed_error, aia_length, aia_times, aia_gradient, VISUALISE='in_cur_win', /filter, /view_og
;  cc_vel, eui_td, [0:25], eui_period, eui_cad, eui_pixel_size, eui_speed, $
;    eui_speed_error, eui_length, eui_times, eui_gradient, VISUALISE='in_cur_win', /filter, /view_og
  
  read,'What maximum height should be used for speed calculation? (in Mm)', aia_max_height
  
  ;Determine wave speeds
  cc_vel, aia_td, [0:15], aia_period, aia_cad, aia_pixel_size, aia_speed, $
    aia_speed_error, aia_length, aia_times, aia_gradient, aia_max_height, VISUALISE='in_cur_win', /filter
  

  ;Choose region for decay length
  if keyword_set(dl_region) then begin
    start=dl_region[0] & finish=dl_region[1]
;    print, 'start = ',start
;    print, 'finish = ',finish 
  endif else begin
    read,'Starting time for decay length (min):',start
    read,'Starting time for decay length (min):',finish
  endelse
  max=n_elements(aia_td[0, *])*aia_pixel_size
  oplot, [start, finish, finish, start, start], [0, 0, max*0.95, max*0.95, 0], color=250, thick=2
  
  cc_vel, eui_td, [0:25], eui_period, eui_cad, eui_pixel_size, eui_speed, $
    eui_speed_error, eui_length, eui_times, eui_gradient, aia_max_height, VISUALISE='in_cur_win', /filter
  
  max=n_elements(eui_td[0, *])*eui_pixel_size
  oplot, [start, finish, finish, start, start], [0, 0, max*0.95, max*0.95, 0], color=250, thick=2
  
  
  ;Determine decay length
  aia_dl_td=aia_td[start*60./aia_cad:finish*60./aia_cad, *]
  expfit_decay_length, aia_dl_td, aia_period, aia_cad, aia_pixel_size, aia_max_index, aia_efold_len, $
    aia_efold_err_plus, aia_efold_err_minus, VISUALISE='in_cur_win'
  
  eui_dl_td=eui_td[start*60./eui_cad:finish*60./eui_cad, *]
  expfit_decay_length, eui_dl_td, eui_period, eui_cad, eui_pixel_size, eui_max_index, eui_efold_len, $
    eui_efold_err_plus, eui_efold_err_minus, VISUALISE='in_cur_win'

  
  ;adding labels to plot
  ydiff=0.167
  aia_xpos = 0.35
  eui_xpos = 0.85
  per_ypos = 0.455
  speed_ypos = per_ypos - ydiff
  dl_ypos = per_ypos - 2*ydiff
  xyouts, aia_xpos, per_ypos, round_num(aia_period, 2)+'!U + '+round_num(aia_perr_p, 2)+'!N min', /normal, charsize=1
  xyouts, aia_xpos, per_ypos, round_num(aia_period, 2)+'!D - '+round_num(aia_perr_m, 2), /normal, charsize=1
  xyouts, eui_xpos, per_ypos, round_num(eui_period, 2)+'!U + '+round_num(eui_perr_p, 2)+'!N min', /normal, charsize=1
  xyouts, eui_xpos, per_ypos, round_num(eui_period, 2)+'!D - '+round_num(eui_perr_m, 2), /normal, charsize=1
  xyouts, aia_xpos, speed_ypos, round_num(aia_speed, 2)+cgsymbol('+-')+round_num(aia_speed_error, 1)+' kms!U-1', /normal, charsize=1
  xyouts, eui_xpos, speed_ypos, round_num(eui_speed, 2)+cgsymbol('+-')+round_num(eui_speed_error, 1)+' kms!U-1', /normal, charsize=1
  xyouts, aia_xpos, dl_ypos, round_num(aia_efold_len, 2)+'!U + '+round_num(aia_efold_err_plus,2)+'!N Mm', /normal, charsize=1
  xyouts, aia_xpos, dl_ypos, round_num(aia_efold_len, 2)+'!D + '+round_num(aia_efold_err_minus,2), /normal, charsize=1
  xyouts, eui_xpos, dl_ypos, round_num(eui_efold_len, 2)+'!U + '+round_num(eui_efold_err_plus,2)+'!N Mm', /normal, charsize=1
  xyouts, eui_xpos, dl_ypos, round_num(eui_efold_len, 2)+'!D + '+round_num(eui_efold_err_minus,2), /normal, charsize=1

  

    
    
    
  ;Propagating erros for speed/dl errors
  c_props = propagate_division(AIA_speed, AIA_speed_error, EUI_speed, EUI_speed_error)
  c_rat = c_props[0] & c_rat_err = c_props[1]
  
  dl_props = propagate_division(aia_efold_len, mean([aia_efold_err_plus, aia_efold_err_minus]), eui_efold_len, mean([eui_efold_err_plus, eui_efold_err_minus]))
  dl_rat = dl_props[0] & dl_rat_err = dl_props[1]

  rat_props = propagate_division(c_rat, c_rat_err, dl_rat, dl_rat_err)
  rat_rat = rat_props[0] & rat_rat_err = rat_props[1]
  
  ;Create txt file that writes all of the determined values -----------------------------------------------------------------
  get_lun, unit
  openw, unit, strtrim(savename,2)+'.txt'             ;create txt file in cd
  printf, unit, strtrim(savename,2)+ ' WAVE PROPERTIES:'
  printf, unit, ' '
  printf, unit, ' AIA Period = '+ strtrim(aia_period,2)+ ' + '+ strtrim(aia_perr_p,2) + ' - '+ strtrim(aia_perr_m,2)+ ' min'
  printf, unit, ' EUI Period = '+ strtrim(eui_period,2)+ ' + '+ strtrim(eui_perr_p,2) + ' - '+ strtrim(eui_perr_m,2)+ ' min'
  printf, unit, ' '
  printf, unit, ' AIA Speed = '+ strtrim(AIA_speed,2)+ ' +- '+ strtrim(AIA_speed_error,2)+ ' km/s'
  printf, unit, ' EUI Speed = '+ strtrim(EUI_speed,2)+ ' +- '+ strtrim(EUI_speed_error,2)+ ' km/s'
  printf, unit, ' AIA/EUI = '+ strtrim(c_rat,2) + ' ± ' + strtrim(c_rat_err,2)
  printf, unit, ' '
  printf, unit, ' AIA e-fold Length = '+ strtrim(aia_efold_len,2)+' Mm'+ ' + '+ strtrim(aia_efold_err_plus,2)+ ' - ' + strtrim(aia_efold_err_minus,2) +' Mm'
  printf, unit, ' EUI e-fold Length = '+ strtrim(eui_efold_len,2)+' Mm'+ ' + '+ strtrim(eui_efold_err_plus,2)+ ' - ' + strtrim(eui_efold_err_minus,2) +' Mm'
  printf, unit, ' AIA/EUI = ' +  strtrim(dl_rat,2) + ' ± ' +  strtrim(dl_rat_err,2)
  printf, unit, ' '
;  printf, unit, ' AIA exp-fit Decay Length = '+ strtrim(aia_efold_len,2)+ ' ± '+ strtrim(aia_efold_err,2)+ ' Mm'
;  printf, unit, ' EUI exp-fit Decay Length = '+ strtrim(eui_efold_len,2)+ ' ± '+ strtrim(eui_efold_err,2)+ ' Mm'
;  printf, unit, ' AIA/EUI = ' +  strtrim(aia_decay_len/eui_decay_len,2)
;  printf, unit, ' '
;  printf, unit, ' AIA Detection Length = '+ strtrim(aia_det_len,2)+' Mm';+ ' ± '+ strtrim(aia_det_err,2)+ ' Mm'
;  printf, unit, ' EUI Detection Length = '+ strtrim(eui_det_len,2)+' Mm';+ ' ± '+ strtrim(eui_det_err,2)+ ' Mm'
;  printf, unit, ' AIA/EUI = ' +  strtrim(aia_det_len/eui_det_len,2)
;  printf, unit, ' '
  printf, unit, 'Speed_ratio/e_fold_ratio = '+strtrim(rat_rat,2)+' ± '+strtrim(rat_rat_err,2)
;  printf, unit, 'Speed_ratio/e_fit_ratio = '+strtrim((AIA_speed/EUI_speed)/(aia_decay_len/eui_decay_len),2)
;  printf, unit, 'Speed_ratio/det_ratio = '+strtrim((AIA_speed/EUI_speed)/(aia_det_len/eui_det_len),2)
  

  close, unit

  
  
;  ;creating a csv file with the determined parameters
;  header=['AIA period', 'AIA perr p', 'AIA perr m', 'EUI_period', 'AIA perr p', 'AIA perr m', $
;    'AIA Speed', 'AIA Speed Error', 'EUI Speed', 'EUI Speed Error',$
;    'AIA e-fold Length', $
;    'AIA exp-fit Decay Length', 'AIA exp-fit err p', 'AIA exp-fit err m', 'EUI exp-fit Decay Length', 'EUI exp-fit err p', 'EUI exp-fit err m', 'AIA/EUI exp-fit'
;  
  
;  if keyword_set(savename) then write_jpeg, strtrim(savename,2)+'.jpg', TVRD(/TRUE), /TRUE, quality=100
  
  
  
  !p.multi=0
END


;;TRAIL RUN
;
;slitnumber=4
;
;;Create td maps
;aia_td=timedistance(aia_data_a, [x1s[2*slitnumber], y1s[2*slitnumber], x2s[2*slitnumber], y2s[2*slitnumber]], width=3)
;eui_td=timedistance(eui_data_a, [x1s[2*slitnumber+1], y1s[2*slitnumber+1], x2s[2*slitnumber+1], y2s[2*slitnumber+1]], width=5)
;
;
;per_spe_dl, aia_td, aia_index_a, aia_pixel_size, 12., $
;  eui_td, eui_index_a, eui_pixel_size, 3., $
;  4, SAVENAME='practice_run'
;  END