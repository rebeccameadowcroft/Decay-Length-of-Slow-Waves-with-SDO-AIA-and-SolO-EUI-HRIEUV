; EXPONENTIAL FIT DECAY LENGTH
;
; This procedure calculates the decay length of a periodic signal by:
; - filtering the data in time,
; - extracting amplitude profiles over distance,
; - fitting an exponential decrease from the max to end to identify the decay length,
; - estimating asymmetric error margins based on amplitude uncertainty.
;
; INPUT PARAMETERS:
; - data:       2D array of signal values [time, distance].
; - period:     Period of the oscillation in minutes.
; - cad_s:      Temporal cadence of the data in seconds.
; - pixel_size: Spatial size of a pixel in megametres (Mm).
;
; OUTPUT PARAMETERS:
; - decay_len:  Calculated decay length (in Mm).
; - dl_err:     Estimated uncertainty on the decay length (in Mm).
;
; OPTIONAL KEYWORDS:
; - VISUALISE:  Controls the plotting behavior:
;               * If set to 'in_cur_win', plots only the amplitude vs. distance panel
;                 in the current graphics window.
;               * Otherwise, opens a new window with three panels: original data, filtered
;                 data, and amplitude vs. distance with fit and error bars.
; - SAVENAME:   If set, saves the final plot as a JPEG with the given name.
; - PRINT:      If set, prints decay length and error to console.
;
pro expfit_decay_length, data, period, cad_s, pixel_size, max_index, fold_len, e_fold_len_plus, e_fold_len_minus, $
  det_len, $
  VISUALISE=visualise, SAVENAME=savename, PRINT=print

  cad_min = cad_s / 60.

  ; FILTER --------------------------------------------------------------
  ff = 1. / (period * 60.)
  fwidth = ff / 10.

  f_data = data
  sz = size(data)
  tlen = sz[1]
  tt = indgen(tlen) * cad_s

  heights = []
  freqs = []
  for i = 0, n_elements(data[0, *]) - 1 do begin
    xx = data[*, i]
    freqs = [freqs, ff]
    heights = [heights, i]
    filtsig, tt, xx, [ff - fwidth / 2, ff + fwidth / 2], xf, sd, filt='gauss'
    f_data[*, i] = xf
  endfor
  
;  s_data=smooth(data, [36./cad_s, 1])-smooth(data,[250./cad_s,1], /edge_truncate)

  ; Amplitude extraction and smoothing
;  plot_image, f_data, /nosquare
  ampslens, f_data, cad_s, meanamps, xlengths, allamps, alllens, mean_level, sig_level, pixel_size
  xlengths_Mm = xlengths * pixel_size
  alllens_Mm = alllens * pixel_size

  ; Error estimate on amplitudes
  meanamps_err = []
  for i = 0, alllens[-1] do begin
    stds = []
    for j = 0, n_elements(alllens) - 1 do begin
      if alllens[j] eq i then stds = [stds, allamps[j]]
    endfor
    if stds ne [] then meanamps_err = [meanamps_err, stddev(stds)]
  endfor

  ; Smooth + interpolate amplitudes
  width = round(2.0/pixel_size)
  smooth_meanamps = smooth(meanamps, width, /edge_truncate)
  smooth_ds = indgen(n_elements(smooth_meanamps) * 100) / 100.
  smooth_inter = interpolate(smooth_meanamps, smooth_ds)
  xlengths_inter = interpolate(xlengths_Mm, smooth_ds)

  ; Locate max amplitude
  A_max = max(smooth_inter)
  max_index = where(smooth_inter eq A_max[0])
  crop_index = round(max_index / 100.)

  ; Crop from max onward
  smooth_inter_crop = smooth_inter[max_index:-1]
  xlengths_inter_crop = xlengths_inter[max_index:-1]
  xlengths_Mm_crop = xlengths_Mm[crop_index:-1]
  meanamps_crop = meanamps[crop_index:-1]
  meanamps_err_crop = meanamps_err[crop_index:-1]

;  ; Fit exponential model
;  weights = 1.0 / (meanamps_err_crop)^2
;  result = COMFIT(xlengths_Mm_crop, meanamps_crop, [200, 0.8, 2], /exponential, weights=weights, sigma=sigmas)
;  a1_err = sigmas[1]
;  a1 = result[1]
;  decay_len = -1 / alog(a1) + xlengths_Mm[crop_index]
;  dl_err = sqrt(((-(1 / a1 / (alog(a1))^2)) * a1_err)^2)
  
;  fit = result[0] * result[1]^xlengths_inter_crop + result[2]
  
;  ; === Detection length with error ===
;  amp_err = mean(meanamps_err)   ; average error across all amplitudes
;  ; central value
;  AA = where(smooth_inter_crop lt sig_level)
;  det_len = xlengths_inter_crop[AA[0]]
;  
;  det_info_plus=smooth_inter_crop+amp_err
;  det_info_minus=smooth_inter_crop-amp_err
;  AA_plus=where(det_info_plus lt sig_level)
;  AA_minus=where(det_info_minus lt sig_level)
;  det_len_plus=xlengths_inter[AA_plus[0]]
;  det_len_minus=xlengths_inter[AA_minus[0]]
;  det_err = [abs(det_len-det_len_minus), abs(det_len_plus-det_len)]

  ; === E-folding length with error ===
  smooth_unmean = smooth_inter - mean_level
  Amax = max(smooth_unmean)

  ; central value
  e_height = Amax/2.71828
  BB = where(smooth_unmean lt e_height)
  e_fold_len = xlengths_inter[BB[0]]
  
  ; upper/lower by shifting Amax by amp_err
  amp_err = mean(meanamps_err)   ; average error across all amplitudes
  e_height_plus = (Amax+amp_err)/2.71828
  e_height_minus = (Amax-amp_err)/2.71828
  BB_plus = where(smooth_unmean lt e_height_plus)
  BB_minus = where(smooth_unmean lt e_height_minus)
  e_fold_len_plus = xlengths_inter[BB_minus[0]] - e_fold_len
  e_fold_len_minus = e_fold_len - xlengths_inter[BB_plus[0]]
  e_fold_err = [e_fold_len_minus, e_fold_len_plus]
  fold_len =  e_fold_len - xlengths_inter[max_index]

  
  
if keyword_set(print) then begin
  print, 'e-folding length = ' + strtrim(e_fold_len,2) + ' +'+strtrim(e_fold_err[1],2)+' / -'+strtrim(e_fold_err[0],2)+' Mm'
  ;print, 'Exp-fit Decay length = ' + strtrim(decay_len,2) + ' ± ' + strtrim(dl_err,2) + ' Mm'
;  print, 'Detection length = ' + strtrim(det_len,2) + ' +'+strtrim(det_err[1],2)+' / -'+strtrim(det_err[0],2)+' Mm'
endif

  
  
  ;PLOTTING -------------------------------------------------------------------------------------------------------------------  
  if keyword_set(visualise) then begin
    if visualise eq 'in_cur_win' then begin
      height = max(meanamps + meanamps_err) * 1.1
      plot, xlengths_Mm_crop, meanamps_crop, psym=2, title='Decay Length Determination', $
        xtitle='Distance [Mm]', ytitle='Amplitude [A.U.]', yrange=[0, height], xrange=[0, max(xlengths_Mm_crop)], /ysty, /xsty
      loadct, 0, /silent
      oplot, xlengths_Mm, meanamps, psym=2, color=150
;      loadct, 62, /silent
;      polyfill, [decay_len-dl_err+crop_index*pixel_size, decay_len+dl_err+pixel_size*crop_index, decay_len+dl_err+pixel_size*crop_index, decay_len-dl_err+pixel_size*crop_index, decay_len-dl_err+pixel_size*crop_index], $
;        [0, 0, height, height, 0], color=50
;      loadct, 0, /silent
      ; --- E-folding shaded box ---
      loadct, 50, /silent   ; choose a colour table index (adjust if needed)
      polyfill, [e_fold_len-e_fold_err[0], e_fold_len+e_fold_err[1], $
        e_fold_len+e_fold_err[1], e_fold_len-e_fold_err[0], $
        e_fold_len-e_fold_err[0]], $
        [0, 0, height, height, 0], color=50, /fill
      loadct, 0, /silent

      
;      ; --- Detection shaded box ---
;      if finite(det_err[0]) and finite(det_err[1]) then begin
;        loadct, 39, /silent  ; choose a different colour table index
;        polyfill, [det_len-det_err[0], det_len+det_err[1], $
;          det_len+det_err[1], det_len-det_err[0], $
;          det_len-det_err[0]], $
;          [0, 0, height, height, 0], color=120, /fill
;        loadct, 0, /silent
;      endif

      oplot, xlengths_Mm_crop, meanamps_crop, psym=2
      loadct, 39, /silent
      oplot, xlengths_inter, smooth_inter, color=250, thick=2
      errplot, xlengths_Mm_crop, meanamps_crop-meanamps_err_crop, meanamps_crop+meanamps_err_crop
;      oplot, [decay_len+pixel_size*crop_index, decay_len+pixel_size*crop_index], [0, height], linestyle=2, color=250
;      oplot, xlengths_inter_crop, fit, color=250, thick=1.5
;      oplot, [0, max(xlengths_Mm_crop)], [sig_level, sig_level], color=120, linestyle=2, thick=2
;      
      oplot, [0, max(xlengths_Mm_crop)], [mean_level, mean_level], color=50, linestyle=2
;      oplot, [0, max(xlengths_Mm_crop)], [e_height+mean_level, e_height+mean_level], color=50, linestyle=2
      oplot, [e_fold_len, e_fold_len], [0, height], color=50
      
;      oplot, [0, max(xlengths_Mm_crop)], [sig_level, sig_level], color=120, linestyle=2
;      oplot, [det_len, det_len], [0, height], color=120
    endif else begin
      ; CASE 2: Open new window with all 3 plots
      wdef, visualise, 600, 800
      !p.multi = [0,1,3]
      aia_lct, wave=171, /load
      plot_image, data, /nosquare, title='OG TD map', scale=[cad_s/60., pixel_size], xtitle='Time [min]', ytitle='Distance [Mm]'
      plot_image, f_data, /nosquare, title='Filtered TD map', scale=[cad_s/60., pixel_size], xtitle='Time [min]', ytitle='Distance [Mm]'
      loadct, 39, /silent
      height = max(meanamps + meanamps_err) * 1.1
      plot, xlengths_Mm_crop, meanamps_crop, psym=2, title='Decay Length by Exponential Fit', $
         xtitle='Distance [Mm]', ytitle='Amplitude [A.U.]', yrange=[0, height], /ysty
      loadct, 0, /silent
      oplot, xlengths_Mm, meanamps, psym=2, color=150
      polyfill, [decay_len-dl_err+pixel_size*crop_index, decay_len+dl_err+pixel_size*crop_index, decay_len+dl_err+pixel_size*crop_index, decay_len-dl_err+pixel_size*crop_index, decay_len-dl_err+pixel_size*crop_index], $
        [0, 0, height, height, 0], color=180
      oplot, xlengths_Mm_crop, meanamps_crop, psym=2
      oplot, xlengths_inter, smooth_inter, color=150
      errplot, xlengths_Mm_crop, meanamps_crop-meanamps_err_crop, meanamps_crop+meanamps_err_crop
      oplot, [decay_len+pixel_size*crop_index, decay_len+pixel_size*crop_index], [0, height], linestyle=2
      loadct, 39, /silent
;      oplot, xlengths_inter_crop, fit, color=250, thick=1.5
      if keyword_set(savename) then write_jpeg, strtrim(savename,2)+'.jpeg', TVRD(/TRUE), /TRUE, quality=100
    endelse
  endif

end
