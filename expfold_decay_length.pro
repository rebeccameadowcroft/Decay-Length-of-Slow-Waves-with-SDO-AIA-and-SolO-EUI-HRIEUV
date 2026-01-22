; EXPONENTIAL FOLD DECAY LENGTH
;
; This procedure calculates the decay length of a periodic signal by:
; - filtering the data in time,
; - extracting amplitude profiles over distance,
; - smoothing the data, subtracting background noise
; - Calculating the distance where the maximum value decreases by a factor of e
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
pro expfold_decay_length, data, period, cad_s, pixel_size, max_index, e_fold_len, e_fold_err, $
  VISUALISE=visualise, SAVENAME=savename, PRINT=print

  cad_min = cad_s / 60.

  ; FILTER --------------------------------------------------------------
  ff = 1. / (period * 60.)
  fwidth = ff / 5.
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
    filtsig, tt, xx, [ff - fwidth / 2, ff + fwidth / 2], xf, sd, filt='ideal'
    f_data[*, i] = xf
  endfor

  ; Amplitude extraction and smoothing
  ampslens, f_data, cad_s, meanamps, xlengths, allamps, alllens, 50./cad_s, mean_level, sig_level, pixel_size
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
  width = round(3/pixel_size)
  smooth_meanamps = smooth(meanamps, width, /edge_mirror)
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

  ;find e-folding length (subtracting mean)
  smooth_unmean=smooth_inter-mean_level
  e_height=max(smooth_unmean)/2.71828
  BB=where(smooth_unmean lt e_height)
  e_fold_len=xlengths_inter[BB[0]]
  ;error of e-fold length
  
  


  if keyword_set(print) then begin
    print, 'e-folding length = ' + strtrim(e_fold_len,2) + ' Mm'
  endif


  ;PLOTTING -------------------------------------------------------------------------------------------------------------------
  if keyword_set(visualise) then begin
    if visualise eq 'in_cur_win' then begin
      height = max(meanamps + meanamps_err) * 1.1
      plot, xlengths_Mm_crop, meanamps_crop, psym=2, title='Decay Length Determination', $
        xtitle='Distance [Mm]', ytitle='Amplitude [A.U.]', yrange=[0, height], xrange=[0, max(xlengths_Mm_crop)], /ysty, /xsty
      loadct, 0, /silent
      oplot, xlengths_Mm, meanamps, psym=2, color=50
      oplot, xlengths_Mm_crop, meanamps_crop, psym=2
      loadct, 39, /silent
      oplot, xlengths_inter, smooth_inter, color=50, thick=2
      errplot, xlengths_Mm_crop, meanamps_crop-meanamps_err_crop, meanamps_crop+meanamps_err_crop
      loadct, 39, /silent
      oplot, [0, max(xlengths_Mm_crop)], [mean_level, mean_level], color=50 , linestyle=2
      oplot, [0, max(xlengths_Mm_crop)], [e_height+mean_level, e_height+mean_level], color=50, linestyle=2
      oplot, [e_fold_len, e_fold_len], [0, height], color=250, thick=2
    endif else begin
      ; CASE 2: Open new window with all 3 plots
      wdef, visualise, 600, 800
      !p.multi = [0,1,3]
      !p.charsize = 2.5
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
