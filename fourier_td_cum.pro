; ---------------------------------------------------------------------------
; CUMULATIVE FOURIER ANALYSIS
;
; Procedure to analyze time-distance (TD) data and generate a cumulative
; Fourier power spectrum. Optionally detrends or autocorrelates the data,
; and identifies the dominant period with error estimates.
;
; INPUTS:
;   data         - 2D array of time-distance data [time, distance].
;   pixel_size   - Spatial scale per pixel in Mm.
;   cadence      - Time between frames in seconds.
;
; OUTPUTS:
;   period       - Dominant period in minutes.
;   perr_p       - Positive error on the dominant period [1σ].
;   perr_m       - Negative error on the dominant period [1σ].
;   waveno       - AIA waveband number for plot coloring.
;
; OPTIONAL KEYWORDS:
;   DETREND      - Window width in seconds for removing long-term trends.
;   AC           - If set, computes the auto-correlation of the TD map.
;   SPECTRUM     - Returns 2D array: [frequencies, fourier_power].
;   VISUALISE    - If numeric, opens a new window with all plots.
;                  If 'in_cur_win', only plots Fourier spectrum in current window.
;   SAVENAME     - String to save Fourier plot as JPEG.
;
; ---------------------------------------------------------------------------

PRO fourier_td_cum, td_map, pixel_size, cadence, period, perr_p, perr_m, waveno, $
  detrend=detrend, AC=AC, spectrum, VISUALISE=visualise, savename=savename, half=half
  
  data=td_map
;  if keyword_set(half) then data=data[*,0:n_elements(data[0,*])/2.]
  data=data[*,0:n_elements(data[0,*])/3.]
  
  ; Compute time array for data sequence
  cad_s = cadence
  cad_min = cad_s / 60.
  n = size(data[*,0], /n_elements)
  Td = n * cad_s
  tt = findgen(n) * Td / (n-1)

  ; Apply detrending if requested
  if arg_present(detrend) then begin
    smooth_window = detrend / cadence
    data = data - smooth(data, [smooth_window,1])
  endif

  data_ac = data

  ; Apply auto-correlation if requested
  if arg_present(AC) then begin
    for i = 0, n_elements(data[0,*])-1 do begin
      att = indgen(n)
      data_ac[*,i] = a_correlate(data[*,i], att)
    endfor
  endif

  ; Fourier Transform: Compute cumulative Fourier spectrum, determine dominant period and error
  xx0 = data_ac[*,0]
  temp0 = fft(xx0)
  mag0 = abs(temp0[0:n/2])^2
  freq = findgen(n/2+1) / (n * Td / (n-1)) ; Frequencies in Hz.
  fourier0 = (n * mag0 / variance(xx0))    ; Normalized power spectrum.

  for i = 1, n_elements(data_ac[0,*])-1 do begin
    temp = fft(data_ac[*,i])
    mag = abs(temp[0:n/2])^2
    fourier0 += (n * mag / variance(data_ac[*,i])) ; Accumulate power.
  endfor
  
  spectrum=[[freq], [fourier0]]
  
  sd = stddev(fourier0)
  sd3=3*sd
  sd5=5*sd
  sd10=10*sd
  
  ; Identify peak frequency in range 0–0.02 Hz
  freq_range = where(freq ge 0 and freq le 0.02, count)
  if count gt 0 then begin
    fourier_sub = fourier0[freq_range]
    freq_sub = freq[freq_range]

    fourier_sub = interpol(fourier_sub, n_elements(fourier_sub)*100)*100
    freq_sub = interpol(freq_sub, n_elements(fourier_sub)*100)*100

    max_value = max(fourier_sub, peak_index)
    peak_freq = freq_sub[peak_index]

    left_half = fourier_sub[0:peak_index]
    left_freq = freq_sub[0:peak_index]
    right_half = fourier_sub[peak_index:*]
    right_freq = freq_sub[peak_index:*]

    half_max = max_value / 2.0
    left_index = max(where(left_half le half_max, left_count))
    right_index = min(where(right_half le half_max, right_count))

    if left_count gt 0 and right_count gt 0 then begin
;      print, left_freq
      lower_freq = left_freq[left_index]
      upper_freq = right_freq[right_index]
      fwhm = upper_freq - lower_freq
    endif else begin
      fwhm = !values.f_nan
    endelse
  endif else begin
    print, 'No valid data in the specified frequency range.'
  endelse

  ; Estimate period and its uncertainties
  period = 1 / peak_freq / 60.
  period_upper = 1 / lower_freq / 60.
  period_lower = 1 / upper_freq / 60.
;  period_upper = 0
;  period_lower = 0
  perr_m = [period_upper - period]
  perr_p = [period - period_lower]
;  print, 'Period: ' + strtrim(period, 2) + ' + ' + strtrim(perr_p, 2) + ' - ' + strtrim(perr_m, 2) + ' min'
  
  
  ; Plot if visualise is set
  if keyword_set(visualise) then begin
    if visualise eq 'in_cur_win' then begin
;      print, 'peak_freq = '+strtrim(peak_freq)
;      print, 'max_value = '+strtrim(max_value)
      ymax=max(fourier0)*1.05
      plot, freq, fourier0, xtitle='Frequency [Hz]', ytitle='Fourier Power [a.u.]', $
        title='Cumulative Fourier Spectrum', xrange=[0, 0.02], yrange=[0, ymax], /ystyl
      loadct, 0, /silent
      f_min=peak_freq-peak_freq/10.
      f_max=peak_freq+peak_freq/10.
      polyfill, [f_min, f_max, f_max, f_min, f_min], [0, 0, ymax, ymax, 0], color=200
      oplot, freq, fourier0
      loadct, 39, /silent
      oplot, [0, 0.02], [sd3, sd3], color=200
      oplot, [0, 0.02], [sd5, sd5], color=120
      oplot, [0, 0.02], [sd, sd], color=250
      oplot, [0, 0.02], [sd10, sd10], color=100
      oplot, [peak_freq], [max_value/100.], psym=2, color=250, symsize=1
      oplot, [lower_freq], [max_value/200.], psym=2, color=250, symsize=1
      oplot, [upper_freq], [max_value/200.], psym=2, color=250, symsize=1
      
    endif else begin
      wdef, visualise, 600, 800
      !p.multi = [0, 1, 3]
      !p.background = 255
      !p.color = 0

      aia_lct, wave=waveno, /load
      plot_image, data, xtitle='Time [s]', ytitle='Distance [Mm]', $
        title='Detrended TD Map', scale=[cad_s/60., pixel_size], /nosquare

      loadct, 39, /silent
      plot_image, data_ac, xtitle='Time [s]', ytitle='Distance [Mm]', $
        title='Auto-correlated TD Map', scale=[cad_s, pixel_size], /nosquare

      plot, freq, fourier0, xtitle='Frequency [Hz]', ytitle='Fourier Power [a.u.]', $
        title='Cumulative Fourier Spectrum', xrange=[0, 0.02]
      oplot, [0, 0.02], [sd3, sd3], color=200
      oplot, [0, 0.02], [sd5, sd5], color=120
      oplot, [0, 0.02], [sd, sd], color=250
      oplot, [0, 0.02], [sd10, sd10], color=100
      oplot, [peak_freq], [max_value], psym=2, color=250, symsize=2
;      oplot, [lower_freq], [max_value/2.], psym=2, color=250, symsize=2
;      oplot, [upper_freq], [max_value/2.], psym=2, color=250, symsize=2
    endelse
  endif

  ; Save JPEG image if requested
  if keyword_set(savename) then $
    write_jpeg, 'fourier_cum_' + savename + '.jpg', tvrd(/true), /true, quality=100

END
