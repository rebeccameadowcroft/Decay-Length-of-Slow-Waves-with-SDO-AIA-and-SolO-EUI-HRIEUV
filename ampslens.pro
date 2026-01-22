;
; Function: maxmin
; Description: Finds the local maxima (or minima if the `minima` keyword is set) 
;              in a given dataset.
;
; Input:
;   xx - Array of numerical data.
;
; Output:
;   index - Array of indices where local extrema occur.
;   value - Array of corresponding values of the extrema.
;
function maxmin, xx, index, value;, minima=minima

  index = []  ; Initialize empty arrays to store extrema indices and values
  value = []

  ; If keyword `minima` is set, flip the sign of the input data to find minima instead of maxima
  if keyword_set(minima) then xx = -xx
  
  ; Check if the dataset is too small to analyze
  if n_elements(xx) le 2 then begin
    print, 'Data too small'
    return, [[index], [value]]  ; Return empty arrays
  endif
  
  ; Iterate through the dataset to find local extrema
  for i=1, n_elements(xx)-2 do begin
    ; Check if the current point is a local maximum
    if ( (xx[i-1] lt xx[i]) AND (xx[i] gt xx[i+1])) then begin
      index =[index, i]
      value=[value, xx[i]]
    endif
  endfor
  return, [[index], [value]]  ; Return extrema indices and values
end

;
; Function: ampslens
; Description: Computes mean amplitude of oscillations from a 2D dataset
;
; Input:
;   data - 2D array of numerical data (time series data for multiple signals)
;   width - Smoothing width for the signal
;   wind - Window ID for optional plotting
;
; Output:
;   meanamps - Array of mean amplitudes for each signal
;   xlengths - Array of corresponding signal indices
;   allamps  - Array of all computed amplitudes
;   alllens  - Array of indices corresponding to `allamps`
;
pro ampslens, data, cad_s, meanamps, xlengths, allamps, alllens, mean_level, sig_level, pixel_size
  meanamps=[]  ; Initialize output arrays
  xlengths=[]
  allamps=[]
  alllens=[]
  stds=[]
  
  sz = size(data)   ; Get data dimensions
  len = sz[2]       ; Number of signals
  
  for i = 0, len-1 do begin
    xxs=data[*,i]
;    xxs = smooth(data[*,i], width)  ; Smooth the signal
;    xxs = xxs-smooth(xxs, 300/cad_s, /edge_mirror)
;    if xxs[0] ne NaN then begin
      ; Find local maxima and minima
      maxes=maxmin(xxs, max_index0, max_value0)
      mines=maxmin(-xxs, min_index0, min_value0)
      min_value0 = -min_value0  ; Convert minima back to original values
      
      ; Calculate amplitudes
      maxs = size(max_value0, /dim)
      mins = size(min_value0, /dim)
      pairs = min(maxs, mins)  ; Ensure equal number of maxima and minima
      
      amps = []
      for n = 0, pairs-2 do begin
        diff = abs(max_value0[n] - min_value0[n])  ; Peak-to-trough difference
        amp = diff  ; Amplitude calculation (can be modified if needed)
        amps = [amps, amp]
        allamps = [allamps, amp]
        alllens = [alllens, i]
      endfor
      
      meanamp = mean(amps)  ; Compute mean amplitude
      stdamp = stddev(amps)
      stds = [stds, stdamp]
      meanamps = [meanamps, meanamp]
      xlengths = [xlengths, i]  ; Store signal index    
;      ;plotting to work out WHAT IS GOING ON HERE
;      if i eq 5 then begin
;        plot, xxs, /ystyl, xtitle='Time [frames]', ytitle='Intensity'
;        loadct, 39, /silent
;        oplot, max_index0, max_value0, psym=4, color=250
;        oplot, min_index0, min_value0, psym=4, color=50
;        print, 'plotted'
;      endif
    ;endif
  endfor
  
;  loadct, 7, /silent
;  plot, alllens, allamps, psym=4, xtitle='Distance [pix]', ytitle='Amplitude'
;  errplot, xlengths, meanamps-stds, meanamps+stds, color=150
;  oplot, xlengths, meanamps, psym=4, color=150, thick=2
  
  N=2.0/pixel_size
  tail = meanamps[-N:-1]
  sm_tail = smooth(tail, 3, /edge_mirror)  ; small running smoothing
  mean_level = median(sm_tail)
  sig_level = mean_level + 3.0 * stddev(sm_tail)
  
;  print, 'top mean = ', top_mean
;  print, 'std = ', std
;  print, 'SIG LEVEL = ',sig_level

end
