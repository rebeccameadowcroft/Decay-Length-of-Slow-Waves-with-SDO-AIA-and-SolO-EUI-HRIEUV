;----------------------------------------------------------------------------------
; Quasi-Fan Analysis Script
;
; This script aligns AIA and EUI datasets in time, compensates for light travel time,
; creates visualizations and movies of the aligned datasets, generates running and
; mean-difference image cubes, selects slits for time-distance (TD) analysis, and 
; computes wave parameters (period, speed, decay length) along those slits using 
; the `per_spe_dl` routine.
;
; Input: AIA and EUI image cubes and corresponding index structures.
; Output: JPEGs for visual inspection and .SAV files for TD maps and analysis.
;----------------------------------------------------------------------------------



function make_tracks, x_coords, y_coords
  npt=5
  x=interpol(x_coords,10*npt,/spline)      ;NOTE: spline keyword plots a cubic splice to the 4 point neighborhood surrounding interval
  y=interpol(y_coords,10*npt,/spline)
  dis1=dblarr(10*npt-1)
  for i=1,10*5-1 do dis1[i-1]=sqrt((x[i]-x[i-1])^2.+(y[i]-y[i-1])^2.)     ;pythagorus of distances between subsequent points
  dis=total(dis1)                                                           ;length of the selected curve line, in the unit of spatial pixel
  xtrack=interpol(x,dis,/spline)                                            ;interpolating for x and y values for each pixel of total line distance
  ytrack=interpol(y,dis,/spline)
  return, [[xtrack], [ytrack]]
end

cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan'

; Load pre-selected AIA data and index
AIA_index = AIA_index_a
;AIA_data = AIA_data_a
EUI_index = EUI_index_a
;EUI_data = EUI_data_a

aia_cad=12.
eui_cad=3.

; Compute AIA and EUI pixel sizes (in radians per pixel)
aia_pixel_size = mean(aia_index.dsun_obs)/1e6 * mean(aia_index.cdelt1) * !dpi/180./3600.
eui_pixel_size = mean(eui_index.dsun_obs)/1e6 * mean(eui_index.cdelt1) * !dpi/180./3600.


;; ALIGNING DATA IN TIME
;;----------------------------------------------------------------------------------
;; This section corrects for start/end time mismatches and accounts for light travel 
;; time differences between AIA and EUI data. It includes code to crop frames and
;; save the aligned data. It is currently commented out since manual alignment was done.
;;Correct for different start/end times------------------------------------------------
;print, ' '
;print, '---Correcting for different start/end times---'
;print,'AIA start time = '+ AIA_index[0].date_d$obs
;print,'EUI start time = '+ EUI_index[0].date_d$obs
;print,'AIA end time = '+ AIA_index[-1].date_d$obs
;print,'EUI end time = '+ EUI_index[-1].date_d$obs
;
;;;Select manually which frames to use
;;;;New start time = 03:50:09
;aia_data_a = aia_data[*,*,*] & aia_index_a=aia_index[*]
;eui_data_a = eui_data[*,*,3:-1] & eui_index_a=eui_index[3:-1]
;
;
;print, ' '
;print,'Corrected AIA start time = '+ AIA_index_a[0].date_d$obs
;print,'Corrected EUI start time = '+ EUI_index_a[0].date_d$obs
;print,'Corrected AIA end time = '+ AIA_index_a[-1].date_d$obs
;print,'Corrected EUI end time = '+ EUI_index_a[-1].date_d$obs
;
;;Light travel time + aligning data---------------------------------------------------
;AIA_dist_m = mean(AIA_index_a.dsun_obs)
;EUI_dist_m = mean(EUI_index_a.dsun_obs)
;ddist = AIA_dist_m-EUI_dist_m
;dtime = ddist/(3.0e+8)
;AIA_dframes =  round(dtime/12.)
;EUI_dframes =  round(dtime/3.)
;print, ' '
;print, '---Correcting for light travel time---'
;print, 'Travel time = '+ strtrim(dtime)
;print, 'AIA Frames to crop = '+strtrim(AIA_dframes)
;print, 'EUI Frames to crop = '+strtrim(EUI_dframes)
;
;aia_data_a = aia_data_a[*,*,AIA_dframes:-1] & aia_index_a=aia_index_a[AIA_dframes:-1]
;eui_data_a = eui_data_a[*,*,0:-1-EUI_dframes] & eui_index_a=eui_index_a[0:-1-EUI_dframes]
;
;print, ' '
;print,'New AIA start time = '+ AIA_index_a[0].date_d$obs
;print,'New EUI start time = '+ EUI_index_a[0].date_d$obs
;print,'New AIA end time = '+ AIA_index_a[-1].date_d$obs
;print,'New EUI end time = '+ EUI_index_a[-1].date_d$obs
;
;save, filename='AIA_data_a.sav', aia_data_a
;save, filename='AIA_index_a.sav', aia_index_a
;save, filename='EUI_index_a.sav', eui_index_a
;save, filename='EUI_data_a.sav', eui_data_a


;; MOVIE CREATION: AIA and EUI side-by-side
;;----------------------------------------------------------------------------------
;; Creates visual comparisons of AIA and EUI images aligned in time.
;wdef, 0, 1300, 650
;!p.multi=[0, 2, 1]
;!p.charsize=1.7
;!p.charthick=1
;!p.background=255
;!p.color=0
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/MOVIES/AIA and EUI'
;;aia_image=comprange(aia_data_a)
;;eui_image=comprange(eui_data_a)
;;cube_difference, aia_data_a, aia_data_f, 'pfilter'
;;cube_difference, eui_data_a, eui_data_f, 'pfilter'
;cube_difference, aia_data_c, aia_data_f, 'pfilter', cad=aia_cad
;cube_difference, eui_data_c, eui_data_f, 'pfilter', cad=eui_cad
;aia_image=aia_data_f
;eui_image=eui_data_f
;
;aia_lct, wave=171, /load
;for i = 0, n_elements(aia_data_a[0,0,*])-2 do begin
;  for j = 0, 3 do begin    
;    aia_plot=comprange(aia_image[*,*,i])
;    plot_image, bytscl(aia_plot, 0.25, 0.75), scale=[aia_pixel_size, aia_pixel_size], $
;      title='AIA RoI '+aia_index_a[i].date_D$OBS+' UT', xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]'
;    eui_plot=comprange(eui_image[*,*,4*i+j])
;    plot_image, bytscl(eui_plot, 0.25, 0.75), scale=[eui_pixel_size, eui_pixel_size], $
;      title='EUI RoI '+eui_index_a[4*i+j].date_D$OBS+' UT', xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]'
;    write_jpeg, 'AIA_and_EUI'+strtrim(4*i+j, 2)+'.jpg', TVRD(/TRUE), /TRUE, quality=100
;  endfor
;endfor


;; CREATE MEAN AND RUNNING DIFFERENCE IMAGES
;;----------------------------------------------------------------------------------
;; Applies 'running' and 'mean' difference techniques to emphasize dynamic features
;; in AIA and EUI image sequences.
;
;cube_difference, aia_data_a, aia_data_running, 'running'
;cube_difference, aia_data_a, aia_data_mean, 'mean'
;
;cube_difference, eui_data_a, eui_data_running, 'running'
;cube_difference, eui_data_a, eui_data_mean, 'mean'
;
;wdef, 0, 1300, 1000
;!p.multi=[0, 3, 2]
;!p.charsize=1.5
;!p.background=255
;!p.color=0
;aia_lct, wave=171, /load
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/MOVIES/Differences'
;;aia_image_running=comprange(aia_data_running)
;;aia_image_mean=comprange(aia_data_mean)
;;
;;eui_image_running=comprange(eui_data_running)
;;eui_image_mean=comprange(eui_data_mean)
;for i = 0, n_elements(aia_data[0,0,*])-2 do begin
;  for j = 0, 3 do begin
;    plot_image, (aia_data[*,*,i]), title='AIA OG '+aia_index_a[i].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    plot_image, (aia_data_running[*,*,i]), title='AIA Running '+aia_index_a[i].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    plot_image, (aia_data_mean[*,*,i]), title='AIA Mean '+aia_index_a[i].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    
;    plot_image, (eui_data[*,*,i*4+j]), title='EUI OG '+eui_index_a[4*i+j].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    plot_image, (eui_data_running[*,*,i*4+j]), title='EUI Running '+eui_index_a[4*i+j].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    plot_image, (eui_data_mean[*,*,i*4+j]), title='EUI Mean '+eui_index_a[4*i+j].date_D$OBS+' UT', xtitle='Solar X [pix]', ytitle='Solar Y [pix]'
;    write_jpeg, 'AIA_and_EUI_differences'+strtrim(4*i+j, 2)+'.jpg', TVRD(/TRUE), /TRUE, quality=100
;  endfor
;endfor


;; INTERACTIVE SLIT SELECTION & ANALYSIS LOOP (FOR 5 SLITS) - STRAIGHT SLITS
;;----------------------------------------------------------------------------------
;
;; Cropped AIA region
aia_crop_x1 = 70
aia_crop_x2 = 180
aia_crop_y1 = 50
aia_crop_y2 = 160
aia_small = aia_data_a[aia_crop_x1:aia_crop_x2, aia_crop_y1:aia_crop_y2, *]
;;
;; Cropped EUI region
eui_crop_x1 = 110
eui_crop_x2 = 300
eui_crop_y1 = 90
eui_crop_y2 = 280
eui_small = eui_data_a[eui_crop_x1:eui_crop_x2, eui_crop_y1:eui_crop_y2, *]
;; Prepare storage
;x1s = [] & y1s = [] & x2s = [] & y2s = []
;
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Analysis'
;
;for i = 0, 4 do begin
;  done = 0
;  while done eq 0 do begin
;    ;-----------------------------
;    ; AIA Slit Selection
;    ;-----------------------------
;    !p.multi = [0, 2, 1]
;    loadct, 39, /silent
;
;    wdef, 0, 1600, 900
;    aia_lct, wave=171, /load
;    plot_image, alog(aia_small[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', title='Zoomed AIA 171'
;
;    print, 'Select start and end points for slit ', i+1, ' in AIA (left click). Right click when happy.'
;    cursor, x1a, y1a, /down, /data
;    oplot, [x1a], [y1a], color=50, psym=2
;    cursor, x2a, y2a, /down, /data
;    oplot, [x2a], [y2a], color=50, psym=2
;    oplot, [x1a, x2a], [y1a, y2a], color=50, thick=2
;    xyouts, x2a-20, y2a, '('+strtrim(i+1,2)+')', charsize=1.5
;
;    ; Offset to full frame
;    x1_full_aia = x1a + aia_crop_x1
;    y1_full_aia = y1a + aia_crop_y1
;    x2_full_aia = x2a + aia_crop_x1
;    y2_full_aia = y2a + aia_crop_y1
;
;    ;-----------------------------
;    ; EUI Slit Selection
;    ;-----------------------------
;    aia_lct, wave=171, /load
;    plot_image, alog(eui_small[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', title='Zoomed EUI 174'
;
;    print, 'Select start and end points for slit ', i+1, ' in EUI (left click). Right click when happy.'
;    cursor, x1e, y1e, /down, /data
;    oplot, [x1e], [y1e], color=50, psym=2
;    cursor, x2e, y2e, /down, /data
;    oplot, [x2e], [y2e], color=50, psym=2
;    oplot, [x1e, x2e], [y1e, y2e], color=50, thick=2
;    xyouts, x2e-20, y2e, '('+strtrim(i+1,2)+')', charsize=1.5
;
;    ; Offset to full frame
;    x1_full_eui = x1e + eui_crop_x1
;    y1_full_eui = y1e + eui_crop_y1
;    x2_full_eui = x2e + eui_crop_x1
;    y2_full_eui = y2e + eui_crop_y1
;
;    ;-----------------------------
;    ; Perform Analysis
;    ;-----------------------------
;    aia_td = timedistance(aia_data_mean, [x1_full_aia, y1_full_aia, x2_full_aia, y2_full_aia], width=3)
;    eui_td = timedistance(eui_data_mean, [x1_full_eui, y1_full_eui, x2_full_eui, y2_full_eui], width=3)
;    
;    ; Run main analysis and display results
;    per_spe_dl, aia_td, aia_index_a, aia_pixel_size, 12., $
;      eui_td, eui_index_a, eui_pixel_size, 3., 4, $
;      SAVENAME='slit_' + strtrim(i+1, 2)
;
;    ;-----------------------------
;    ; Check if user is happy
;    ;-----------------------------
;    print, 'Right-click to confirm this slit. Left-click to reselect.'
;    cursor, junkx, junky, /data
;    if !mouse.button eq 4 then begin
;      ; Right-click = accept
;      x1s = [x1s, x1_full_aia]
;      y1s = [y1s, y1_full_aia]
;      x2s = [x2s, x2_full_aia]
;      y2s = [y2s, y2_full_aia]
;      x1s = [x1s, x1_full_eui]
;      y1s = [y1s, y1_full_eui]
;      x2s = [x2s, x2_full_eui]
;      y2s = [y2s, y2_full_eui]
;      done = 1
;    endif else begin
;      print, 'Reselecting slit ', i+1
;    endelse
;  endwhile
;endfor
;;-------------------------------------
;; Plot AIA and EUI FoV with all 5 selected STRAIGHT slits
;;-------------------------------------
;wdef, 2, 1600, 800
;!p.multi = [0, 2, 1]  ; Two horizontal panels
;
;; AIA Plot
;aia_lct, wave=171, /load
;plot_image, alog(aia_small[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='Zoomed AIA 171 with Slits';, scale=[aia_pixel_size, aia_pixel_size]
;for i = 0, 4 do begin
;  x1 = x1s[2*i] - aia_crop_x1
;  y1 = y1s[2*i] - aia_crop_y1
;  x2 = x2s[2*i] - aia_crop_x1
;  y2 = y2s[2*i] - aia_crop_y1
;  loadct, 39, /silent
;  oplot, [x1, x2], [y1, y2], color=100, thick=2
;  xyouts, x2+1, y2, '('+strtrim(i+1,2)+')', charsize=1.5, color=255, charthick=2
;endfor
;
;; EUI Plot
;aia_lct, wave=171, /load
;plot_image, alog(eui_small[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='Zoomed EUI 174 with Slits';, scale=[eui_pixel_size, eui_pixel_size]
;for i = 0, 4 do begin
;  x1 = x1s[2*i+1] - eui_crop_x1
;  y1 = y1s[2*i+1] - eui_crop_y1
;  x2 = x2s[2*i+1] - eui_crop_x1
;  y2 = y2s[2*i+1] - eui_crop_y1
;  loadct, 39, /silent
;  oplot, [x1, x2], [y1, y2], color=100, thick=2
;  xyouts, x2+2, y2, '('+strtrim(i+1,2)+')', charsize=1.5, color=255, charthick=2
;endfor
;
;!p.multi = 0


;; Slit selection - curved slits
;;-------------------------------------
;; Prepare storage
;fixed_aia_xs = [] & fixed_aia_ys = []
;fixed_eui_xs = [] & fixed_eui_ys = []
;aia_td_1=[] & eui_td_1=[] & aia_xs_1=[] & aia_ys_1=[] & eui_xs_1=[] & eui_ys_1=[] 
;aia_td_2=[] & eui_td_2=[] & aia_xs_2=[] & aia_ys_2=[] & eui_xs_2=[] & eui_ys_2=[] 
;aia_td_3=[] & eui_td_3=[] & aia_xs_3=[] & aia_ys_3=[] & eui_xs_3=[] & eui_ys_3=[] 
;aia_td_4=[] & eui_td_4=[] & aia_xs_4=[] & aia_ys_4=[] & eui_xs_4=[] & eui_ys_4=[] 
;aia_td_5=[] & eui_td_5=[] & aia_xs_5=[] & aia_ys_5=[] & eui_xs_5=[] & eui_ys_5=[] 
;
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Analysis'
;for i = 1, 4 do begin
;  done = 0
;  while done eq 0 do begin
;    ;-----------------------------
;    ; AIA Slit Selection
;    ;-----------------------------
;    stplot_curve, aia_data_a, sm_scl=250./aia_cad, ncut=1, width=3, winnum=0
;    restore, filename='stplot_cut=0.sav'
;    restore, filename='stplot_cutslocation.sav'
;    aia_xs=xi
;    aia_ys=yi
;    aia_td=transpose(int_slice)
;
;    ;-----------------------------
;    ; EUI Slit Selection
;    ;-----------------------------
;    stplot_curve, eui_data_a, sm_scl=250./eui_cad, ncut=1, width=5, winnum=2
;    restore, filename='stplot_cut=0.sav'
;    restore, filename='stplot_cutslocation.sav'
;    eui_xs=xi
;    eui_ys=yi
;    eui_td=transpose(int_slice)
;
;    ;-----------------------------
;    ; Perform Analysis
;    ;-----------------------------
;    per_spe_dl, aia_td, aia_index_a, aia_pixel_size, 12., $
;      eui_td, eui_index_a, eui_pixel_size, 3., 4, $
;      SAVENAME='slit_' + strtrim(i+1, 2)
;
;    ;-----------------------------
;    ; Check if user is happy
;    ;-----------------------------
;    print, 'Right-click to confirm this slit. Left-click to reselect.'
;    cursor, junkx, junky, /data
;    if !mouse.button eq 4 then begin
;      ; Right-click = accept
;      done = 1
;      fixed_aia_xs=[fixed_aia_xs, aia_xs] & fixed_aia_ys=[fixed_aia_ys, aia_ys]
;      fixed_eui_xs=[fixed_eui_xs, eui_xs] & fixed_eui_ys=[fixed_eui_ys, eui_ys]
;    endif else begin
;      print, 'Reselecting slit ', i+1
;    endelse
;  endwhile
;endfor


;; ALIGNING DATA IN SPACE
;;----------------------------------------------------------------------------------
;
;;Choosing AIA crop region in pixels
aia_width=120
aia_crop_x1=65
aia_crop_x2=aia_crop_x1+aia_width
aia_crop_y1=35
aia_crop_y2=aia_crop_y1+aia_width
aia_data_c=aia_data_a[aia_crop_x1:aia_crop_x2, aia_crop_y1:aia_crop_y2,*]

ratio=aia_pixel_size/eui_pixel_size

eui_width=(aia_crop_x2-aia_crop_x1)*ratio
eui_crop_x1=50
eui_crop_x2=eui_crop_x1+eui_width
eui_crop_y1=20
eui_crop_y2=eui_crop_y1+eui_width
eui_data_c=eui_data_a[eui_crop_x1:eui_crop_x2, eui_crop_y1:eui_crop_y2,*]
;
;
;; Analysis with set slits
;;-------------------------------------
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Analysis'
;restore, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Final Figures/fixed_slit_positions.sav'
fixed_aia_xs=[113.113, 116.875, 119.508, 121.013, 120.637,$
  115.74627, 121.01315, 126.65624, 131.17071, 137.94241,$
  125.15141, 130.79450, 136.06138, 141.32826, 145.84273,$
  127.40865, 134.55656, 140.57585, 148.85238, 156.75270]

fixed_aia_ys=[113.489, 116.875, 122.142, 130.042, 138.695, $
  113.48904, 117.62730, 124.39901, 130.41830, 139.82344,$
  92.797725, 80.759141, 69.849175, 58.186797, 43.138568,$
  96.559782, 89.411873, 83.768787, 73.987438, 64.206089]

fixed_eui_xs=[165.829, 179.368, 190.327, 200.642, 207.734,$
  169.05240, 178.72282, 189.03793, 198.70835, 209.02346,$
  170.98649, 183.23568, 198.06366, 212.89163, 226.43021,$
  183.88038, 199.35305, 214.18102, 230.94308, 248.34983]

fixed_eui_ys=[194.840, 201.932, 207.734, 218.049, 232.232,$
  194.19550, 194.84019, 198.06366, 203.22122, 212.24694,$
  160.67139, 145.19872, 125.85789, 106.51705, 85.242138,$
  156.80322, 144.55403, 130.37075, 112.31930, 96.201943]

crop_aia_xs=fixed_aia_xs-aia_crop_x1
crop_aia_ys=fixed_aia_ys-aia_crop_y1
crop_eui_xs=fixed_eui_xs-eui_crop_x1
crop_eui_ys=fixed_eui_ys-eui_crop_y1
;
;;Dont need to run below part, just restore 'all_td_maps.sav'
;!p.background=255
;!p.color=0
;!p.charthick=1
;aia_lct, wave=171, /load
;for i = 0 , 4 do begin
;  ;create td map
;  xpoints=crop_aia_xs[i*5: i*5+4]
;  ypoints=crop_aia_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, aia_data_c, int_slice, sm_scl=250./aia_cad, ncut=1, width=3, winnum=0, points=points
;  aia_td=transpose(int_slice)
;  xpoints=crop_eui_xs[i*5: i*5+4]
;  ypoints=crop_eui_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, eui_data_c, int_slice, sm_scl=250./eui_cad, ncut=1, width=5, winnum=0, points=points
;  eui_td=transpose(int_slice)
;  save, filename='aia_eui_td_maps_'+strtrim(i+1,2)+'.sav', aia_td, eui_td
;endfor
;
;;restore and save all td map info
;path='/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Analysis/'
;restore, 'aia_eui_td_maps_1.sav'
;aia_td_1=[aia_td] & eui_td_1=[eui_td]
;restore, 'aia_eui_td_maps_2.sav'
;aia_td_2=[aia_td] & eui_td_2=[eui_td]
;restore, 'aia_eui_td_maps_3.sav'
;aia_td_3=[aia_td] & eui_td_3=[eui_td]
;restore, 'aia_eui_td_maps_4.sav'
;aia_td_4=[aia_td] & eui_td_4=[eui_td]
;restore, 'aia_eui_td_maps_5.sav'
;aia_td_5=[aia_td] & eui_td_5=[eui_td]
;
;save, filename='all_td_maps.sav', aia_td_1, aia_td_2, aia_td_3, aia_td_4, aia_td_5, $
;      eui_td_1, eui_td_2, eui_td_3, eui_td_4, eui_td_5
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Analysis'
;;restore, 'all_td_maps.sav'
;crop_options=[1,1,2,1,3]
;dl_regions=[5,25, 20,40, 5,25, 0,20, 10,30 ]
;
;is=[1, 2, 3, 5]
;;is =[2]
;for j = 0, 3 do begin
;  i=is[j]
;  aia_name='aia_td_' + STRTRIM(i, 2) & eui_name='eui_td_'+strtrim(i,2)
;  aia_map=scope_varfetch(aia_name) & eui_map=scope_varfetch(eui_name)
;  per_spe_dl, aia_map, aia_index_c, aia_pixel_size, 12., $
;              eui_map, eui_index_c, eui_pixel_size, 3., 4, $
;              SAVENAME='slit_' + strtrim(i, 2), $
;              crop_option=crop_options[i-1], dl_region = dl_regions[(i-1)*2:(i-1)*2+1]
;;  save, aia_map, eui_map, filename='Slit_'+strtrim(j,2)+'_maps.sav'
;endfor
;
;print, 'ANALYSIS DONE!'

;; Plot AIA and EUI FoV with all 5 selected CURVED slits
;;-------------------------------------
;fixed_aia_xs=[113.113, 116.875, 119.508, 121.013, 120.637,$
;              115.74627, 121.01315, 126.65624, 131.17071, 137.94241,$
;              125.15141, 130.79450, 136.06138, 141.32826, 145.84273,$
;              113.11283, 109.35077, 106.71733, 105.58871, 102.57907,$
;              127.40865, 134.55656, 140.57585, 148.85238, 156.75270]
;
;fixed_aia_ys=[113.489, 116.875, 122.142, 130.042, 138.695, $
;              113.48904, 117.62730, 124.39901, 130.41830, 139.82344,$
;              92.797725, 80.759141, 69.849175, 58.186797, 43.138568,$
;              105.96493, 99.569428, 92.045313, 84.521198, 75.492261,$
;              96.559782, 89.411873, 83.768787, 73.987438, 64.206089]
;
crop_aia_xs=fixed_aia_xs-aia_crop_x1
crop_aia_ys=fixed_aia_ys-aia_crop_y1
;
;fixed_eui_xs=[165.829, 179.368, 190.327, 200.642, 207.734,$
;              169.05240, 178.72282, 189.03793, 198.70835, 209.02346,$
;              170.98649, 183.23568, 198.06366, 212.89163, 226.43021,$
;              152.93504, 158.73729, 161.96077, 169.05240, 175.49935,$
;              183.88038, 199.35305, 214.18102, 230.94308, 248.34983]
;              
;fixed_eui_ys=[194.840, 201.932, 207.734, 218.049, 232.232,$
;              194.19550, 194.84019, 198.06366, 203.22122, 212.24694,$
;              160.67139, 145.19872, 125.85789, 106.51705, 85.242138,$
;              174.20997, 163.25016, 154.22444, 136.81769, 120.05564,$
;              156.80322, 144.55403, 130.37075, 112.31930, 96.201943]
;;Without slit 4
;fixed_aia_xs=[113.113, 116.875, 119.508, 121.013, 120.637,$
;  115.74627, 121.01315, 126.65624, 131.17071, 137.94241,$
;  125.15141, 130.79450, 136.06138, 141.32826, 145.84273,$
;  127.40865, 134.55656, 140.57585, 148.85238, 156.75270]
;
;fixed_aia_ys=[113.489, 116.875, 122.142, 130.042, 138.695, $
;  113.48904, 117.62730, 124.39901, 130.41830, 139.82344,$
;  92.797725, 80.759141, 69.849175, 58.186797, 43.138568,$
;  96.559782, 89.411873, 83.768787, 73.987438, 64.206089]
;
crop_aia_xs=fixed_aia_xs-aia_crop_x1
crop_aia_ys=fixed_aia_ys-aia_crop_y1
;
;fixed_eui_xs=[165.829, 179.368, 190.327, 200.642, 207.734,$
;  169.05240, 178.72282, 189.03793, 198.70835, 209.02346,$
;  170.98649, 183.23568, 198.06366, 212.89163, 226.43021,$
;  183.88038, 199.35305, 214.18102, 230.94308, 248.34983]
;
;fixed_eui_ys=[194.840, 201.932, 207.734, 218.049, 232.232,$
;  194.19550, 194.84019, 198.06366, 203.22122, 212.24694,$
;  160.67139, 145.19872, 125.85789, 106.51705, 85.242138,$
;  156.80322, 144.55403, 130.37075, 112.31930, 96.201943]
;
;
;crop_eui_xs=fixed_eui_xs-eui_crop_x1
;crop_eui_ys=fixed_eui_ys-eui_crop_y1
;
;
;; AIA Plot - full data_a FoV
;aia_lct, wave=171, /load
;wdef, 2, 1300, 700
;!p.multi = [0, 2, 1]  ; Two horizontal panels
;!p.color=0
;!p.charsize=2
;!p.BACKGROUND=255
;plot_image, alog(aia_data_a[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='AIA 171 Slits Selection';, scale=[aia_pixel_size, aia_pixel_size]
;for i = 0, 4 do begin
;  aia_x=fixed_aia_xs[i*5:(i+1)*5-1]
;  aia_y=fixed_aia_ys[i*5:(i+1)*5-1]
;  aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
;  aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
;  aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
;  aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
;  oplot, aia_xtrack, aia_ytrack, color=0, linestyle=0, thick=2
;  oplot, aia_xtrack_m[0:-2], aia_ytrack_m[0:-2], color=0, linestyle=2, thick=1
;  oplot, aia_xtrack_p[0:-2], aia_ytrack_p[0:-2], color=0, linestyle=2, thick=1
;  oplot, [aia_crop_x1, aia_crop_x2, aia_crop_x2, aia_crop_x1, aia_crop_x1], $
;         [aia_crop_y1, aia_crop_y1, aia_crop_y2, aia_crop_y2, aia_crop_y1], color=255, thick=2
;  
;;  xyouts, x[-1]+1, y[-1], '('+strtrim(i+1,2)+')', charsize=1.5, color=255, charthick=2
;endfor
;
;; EUI Plot
;aia_lct, wave=171, /load
;plot_image, alog(eui_data_a[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='EUI 174 Slits Selection';, scale=[eui_pixel_size, eui_pixel_size]
;for i = 0, 4 do begin
;  eui_x=fixed_eui_xs[i*5:(i+1)*5-1]
;  eui_y=fixed_eui_ys[i*5:(i+1)*5-1]
;  eui_tracks=make_tracks(eui_x, eui_y) & eui_xtrack=eui_tracks[*,0] & eui_ytrack=eui_tracks[*,1]
;  eui_track_m1=perp_direc(eui_xtrack, eui_ytrack, -1) & eui_track_p1=perp_direc(eui_xtrack, eui_ytrack, +1)
;  eui_xtrack_m1=(eui_track_m1[*,0]) & eui_ytrack_m1=(eui_track_m1[*,1])
;  eui_xtrack_p1=(eui_track_p1[*,0]) & eui_ytrack_p1=(eui_track_p1[*,1])
;  eui_track_m2=perp_direc(eui_xtrack, eui_ytrack, -2) & eui_track_p2=perp_direc(eui_xtrack, eui_ytrack, +2)
;  eui_xtrack_m2=(eui_track_m2[*,0]) & eui_ytrack_m2=(eui_track_m2[*,1])
;  eui_xtrack_p2=(eui_track_p2[*,0]) & eui_ytrack_p2=(eui_track_p2[*,1])
;  oplot, eui_xtrack, eui_ytrack, color=0, linestyle=0, thick=2
;  oplot, eui_xtrack_m1[0:-2], eui_ytrack_m1[0:-2], color=0, linestyle=2, thick=1
;  oplot, eui_xtrack_p1[0:-2], eui_ytrack_p1[0:-2], color=0, linestyle=2, thick=1
;  oplot, eui_xtrack_m2[0:-2], eui_ytrack_m2[0:-2], color=0, linestyle=2, thick=1
;  oplot, eui_xtrack_p2[0:-2], eui_ytrack_p2[0:-2], color=0, linestyle=2, thick=1
;  oplot, [eui_crop_x1, eui_crop_x2, eui_crop_x2, eui_crop_x1, eui_crop_x1], $
;    [eui_crop_y1, eui_crop_y1, eui_crop_y2, eui_crop_y2, eui_crop_y1], color=255, thick=2
;;  xyouts, x[-1]+1, y[-1], '('+strtrim(i+1,2)+')', charsize=1.5, color=255, charthick=2
;endfor
;
;
;
;; AIA Plot - crop data_c FoV
;aia_lct, wave=171, /load
;wdef, 3, 1300, 700
;!p.multi = [0, 2, 1]  ; Two horizontal panels
;!p.color=0
;!p.charsize=2
;!p.BACKGROUND=255
;plot_image, alog(aia_data_c[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='AIA 171 Slits Selection';, scale=[aia_pixel_size, aia_pixel_size]
;for i = 0, 4 do begin
;  aia_x=crop_aia_xs[i*5:(i+1)*5-1]
;  aia_y=crop_aia_ys[i*5:(i+1)*5-1]
;  aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
;  aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
;  aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
;  aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
;  loadct, 39, /silent
;  slit_col=75
;  oplot, aia_xtrack, aia_ytrack, color=slit_col, linestyle=0, thick=2
;  oplot, aia_xtrack_m[0:-2], aia_ytrack_m[0:-2], color=slit_col, linestyle=2, thick=1
;  oplot, aia_xtrack_p[0:-2], aia_ytrack_p[0:-2], color=slit_col, linestyle=2, thick=1
;  xyouts, aia_x[-1]+1, aia_y[-1]-0.5, '('+strtrim(i+1,2)+')', charsize=2, color=255, charthick=2
;endfor
;
;; EUI Plot
;aia_lct, wave=171, /load
;plot_image, alog(eui_data_c[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
;  title='EUI 174 Slits Selection';, scale=[eui_pixel_size, eui_pixel_size]
;for i = 0, 4 do begin
;  eui_x=crop_eui_xs[i*5:(i+1)*5-1]
;  eui_y=crop_eui_ys[i*5:(i+1)*5-1]
;  eui_tracks=make_tracks(eui_x, eui_y) & eui_xtrack=eui_tracks[*,0] & eui_ytrack=eui_tracks[*,1]
;  eui_track_m1=perp_direc(eui_xtrack, eui_ytrack, -1) & eui_track_p1=perp_direc(eui_xtrack, eui_ytrack, +1)
;  eui_xtrack_m1=(eui_track_m1[*,0]) & eui_ytrack_m1=(eui_track_m1[*,1])
;  eui_xtrack_p1=(eui_track_p1[*,0]) & eui_ytrack_p1=(eui_track_p1[*,1])
;  eui_track_m2=perp_direc(eui_xtrack, eui_ytrack, -2) & eui_track_p2=perp_direc(eui_xtrack, eui_ytrack, +2)
;  eui_xtrack_m2=(eui_track_m2[*,0]) & eui_ytrack_m2=(eui_track_m2[*,1])
;  eui_xtrack_p2=(eui_track_p2[*,0]) & eui_ytrack_p2=(eui_track_p2[*,1])
;  loadct, 39, /silent
;  slit_col=75
;  oplot, eui_xtrack, eui_ytrack, color=slit_col, linestyle=0, thick=2
;  oplot, eui_xtrack_m1[0:-2], eui_ytrack_m1[0:-2], color=slit_col, linestyle=2, thick=1
;  oplot, eui_xtrack_p1[0:-2], eui_ytrack_p1[0:-2], color=slit_col, linestyle=2, thick=1
;  oplot, eui_xtrack_m2[0:-2], eui_ytrack_m2[0:-2], color=slit_col, linestyle=2, thick=1
;  oplot, eui_xtrack_p2[0:-2], eui_ytrack_p2[0:-2], color=slit_col, linestyle=2, thick=1
;  xyouts, eui_x[-1]+2, eui_y[-1]-1, '('+strtrim(i+1,2)+')', charsize=2, color=255, charthick=2
;endfor
;
;!p.multi = 0


;;;-----------------------QUICK RUN ANALYSIS HERE!!!!!!-----------------------------------------
;;restore, 'all_td_maps.sav'
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response'
;crop_options=[1,1,2,3]
;;dl_regions=[5,25, 20,40, 5,25, 0,20, 10,30 ]
;dl_regions=[5,25, 5,25, 5,25, 5,25]
;aia_pixel_size=0.431 & aia_cad=12.
;eui_pixel_size=0.194 & eui_cad=3.
;
;for i = 1, 4 do begin
;  crop_option=crop_options[i-1]
;  aia_name='event1_aia_td_' + STRTRIM(i, 2) & eui_name='event1_eui_td_'+strtrim(i,2)
;  aia_map=scope_varfetch(aia_name) & eui_map=scope_varfetch(eui_name)
;  per_spe_dl, aia_map, aia_index, aia_pixel_size, aia_cad, $
;              eui_map, eui_index, eui_pixel_size, eui_cad, 4, dl_region=dl_regions, crop_option=crop_option;,$
;;              SAVENAME='Review_slit_' + strtrim(i, 2)
;endfor
;;
;print, 'ANALYSIS DONE!'


;MAKING ALL PLOTS FOR PAPER-------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Final Figures'
;path='/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Final Figures'

;Fig. 1 - Setting the Scene  (.eps)
file171full = file_search(path+'/full_disc.fits')
fileHMI = file_search(path+'/hmi.fits')
file174full = file_search(path+'/solo_L2_eui-hrieuvopn-image_20241024T035009214_V01.fits')
help, file174full
print, file174full
read_sdo, file171full, index171full, data171full, /use_shared, /uncomp_delete
read_sdo, fileHMI, indexHMI, dataHMI, /use_shared, /uncomp_delete
mreadfits_tilecomp, file174full, index174full, data174full
;eui_readfits, file174full, index174full, data174full, quiet=quiet
hmi_pixel_size=mean(indexHMI.dsun_obs)/1e6 * mean(indexHMI.cdelt1) * !dpi/180./3600.
aia2hmi=aia_pixel_size/hmi_pixel_size
dataHMI=dataHMI[aia_crop_x1*aia2hmi:aia_crop_x2*aia2hmi, aia_crop_y1*aia2hmi:aia_crop_y2*aia2hmi]
data171=aia_data_c[*,*,0]
data174=eui_data_c[*,*,0]

left=0.07
right=0.05
bg=0.0
tg=0.05
mg=0.08
bottom=0.07
top=0.03
wx=(1-right-left-2*bg)/3.
bigwx=(1-left-right-tg)/2.

height=bottom+top+wx+bigwx+mg
wy=wx/height
bigwy=bigwx/height
;print, 'plot height = '+ strtrim(height, 2)


;checking with usual plotting
;set_plot, 'X'
;plot_width=1000
;wdef, 4, plot_width, plot_width*height
;!p.charsize=3
;!p.charthick=1
;!p.thick=2

;creating good quality plot for paper
set_plot,'ps'
!p.charthick=2.5
!p.charsize=1.5
!p.background=255
!p.font=-1
!p.multi=[0,1,5]
!p.color=0
;Use DEVICE to set PostScript device options
width_eps=17
height_eps=width_eps*height
device,filename='Setting the Scene.eps',$
  xsize=width_eps,ysize=height_eps,$
  xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$
  /color,bits=8,/Helvetica,isolatin1=1;,/bold

aia_pixel_size=0.433
eui_pixel_size=0.194


;Plot full sun 171
    aia_lct, wave=171, /load
    plot_image, (data171full)^0.4, xtickformat="(A1)", ytickformat="(A1)",$; xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]';, $
      position = [left, bottom+wy+mg, left+bigwx, bottom+wy+mg+bigwy], /nosquare, title='AIA 171 Full FoV'
;    xyouts, side+0.005, bottom+4*siz+gap+0.018, strtrim(index_1[0].date_obs,2), charsize=0.8, , color=255, /normal
    xyouts, left+0.01, bottom+wy+mg+bigwy-0.021, '(a) '+strtrim(index171full[0].date_obs,2), charsize=0.8, color=255, /normal

;Plot full HRIEUV
    plot_image, alog(data174full), ytickformat="(A1)", xtickformat="(A1)",$
      position = [left+bigwx+tg, bottom+wy+mg, 1-right, bottom+wy+mg+bigwy], /nosquare, title='HRIEUV Full FoV'
    xyouts, left+bigwx+tg+0.01, bottom+wy+mg+bigwy-0.021, '(b) '+strtrim(index174full[0].date_obs,2), charsize=0.8, color=255, /normal

;Plot HMI crop region
    loadct, 0, /silent
    plot_image, alog(dataHMI), xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]',$
      position = [left, bottom, left+wx, bottom+wy], scale=[hmi_pixel_size, hmi_pixel_size], $
      title='HMI Continuum ROI', xtickv=[0, 10, 20, 30, 40]
    xyouts, left+0.01,bottom+wy-0.021, '(c)', charsize=0.8, color=0, /normal
    for i = 0, 3 do begin
      aia_x=crop_aia_xs[i*5:(i+1)*5-1]
      aia_y=crop_aia_ys[i*5:(i+1)*5-1]
      aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
      aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
      aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
      aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
      loadct, 39, /silent
      slit_col=97
      oplot, aia_xtrack*aia_pixel_size, aia_ytrack*aia_pixel_size, color=slit_col, linestyle=0, thick=3
      oplot, aia_xtrack_m[0:-2]*aia_pixel_size, aia_ytrack_m[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, aia_xtrack_p[0:-2]*aia_pixel_size, aia_ytrack_p[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
      xyouts, (aia_x[-1]+1)*aia_pixel_size, (aia_y[-1]-0.5)*aia_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=0, charthick=3
    endfor


;Plot 171 crop region
    aia_lct, wave=171, /load
    plot_image, alog(data171), xtitle='Solar X [Mm]', ytickformat="(A1)",$
     position = [left+wx+bg, bottom, left+2*wx+bg, bottom+wy], scale=[aia_pixel_size, aia_pixel_size], $
      title='AIA 171 ROI', xtickv=[0, 10, 20, 30, 40]
    xyouts, left+wx+bg+0.01,bottom+wy-0.021, '(d)', charsize=0.8, color=255, /normal
    for i = 0, 3 do begin
      aia_x=crop_aia_xs[i*5:(i+1)*5-1]
      aia_y=crop_aia_ys[i*5:(i+1)*5-1]
      aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
      aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
      aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
      aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
      loadct, 39, /silent
      oplot, aia_xtrack*aia_pixel_size, aia_ytrack*aia_pixel_size, color=slit_col, linestyle=0, thick=3
      oplot, aia_xtrack_m[0:-2]*aia_pixel_size, aia_ytrack_m[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, aia_xtrack_p[0:-2]*aia_pixel_size, aia_ytrack_p[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
      xyouts, (aia_x[-1]+1)*aia_pixel_size, (aia_y[-1]-0.5)*aia_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=255, charthick=3
    endfor
    

;Plot 174 crop region
    aia_lct, wave=171, /load
    plot_image, alog(data174), xtitle='Solar X [Mm]', ytickformat="(A1)",$
      position = [left+2*wx+2*bg, bottom, left+3*wx+2*bg, bottom+wy], scale =[eui_pixel_size, eui_pixel_size], $
      title='HRIEUV 174 ROI', xtickname=[0, 10, 20, 30, 40]
    xyouts, left+2*wx+2*bg+0.01,bottom+wy-0.021, '(e)', charsize=0.8, color=255, /normal
    for i = 0, 3 do begin
      eui_x=crop_eui_xs[i*5:(i+1)*5-1]
      eui_y=crop_eui_ys[i*5:(i+1)*5-1]
      eui_tracks=make_tracks(eui_x, eui_y) & eui_xtrack=eui_tracks[*,0] & eui_ytrack=eui_tracks[*,1]
      eui_track_m1=perp_direc(eui_xtrack, eui_ytrack, -1) & eui_track_p1=perp_direc(eui_xtrack, eui_ytrack, +1)
      eui_xtrack_m1=(eui_track_m1[*,0]) & eui_ytrack_m1=(eui_track_m1[*,1])
      eui_xtrack_p1=(eui_track_p1[*,0]) & eui_ytrack_p1=(eui_track_p1[*,1])
      eui_track_m2=perp_direc(eui_xtrack, eui_ytrack, -2) & eui_track_p2=perp_direc(eui_xtrack, eui_ytrack, +2)
      eui_xtrack_m2=(eui_track_m2[*,0]) & eui_ytrack_m2=(eui_track_m2[*,1])
      eui_xtrack_p2=(eui_track_p2[*,0]) & eui_ytrack_p2=(eui_track_p2[*,1])
      
      eui_track_m3=perp_direc(eui_xtrack, eui_ytrack, -3) & eui_track_p3=perp_direc(eui_xtrack, eui_ytrack, +3)
      eui_xtrack_m3=(eui_track_m3[*,0]) & eui_ytrack_m3=(eui_track_m3[*,1])
      eui_xtrack_p3=(eui_track_p3[*,0]) & eui_ytrack_p3=(eui_track_p3[*,1])
      loadct, 7, /silent
      slit_col=140
      oplot, eui_xtrack*eui_pixel_size, eui_ytrack*eui_pixel_size, color=slit_col, linestyle=0, thick=3
      oplot, eui_xtrack_m1[0:-2]*eui_pixel_size, eui_ytrack_m1[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, eui_xtrack_p1[0:-2]*eui_pixel_size, eui_ytrack_p1[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, eui_xtrack_m2[0:-2]*eui_pixel_size, eui_ytrack_m2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, eui_xtrack_p2[0:-2]*eui_pixel_size, eui_ytrack_p2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, eui_xtrack_m3[0:-2]*eui_pixel_size, eui_ytrack_m2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      oplot, eui_xtrack_p3[0:-2]*eui_pixel_size, eui_ytrack_p2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
      xyouts, (eui_x[-1]+2)*eui_pixel_size, (eui_y[-1]-1)*eui_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=255, charthick=3
    endfor

    ;adding box for window expansion
    ;EUI FoV on AIA FoV 
    xx1 = 0.220 & yy1 = 0.632
    xx2 = 0.315 & yy2 = 0.741
    loadct, 7, /silent
    boxcol=140
    plots, [xx1, xx2, xx2, xx1, xx1], [yy1, yy1, yy2, yy2, yy1], color=boxcol, thick = 4, linestyle=0, /normal
    plots, [xx2, left+bigwx+tg], [yy2, bottom+wy+mg+bigwy], color=boxcol, thick = 4, linestyle=0, /normal
    plots, [xx2, left+bigwx+tg], [yy1, bottom+wy+mg], color=boxcol, thick = 4, linestyle=0, /normal
    
    ;HRIEUV
    xx1 = 0.776 & yy1 = 0.607
    xx2 = 0.834 & yy2 = 0.669
    plots, [xx1, xx2, xx2, xx1, xx1], [yy1, yy1, yy2, yy2, yy1], color=boxcol, thick = 4, linestyle=0, /normal
    plots, [xx1, left+2*wx], [yy1, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
    plots, [xx2, left+3*wx], [yy1, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
    
    
    ;HMI and 171
    loadct, 39, /silent
    boxcol=97
    xx1 = 0.270 & yy1 = 0.660
    xx2 = 0.285 & yy2 = 0.676
    plots, [xx1, xx2, xx2, xx1, xx1], [yy1, yy1, yy2, yy2, yy1], color=boxcol, thick = 4, linestyle=0, /normal
    plots, [xx1, left], [yy1, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx2, left+wx], [yy1, bottom+wy], color=100, thick = 4, linestyle=0, /normal
;    plots, [xx1, left+wx], [yy1, bottom+wy], color=100, thick = 4, linestyle=0, /normal
    plots, [xx2, left+2*wx], [yy1, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
    
    

device, /close
set_plot, 'X'

;;;Fig. 2 - TD Analysis (.eps)
;;;To use this section, need to comment out any plotting/window making with wdef and opening alternative windows
;;set_plot,'ps'
;;!p.multi=[0, 2, 6]
;;!p.charthick=2.5
;;!p.charsize=1.35
;;width_eps=20
;;height_eps=20
;;device,filename='TD Analysis.eps',$
;;  xsize=width_eps,ysize=height_eps,$
;;  xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$
;;  /color,bits=8,/Helvetica,isolatin1=1;,/bold
;; 
;;aia_lct, wave=171, /load
;;aia_use=aia_td_1
;;eui_use=eui_td_1
;;per_spe_dl, aia_use, aia_index_c, aia_pixel_size, 12., $
;;  eui_use, eui_index_c, eui_pixel_size, 3., 4, $
;;  SAVENAME='slit_' + strtrim(i, 2), $
;;  crop_option=1, dl_region = [5,25]
;;
;;label_xpos=0.08
;;alabel_ypos = 0.955
;;ydiff=0.167
;;xyouts, label_xpos, alabel_ypos, '(a)', /normal, charsize=1, charthick=3, color=255
;;xyouts, label_xpos, alabel_ypos-ydiff, '(b)', /normal, charsize=1, charthick=3, color=255
;;xyouts, label_xpos, alabel_ypos-2*ydiff, '(c)', /normal, charsize=1
;;xyouts, label_xpos, alabel_ypos-3*ydiff, '(d)', /normal, charsize=1
;;xyouts, label_xpos, alabel_ypos-4*ydiff, '(e)', /normal, charsize=1
;;xyouts, label_xpos, alabel_ypos-5*ydiff, '(f)', /normal, charsize=1
;;
;;device, /close
;;set_plot, 'X'

; Fig. 3 - Decay Lengths (.eps)
; This script produces an 8‑panel figure (AIA + EUI for 4 slits).
; Panels are arranged manually using normalized coordinates.
; The goal of this commented version is to make it easier to
; change panel spacing, positions, and overall layout.
; ================================

;; ----- GENERAL PLOTTING SETUP -----
;set_plot,'ps'
;;wdef, 13, 1000, 500
;!p.multi=[0, 2, 4]
;!p.charsize=1.3
;!p.color=0
;!p.background=255
;!p.charthick=2.5
;width_eps=20
;height_eps=9.5
;
;device,filename='Decay Length.eps',$ 
;   xsize=width_eps,ysize=height_eps,$ 
;   xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$ 
;  /color,bits=8,/Helvetica,isolatin1=1;  ,/bold
;
;; ----- PANEL LAYOUT CONTROL -----
;; You can change these four variables to adjust the panel spacing:
;; gap  = vertical spacing between rows (in normalized coordinates)
;; bgap = bottom margin
;; sgap = left margin
;; h    = panel height (computed from gap)
;; w    = panel width  (computed from sgap)
;;
;; To shift everything up/down or left/right:
;;   - increase/decrease bgap (vertical shift)
;;   - increase/decrease sgap (horizontal shift)
;;
;; To increase spacing between rows or columns:
;;   - increase gap  (vertical spacing)
;;   - change sgap or formula for w for horizontal spacing
;;
;; To make panels taller or wider:
;;   - modify the formulas for h and w directly.
;
;gap=0.06      ; spacing between panel rows
;cgap=0.02     ; NEW: horizontal gap between left and right panels      ; spacing between panel rows
;bgap=0.1      ; bottom margin
;sgap=0.1      ; left margin
;h=(1-4*gap-bgap)/4.   ; automatic panel height
;w=(1-2*sgap-cgap)/2.        ; automatic panel width (2 columns)
;
;
;; ----- DECAY LENGTH & PERIOD INPUTS -----
;aia_dls=[4.9, 3.7, 2.5, 2.2]
;aia_dl_errs=[0.4, 0.3, 0.1, 0.5]
;aia_dets=[8.7, 3.9, 8.2, 3.4]
;aia_periods=[2.5, 2.5, 3.0, 2.5]
;
;eui_dls=[3.7, 2.4, 5.8, 3.7]
;eui_dl_errs=[0.7, 0.2, 0.9, 0.3]
;eui_dets=[8.9, 3.5, 4.6, 2.9]
;eui_periods=[2.4, 2.4, 3.0, 2.4]
;
;
;; ----- MAIN PANEL LOOP -----
;for j =0, 3 do begin
;  ; Fetch maps by name
;  aia_name='aia_map_' + STRTRIM(j+1, 2)+'_crop'
;  eui_name='eui_map_'+strtrim(j+1,2)+'_crop'
;  aia_map=scope_varfetch(aia_name)
;  eui_map=scope_varfetch(eui_name)
;
;  ; Crop maps to 8 Mm vertically
;  aia_map=aia_map[*, 0:8./aia_pixel_size]
;  eui_map=eui_map[*, 0:8./eui_pixel_size]
;
;  ; ----- AXIS LABEL FORMAT -----
;  if j eq 3 then xtitle_str = 'Time [min]' else xtitle_str = ''
;  if j eq 3 then xtickformat = "(I0)" else xtickformat="(A1)"
;
;  ; ----- AIA PANEL -----
;  ; Compute dynamic panel position using sgap, gap, bgap, h, w
;  ; Row index is (3-j) so that slit 1 is row 4, slit 4 is row 1
;  aia_lct, wave=171, /load
;  plot_image, bright_spot(aia_map, 1500), /nosquare, $
;    xtitle=xtitle_str, ytitle='Distance [Mm]', $
;    title='AIA, Slit '+strtrim(j+1,2), $
;    scale=[aia_cad/60., aia_pixel_size], $
;    position=[sgap, bgap+gap*(3-j)+(3-j)*h, sgap+w, bgap+gap*(3-j)+(4-j)*h], $
;    xtickformat=xtickformat
;
;  ; Plot decay length line & error box
;  loadct, 39, /silent
;  oplot, [0, n_elements(aia_map[*,0])], [aia_dls[j],aia_dls[j]], color=255, linestyle=2, thick=4
;  x1=0 & x2=n_elements(aia_map[*,0])
;  y1=aia_dls[j]-aia_dl_errs[j] & y2=aia_dls[j]+aia_dl_errs[j]
;  oplot, [x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], color=255, thick=1.5
;
;
;  ; ----- EUI PANEL -----
;  aia_lct, wave=171, /load
;  plot_image, bright_spot(eui_map, 1000), /nosquare, $
;    xtitle=xtitle_str, $
;    title='EUI, Slit '+strtrim(j+1,2), $
;    scale=[eui_cad/60., eui_pixel_size], $
;    position=[sgap+w+cgap, bgap+gap*(3-j)+(3-j)*h, sgap+2*w+cgap, bgap+gap*(3-j)+(4-j)*h], $
;    ytickformat="(A1)", xtickformat=xtickformat
;
;  ; Plot decay length and error box for EUI
;  loadct, 39, /silent
;  oplot, [0, n_elements(eui_map[*,0])], [eui_dls[j],eui_dls[j]], color=255, linestyle=2, thick=4
;  x1=0 & x2=n_elements(eui_map[*,0])
;  y1=eui_dls[j]-eui_dl_errs[j] & y2=eui_dls[j]+eui_dl_errs[j]
;  oplot, [x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], color=255, thick=1.5
;endfor
;
;
;; ----- PANEL LABELS (a,i), (a,ii), ... -----
;xposi=0.105     ; left‑column label horizontal position
;xposii=0.515    ; right‑column label horizontal position
;ypos=0.9        ; top label position
;ydiff=0.225     ; vertical spacing between label rows
;
;labels=['(a,i)', '(a,ii)', '(b,i)', '(b,ii)', '(c,i)', '(c,ii)', '(d,i)', '(d,ii)']
;
;for X = 0, 3 do begin
;  xyouts, xposi, ypos-X*ydiff, labels[2*X], color=255, charthick=3, /normal, charsize=0.9
;  xyouts, xposii, ypos-X*ydiff, labels[2*X+1], color=255, charthick=3, /normal, charsize=0.9
;endfor
;
;device, /close
;set_plot, 'X'


END


;
;
;END