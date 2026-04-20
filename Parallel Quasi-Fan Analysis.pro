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

function nearest_index, tarr, t
  ; returns index of tarr closest to t
  return, (where(abs(tarr - t) eq min(abs(tarr - t)), cnt))[0]
end

cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Correct Parallel'

;; Load pre-selected AIA data and index
;AIA_index = aia_index
;AIA_data  = aia_data
;EUI_index = eui_index_spacealigned
;EUI_data  = eui_data_spacealigned_corrected
;
;aia_cad=12.
;eui_cad=5.
;
;; Compute AIA and EUI pixel sizes (in radians per pixel)
;aia_pixel_size = mean(aia_index.dsun_obs)/1e6 * mean(aia_index.cdelt1) * !dpi/180./3600.
;eui_pixel_size = mean(eui_index.dsun_obs)/1e6 * mean(eui_index.cdelt1) * !dpi/180./3600.

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
;print,'AIA length [s] = '+ strtrim(n_elements(aia_data[0,0,*])*aia_cad)
;print,'EUI length [s] = '+ strtrim(n_elements(eui_data[0,0,*])*eui_cad)
;
;;Light travel time + aligning data---------------------------------------------------
;AIA_dist_m = mean(AIA_index.dsun_obs)
;EUI_dist_m = mean(EUI_index.dsun_obs)
;ddist = AIA_dist_m-EUI_dist_m
;dtime = ddist/(3.0e+8)
;AIA_dframes =  round(dtime/aia_cad)
;EUI_dframes =  round(dtime/eui_cad)
;print, ' '
;print, '---Correcting for light travel time---'
;print, 'Travel time = '+ strtrim(dtime)
;print, 'AIA Frames to crop = '+strtrim(AIA_dframes)
;print, 'EUI Frames to crop = '+strtrim(EUI_dframes)
;
;;-2 from start of EUI-> same length in time
;;-6 from start of AIA for same observing
;;-14 from end of EUI
;
;aia_data_ready = aia_data[*,*, AIA_dframes:-1]
;eui_data_ready = eui_data[*,*, 2:-EUI_dframes]
;
;save, filename='aia_data_ready.sav', aia_data_ready
;save, filename='eui_data_ready.sav', eui_data_ready


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
;aia_data_final = aia_data_cropped
;eui_data_final = eui_data_cropped
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Correct Parallel/MOVIE'
;;aia_image=comprange(aia_data_a)
;;eui_image=comprange(eui_data_a)
;;cube_difference, aia_data_a, aia_data_f, 'pfilter'
;;cube_difference, eui_data_a, eui_data_f, 'pfilter'
;;cube_difference, aia_data_c, aia_data_f, 'pfilter', cad=aia_cad
;;cube_difference, eui_data_c, eui_data_f, 'pfilter', cad=eui_cad
;;aia_image=aia_data_f
;;eui_image=eui_data_f
;
;aia_lct, wave=171, /load
;n_aia=n_elements(aia_data_final[0,0,*])
;n_eui=n_elements(eui_data_final[0,0,*])
;t_aia = findgen(n_aia)*aia_cad
;t_eui = findgen(n_eui)*eui_cad
;dt_movie=1.0
;t0=0.0
;t1 = min([t_aia[n_aia-1], t_eui[n_eui-1]])
;n_movie = fix((t1 - t0)/dt_movie) + 1
;t_mov = t0 + findgen(n_movie)*dt_movie
;
;for k=0, n_movie-1 do begin
;  ia = nearest_index(t_aia, t_mov[k])
;  ie = nearest_index(t_eui, t_mov[k])
;
;  aia_image = aia_data_final[*,*,ia]
;  eui_image = eui_data_final[*,*,ie]
;  
;  aia_plot=comprange(aia_image)
;  plot_image, bytscl(aia_plot, 0.2, 0.8), scale=[aia_pixel_size, aia_pixel_size], $
;    title='AIA RoI', xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]'
;  eui_plot=comprange(eui_image)
;  plot_image, bytscl(eui_plot, 0.2, 0.8), scale=[eui_pixel_size, eui_pixel_size], $
;    title='EUI RoI', xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]'
;  write_jpeg, 'AIA_and_EUI'+strtrim(k, 2)+'.jpg', TVRD(/TRUE), /TRUE, quality=100
;endfor
;
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Correct Parallel'

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


;;Choosing positions, select just in AIA and then convert to EUI scale
;wdef, 0, 1400, 800
;!p.multi=[0, 2, 1]
;!p.background=255
;!p.color=0
;!p.charsize=1
;aia_lct, wave=171, /load
;plot_image, aia_data_final[*,*,0], xtitle='Solar X [pix]', ytitle='Solar Y [pix]', title = 'AIA - Select 3 lots of 5 points for td maps'
;fixed_aia_xs = []
;fixed_aia_ys = []
;for i = 0 , 2 do begin
;  print, 'Slit '+strtrim(i+1,2)+': select 5 points'
;  for j=0, 4 do begin
;    cursor, x, y, /down
;    oplot, [x], [y], color=255, psym=4, thick=2, symsize=1.5
;    fixed_aia_xs = [fixed_aia_xs, x]
;    fixed_aia_ys = [fixed_aia_ys, y]
;  endfor
;endfor
;fixed_eui_xs = fixed_aia_xs*aia_pixel_size/eui_pixel_size
;fixed_eui_ys = fixed_aia_ys*aia_pixel_size/eui_pixel_size
;plot_image, eui_data_final[*,*,0], xtitle='Solar X [pix]', ytitle='Solar Y [pix]', title = 'EUI - Points converted to EUI'
;oplot, fixed_eui_xs, fixed_eui_ys, color=255, psym=4, thick=2, symsize=1.5


;create td maps 
;; DEFINING CURVED SLIT POSITIONS
fixed_aia_xs=[116.41183, 125.28643, 133.38932, 143.03562, 153.45363, $
              116.41183, 122.19961, 127.21569, 133.38932, 139.94881, $
              120.10957, 123.22211, 127.50185, 131.00346, 135.28320]


fixed_aia_ys=[122.28644, 122.28644, 122.28644, 122.67230, 126.53082, $
              124.21571, 126.14497, 130.00349, 134.24786, 139.26394, $
              127.50346, 131.78321, 136.45202, 141.12083, 146.56777]

              
fixed_eui_xs=[168.79228, 183.12990, 198.57041, 211.25369, 220.62829, $
              165.48360, 174.85820, 182.57845, 191.40160, 200.77620, $
              167.49162, 173.61059, 180.28583, 186.40480, 195.30512]


fixed_eui_ys=[173.21834, 172.66689, 172.66689, 175.97557, 180.93859, $
              172.66689, 177.62991, 183.69583, 190.31319, 198.58489, $
              176.51091, 183.18615, 189.30512, 196.53663, 210.44338]


;save, fixed_aia_xs, fixed_aia_ys, fixed_eui_xs, fixed_eui_ys, filename='fixed_slit_positions.sav'
;restore, '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Final Figures/fixed_slit_positions.sav'


;;CROPPING TO SMALLER FOV-
;;-----------------------------------------------------------
;;Choosing AIA crop region in pixels
;aia_width=130
;aia_crop_x1=50
;aia_crop_x2=aia_crop_x1+aia_width
;aia_crop_y1=80
;aia_crop_y2=aia_crop_y1+aia_width
;aia_data_cropped=aia_data_a[aia_crop_x1:aia_crop_x2, aia_crop_y1:aia_crop_y2,*]
;
;ratio=aia_pixel_size/eui_pixel_size
;
;eui_width=(aia_crop_x2-aia_crop_x1)*ratio
;eui_crop_x1=aia_crop_x1*ratio
;eui_crop_x2=eui_crop_x1+eui_width
;eui_crop_y1=aia_crop_y1*ratio
;eui_crop_y2=eui_crop_y1+eui_width
;eui_data_cropped=eui_data_a[eui_crop_x1:eui_crop_x2, eui_crop_y1:eui_crop_y2,*]
;
;crop_aia_xs=fixed_aia_xs-aia_crop_x1
;crop_aia_ys=fixed_aia_ys-aia_crop_y1
;crop_eui_xs=fixed_eui_xs-eui_crop_x1
;crop_eui_ys=fixed_eui_ys-eui_crop_y1

;;Dont need to run below ***...*** part, just restore 'all_td_maps.sav' after first time
;!p.background=255
;!p.color=0
;!p.charthick=1
;aia_lct, wave=171, /load
;;***
;for i = 0 , 2 do begin
;  ;create td map
;  xpoints=fixed_aia_xs[i*5: i*5+4]
;  ypoints=fixed_aia_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, aia_data_final, int_slice, sm_scl=250./aia_cad, ncut=1, width=3, winnum=0, points=points
;  aia_td=transpose(int_slice)
;  xpoints=fixed_eui_xs[i*5: i*5+4]
;  ypoints=fixed_eui_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, eui_data_final, int_slice, sm_scl=250./eui_cad, ncut=1, width=5, winnum=0, points=points
;  eui_td=transpose(int_slice)
;  save, filename='aia_eui_td_maps_'+strtrim(i+1,2)+'.sav', aia_td, eui_td
;endfor
;;restore and save all td map info
;path='/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Correct Parallel'
;restore, 'aia_eui_td_maps_1.sav'
;aia_td_1=[aia_td] & eui_td_1=[eui_td]
;restore, 'aia_eui_td_maps_2.sav'
;aia_td_2=[aia_td] & eui_td_2=[eui_td]
;restore, 'aia_eui_td_maps_3.sav'
;aia_td_3=[aia_td] & eui_td_3=[eui_td]
;save, filename='all_td_maps.sav', aia_td_1, aia_td_2, aia_td_3, $
;      eui_td_1, eui_td_2, eui_td_3
;;***

;;;-----------------------QUICK RUN ANALYSIS HERE!!!!!!-----------------------------------------
;;restore, 'all_td_maps.sav'
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response/Parallel'
;crop_options=[0,0,0]
;;dl_regions=[5,25, 20,40, 5,25, 0,20, 10,30 ]
;dl_regions=[0,10, 0,10, 0,10]
;aia_pixel_size=0.431 & aia_cad=12.
;eui_pixel_size=0.304 & eui_cad=5.
;
;for i = 1, 3 do begin
;  crop_option=crop_options[i-1]
;  aia_name='event2_aia_td_' + STRTRIM(i, 2) & eui_name='event2_eui_td_'+strtrim(i,2)
;  aia_map=scope_varfetch(aia_name) & eui_map=scope_varfetch(eui_name)
;  per_spe_dl, aia_map, aia_index, aia_pixel_size, aia_cad, $
;              eui_map, eui_index, eui_pixel_size, eui_cad, 4, dl_region=dl_regions, crop_option=crop_option,$
;              SAVENAME='Review_slit_' + strtrim(i, 2)
;endfor
;;
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
;crop_aia_xs=fixed_aia_xs-aia_crop_x1
;crop_aia_ys=fixed_aia_ys-aia_crop_y1
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
;crop_aia_xs=fixed_aia_xs-aia_crop_x1
;crop_aia_ys=fixed_aia_ys-aia_crop_y1
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
;plot_image, alog(aia_data_cropped[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
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
;plot_image, alog(eui_data_cropped[*,*,0]), xtitle='Solar X [pix]', ytitle='Solar Y [pix]', $
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



;MAKING ALL PLOTS FOR PAPER-------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response/Parallel'
;path='/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan/Final Figures'

;;Fig. 1 - Setting the Scene  (.eps)
;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response/Parallel'
;path = '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Correct Parallel/Figures'
;file171full = file_search(path+'/full_disc.fits')
;fileHMI = file_search(path+'/hmi.fits')
;file174full = file_search(path+'/full_hrieuv.fits')
;;help, file174full
;;print, file174full
;read_sdo, file171full, index171full, data171full, /use_shared, /uncomp_delete
;read_sdo, fileHMI, indexHMI, dataHMI, /use_shared, /uncomp_delete
;mreadfits_tilecomp, file174full, index174full, data174full
;;eui_readfits, file174full, index174full, data174full, quiet=quiet
;hmi_pixel_size=mean(indexHMI.dsun_obs)/1e6 * mean(indexHMI.cdelt1) * !dpi/180./3600.
;aia2hmi=aia_pixel_size/hmi_pixel_size
;crop_x1 = 40 & crop_y1 = 35 
;crop_x2 = 80 & crop_y2 = 75
;hmi_crop_x1 = 47 & hmi_crop_y1 = 44
;hmi_crop_x2 = 87 & hmi_crop_y2 = 84
;dataHMI=dataHMI[hmi_crop_x1/hmi_pixel_size: hmi_crop_x2/hmi_pixel_size, hmi_crop_y1/hmi_pixel_size: hmi_crop_y2/hmi_pixel_size]
;data171=aia_data_cropped[*,*,0];[crop_x1/aia_pixel_size: crop_x2/aia_pixel_size, crop_y1/aia_pixel_size: crop_y2/aia_pixel_size, 0]
;data174=eui_data_cropped[*,*,0];[crop_x1/eui_pixel_size: crop_x2/eui_pixel_size, crop_y1/eui_pixel_size: crop_y2/eui_pixel_size, 0]
;
;left=0.07
;right=0.05
;bg=0.0
;tg=0.05
;mg=0.08
;bottom=0.07
;top=0.03
;wx=(1-right-left-2*bg)/3.
;bigwx=(1-left-right-tg)/2.
;
;height=bottom+top+wx+bigwx+mg
;wy=wx/height
;bigwy=bigwx/height
;;print, 'plot height = '+ strtrim(height, 2)
;
;
;;checking with usual plotting
;;set_plot, 'X'
;plot_width=1000
;;wdef, 4, plot_width, plot_width*height
;;!p.charsize=3
;;!p.charthick=1
;;!p.thick=2
;
;;creating good quality plot for paper
;set_plot,'ps'
;!p.charthick=2.5
;!p.charsize=1.5
;!p.background=255
;!p.font=-1
;!p.multi=[0,1,5]
;!p.color=0
;;Use DEVICE to set PostScript device options
;width_eps=17
;height_eps=width_eps*height
;device,filename='Setting the Scene.eps',$
;  xsize=width_eps,ysize=height_eps,$
;  xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$
;  /color,bits=8,/Helvetica,isolatin1=1;,/bold
;
;
;
;;Plot full sun 171
;    aia_lct, wave=171, /load
;    plot_image, (data171full)^0.4, xtickformat="(A1)", ytickformat="(A1)",$; xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]';, $
;      position = [left, bottom+wy+mg, left+bigwx, bottom+wy+mg+bigwy], /nosquare, title='AIA 171 Full FoV'
;;    xyouts, side+0.005, bottom+4*siz+gap+0.018, strtrim(index_1[0].date_obs,2), charsize=0.8, , color=255, /normal
;    xyouts, left+0.01, bottom+wy+mg+bigwy-0.021, '(a) '+strtrim(index171full[0].date_obs,2), charsize=0.8, color=255, /normal
;
;;Plot full HRIEUV
;    plot_image, alog(data174full), ytickformat="(A1)", xtickformat="(A1)",$
;      position = [left+bigwx+tg, bottom+wy+mg, 1-right, bottom+wy+mg+bigwy], /nosquare, title='HRIEUV Full FoV'
;    xyouts, left+bigwx+tg+0.01, bottom+wy+mg+bigwy-0.021, '(b) '+strtrim(index174full[0].date_obs,2), charsize=0.8, color=255, /normal
;
;;Plot HMI crop region
;    loadct, 0, /silent
;    plot_image, alog(dataHMI), xtitle='Solar X [Mm]', ytitle='Solar Y [Mm]',$
;      position = [left, bottom, left+wx, bottom+wy], scale=[hmi_pixel_size, hmi_pixel_size], $
;      title='HMI Continuum ROI', xticks= 4, xtickv=[0, 10, 20, 30]
;    xyouts, left+0.01,bottom+wy-0.021, '(c)', charsize=0.8, color=0, /normal
;    for i = 0, 2 do begin
;      aia_x=fixed_aia_xs[i*5:(i+1)*5-1]-crop_x1/aia_pixel_size
;      aia_y=fixed_aia_ys[i*5:(i+1)*5-1]-crop_y1/aia_pixel_size
;      aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
;      aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
;      aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
;      aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
;      loadct, 39, /silent
;      slit_col=97
;      oplot, aia_xtrack*aia_pixel_size, aia_ytrack*aia_pixel_size, color=slit_col, linestyle=0, thick=3
;      oplot, aia_xtrack_m[0:-2]*aia_pixel_size, aia_ytrack_m[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
;      oplot, aia_xtrack_p[0:-2]*aia_pixel_size, aia_ytrack_p[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
;      xyouts, (aia_x[-1]+1)*aia_pixel_size, (aia_y[-1]-0.5)*aia_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=0, charthick=3
;    endfor
;
;
;;Plot 171 crop region
;    aia_lct, wave=171, /load
;    plot_image, alog(data171), xtitle='Solar X [Mm]', ytickformat="(A1)",$
;     position = [left+wx+bg, bottom, left+2*wx+bg, bottom+wy], scale=[aia_pixel_size, aia_pixel_size], $
;      title='AIA 171 ROI', xticks= 4, xtickv=[0, 10, 20, 30]
;    xyouts, left+wx+bg+0.01,bottom+wy-0.021, '(d)', charsize=0.8, color=255, /normal
;    for i = 0, 2 do begin
;      aia_x=fixed_aia_xs[i*5:(i+1)*5-1]-crop_x1/aia_pixel_size
;      aia_y=fixed_aia_ys[i*5:(i+1)*5-1]-crop_y1/aia_pixel_size
;      aia_tracks=make_tracks(aia_x, aia_y) & aia_xtrack=aia_tracks[*,0] & aia_ytrack=aia_tracks[*,1]
;      aia_track_m=perp_direc(aia_xtrack, aia_ytrack, -1) & aia_track_p=perp_direc(aia_xtrack, aia_ytrack, +1)
;      aia_xtrack_m=aia_track_m[*,0] & aia_ytrack_m=aia_track_m[*,1]
;      aia_xtrack_p=aia_track_p[*,0] & aia_ytrack_p=aia_track_p[*,1]
;      loadct, 39, /silent
;      oplot, aia_xtrack*aia_pixel_size, aia_ytrack*aia_pixel_size, color=slit_col, linestyle=0, thick=3
;      oplot, aia_xtrack_m[0:-2]*aia_pixel_size, aia_ytrack_m[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
;      oplot, aia_xtrack_p[0:-2]*aia_pixel_size, aia_ytrack_p[0:-2]*aia_pixel_size, color=slit_col, linestyle=2, thick=3
;      xyouts, (aia_x[-1]+1)*aia_pixel_size, (aia_y[-1]-0.5)*aia_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=255, charthick=3
;    endfor
;    
;
;;Plot 174 crop region
;    aia_lct, wave=171, /load
;    plot_image, alog(data174), xtitle='Solar X [Mm]', ytickformat="(A1)",$
;      position = [left+2*wx+2*bg, bottom, left+3*wx+2*bg, bottom+wy], scale =[eui_pixel_size, eui_pixel_size], $
;      title='HRIEUV 174 ROI', xticks= 5, xtickv=[0, 10, 20, 30, 40]
;    xyouts, left+2*wx+2*bg+0.01,bottom+wy-0.021, '(e)', charsize=0.8, color=255, /normal
;    for i = 0, 2 do begin
;      eui_x=fixed_eui_xs[i*5:(i+1)*5-1]-crop_x1/eui_pixel_size
;      eui_y=fixed_eui_ys[i*5:(i+1)*5-1]-crop_y1/eui_pixel_size
;      eui_tracks=make_tracks(eui_x, eui_y) & eui_xtrack=eui_tracks[*,0] & eui_ytrack=eui_tracks[*,1]
;      eui_track_m1=perp_direc(eui_xtrack, eui_ytrack, -1) & eui_track_p1=perp_direc(eui_xtrack, eui_ytrack, +1)
;      eui_xtrack_m1=(eui_track_m1[*,0]) & eui_ytrack_m1=(eui_track_m1[*,1])
;      eui_xtrack_p1=(eui_track_p1[*,0]) & eui_ytrack_p1=(eui_track_p1[*,1])
;      eui_track_m2=perp_direc(eui_xtrack, eui_ytrack, -2) & eui_track_p2=perp_direc(eui_xtrack, eui_ytrack, +2)
;      eui_xtrack_m2=(eui_track_m2[*,0]) & eui_ytrack_m2=(eui_track_m2[*,1])
;      eui_xtrack_p2=(eui_track_p2[*,0]) & eui_ytrack_p2=(eui_track_p2[*,1])
;      loadct, 7, /silent
;      slit_col=140
;      oplot, eui_xtrack*eui_pixel_size, eui_ytrack*eui_pixel_size, color=slit_col, linestyle=0, thick=3
;      oplot, eui_xtrack_m1[0:-2]*eui_pixel_size, eui_ytrack_m1[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
;      oplot, eui_xtrack_p1[0:-2]*eui_pixel_size, eui_ytrack_p1[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
;      oplot, eui_xtrack_m2[0:-2]*eui_pixel_size, eui_ytrack_m2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
;      oplot, eui_xtrack_p2[0:-2]*eui_pixel_size, eui_ytrack_p2[0:-2]*eui_pixel_size, color=slit_col, linestyle=2, thick=3
;      xyouts, (eui_x[-1]+2)*eui_pixel_size, (eui_y[-1]-1)*eui_pixel_size, '('+strtrim(i+1,2)+')', charsize=1, color=255, charthick=3
;    endfor
;
;    ;adding box for window expansion
;    ;EUI FoV on AIA FoV 
;    xx1 = 0.227 & yy1 = 0.843
;    xx2 = 0.363 & yy2 = 0.783
;    xx3 = 0.3097 & yy3 = 0.6299043152
;    xx4 = 0.1737 & yy4 = 0.6899043152
;
;    loadct, 7, /silent
;    boxcol=140
;    plots, [xx1, xx2, xx3, xx4, xx1], [yy1, yy2, yy3, yy4, yy1], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx2, left+bigwx+tg], [yy2, bottom+wy+mg+bigwy], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx3, left+bigwx+tg], [yy3, bottom+wy+mg], color=boxcol, thick = 4, linestyle=0, /normal
;    
;    ;HRIEUV
;;    xx1 = 0.631 & yy1 = 0.73
;;    xx2 = 0.6532 & yy2 = 0.705
;    xx1 = 0.63581561 & yy1 = 0.71630213
;    xx2 = 0.65648494 & yy2 = 0.72542095
;    xx3 = 0.64838439 & yy3 = 0.74868849
;    xx4 = 0.62771506 & yy4 = 0.73956967
;    plots, [xx1, xx2, xx3, xx4, xx1], [yy1, yy2, yy3, yy4, yy1], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx1, left+2*wx], [yy1, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx2, left+3*wx], [yy2, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
;    
;    
;    ;HMI and 171
;    loadct, 39, /silent
;    boxcol=97
;    xx1 = 0.231 & yy1 = 0.761
;    xx2 = 0.24165 & yy2 = 0.749
;    plots, [xx1, xx2, xx2, xx1, xx1], [yy1, yy1, yy2, yy2, yy1], color=boxcol, thick = 4, linestyle=0, /normal
;    plots, [xx1, left], [yy2, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
;;    plots, [xx2, left+wx], [yy1, bottom+wy], color=100, thick = 4, linestyle=0, /normal
;;    plots, [xx1, left+wx], [yy1, bottom+wy], color=100, thick = 4, linestyle=0, /normal
;    plots, [xx2, left+2*wx], [yy2, bottom+wy], color=boxcol, thick = 4, linestyle=0, /normal
;    
;    
;
;device, /close
;set_plot, 'X'

;;Fig. 2 - TD Analysis (.eps)
;;To use this section, need to comment out any plotting/window making with wdef and opening alternative windows
;set_plot,'ps'
;!p.multi=[0, 2, 6]
;!p.charthick=2.5
;!p.charsize=1.35
;width_eps=20
;height_eps=20
;device,filename='TD Analysis.eps',$
;  xsize=width_eps,ysize=height_eps,$
;  xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$
;  /color,bits=8,/Helvetica,isolatin1=1;,/bold
; 
;aia_lct, wave=171, /load
;aia_use=aia_td_1
;eui_use=eui_td_1
;per_spe_dl, aia_use, aia_index_c, aia_pixel_size, 12., $
;  eui_use, eui_index_c, eui_pixel_size, 3., 4, $
;  SAVENAME='slit_' + strtrim(i, 2), $
;  crop_option=1, dl_region = [5,25]
;
;label_xpos=0.08
;alabel_ypos = 0.955
;ydiff=0.167
;xyouts, label_xpos, alabel_ypos, '(a)', /normal, charsize=1, charthick=3, color=255
;xyouts, label_xpos, alabel_ypos-ydiff, '(b)', /normal, charsize=1, charthick=3, color=255
;xyouts, label_xpos, alabel_ypos-2*ydiff, '(c)', /normal, charsize=1
;xyouts, label_xpos, alabel_ypos-3*ydiff, '(d)', /normal, charsize=1
;xyouts, label_xpos, alabel_ypos-4*ydiff, '(e)', /normal, charsize=1
;xyouts, label_xpos, alabel_ypos-5*ydiff, '(f)', /normal, charsize=1
;
;device, /close
;set_plot, 'X'

;; Fig. 3 - Decay Lengths (.eps)
;; This script produces an 8‑panel figure (AIA + EUI for 4 slits).
;; Panels are arranged manually using normalized coordinates.
;; The goal of this commented version is to make it easier to
;; change panel spacing, positions, and overall layout.
;; ================================
;;event1_aia_td_1 = aia_td_1
;;event1_aia_td_2 = aia_td_2
;;event1_aia_td_3 = aia_td_3
;;event1_aia_td_4 = aia_td_5
;;event1_eui_td_1 = eui_td_1
;;event1_eui_td_2 = eui_td_2
;;event1_eui_td_3 = eui_td_3
;;event1_eui_td_4 = eui_td_5

; ----- GENERAL PLOTTING SETUP -----
set_plot,'ps'
!p.multi=[0, 2, 7]      ; 2 columns, 7 rows
!p.charsize=1.4
!p.color=0
!p.background=255
!p.charthick=2.5

width_eps=20
height_eps=19.0        ; ~2x original height (was 9.5)

device,filename='Decay Length Parallel NEW.eps',$ 
   xsize=width_eps,ysize=height_eps,$ 
   xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$ 
  /color,bits=8,/Helvetica,isolatin1=1

;wdef, 0, 600, 700


; ----- PANEL LAYOUT CONTROL -----
gap=0.025           ; regular spacing between adjacent rows
group_gap=0.07     ; extra gap between Event 1 (top 4) and Event 2 (bottom 3)

cgap=0.02          ; horizontal gap between left and right panels
bgap=0.08          ; bottom margin
tgap=0.06          ; top margin
sgap=0.10          ; left margin

nrows_total  = 7
nrows_event2 = 3
nrows_event1 = 4

; panel height from available vertical space
avail = 1.0 - bgap - tgap - (nrows_total-1)*gap - group_gap
h = avail / float(nrows_total)

; panel width for 2 columns
w = (1.0 - 2*sgap - cgap) / 2.0


; ----- DECAY LENGTH & PERIOD INPUTS -----
; EVENT 1 (4 slits) 
event1_aia_dls     = [4.5, 3.8, 2.2, 1.5]
event1_aia_dl_errs = [0.3, 0.6, 0.1, 0.1]
event1_eui_dls     = [2.6, 3.6, 5.0, 3.0]
event1_eui_dl_errs = [0.3, 0.5, 1.0, 0.2]

; EVENT 2 (3 slits) 
event2_aia_dls     = [2.1, 2.1, 1.7]
event2_aia_dl_errs = [0.1, 0.1, 0.1]
event2_eui_dls     = [2.1, 1.8, 1.5]
event2_eui_dl_errs = [0.1, 0.1, 0.1]


; ============================================================
; ======================= EVENT 1 (TOP) ======================
; rows 6..3 (top to bottom), slits 1..4
; ============================================================
eui_cad = 3.
aia_pixel_size=0.433
eui_pixel_size=0.194
for j = 0, 3 do begin
  print, 'j = '+strtrim(j)
  ; Row index in the 7-row stack:
  ; j=0 (slit1) -> row=6 (top)
  ; j=3 (slit4) -> row=3
  row = (nrows_total-1) - j

  ; y0 calculation with group_gap applied for rows >= nrows_event2
  ; (Event 1 rows are 3..6, so they all get the group_gap offset)
  y0 = bgap + row*(h + gap) + group_gap
  y1 = y0 + h

  ; Fetch TD maps by name (Event 1) - NO "_crop"
  aia_name = 'event1_aia_td_' + strtrim(j+1,2)
  print, aia_name
  eui_name = 'event1_eui_td_' + strtrim(j+1,2)
  aia_map = scope_varfetch(aia_name)
  eui_map = scope_varfetch(eui_name)

  ; Crop TD maps to 8 Mm vertically
  aia_crops = [1.5, 1.0, 0, 0]
  eui_crops = [0, 0, 3.2, 0]
  aia_map = aia_map[*, aia_crops[j]:aia_crops[j]+8./aia_pixel_size]
  eui_map = eui_map[*, eui_crops[j]:eui_crops[j]+8./eui_pixel_size]
  
  ; ------------------------------------------------------------
  ; DETREND using SMOOTH over 250 s with mirrored edges
  ; (smooth along TIME axis, assumed to be the FIRST dimension)
  ; ------------------------------------------------------------
  win_s = 250.0
  
  ; AIA: window in samples
  w_aia = long(round(win_s / aia_cad))
  if w_aia lt 3 then w_aia = 3
  if (w_aia mod 2) eq 0 then w_aia = w_aia + 1   ; force odd
  
  ; EUI: window in samples
  w_eui = long(round(win_s / eui_cad))
  if w_eui lt 3 then w_eui = 3
  if (w_eui mod 2) eq 0 then w_eui = w_eui + 1   ; force odd
  
  ; Smooth along TIME (dimension 1 in SMOOTH's box size)
  aia_bg  = smooth(aia_map, [w_aia, 1], /edge_mirror)
  eui_bg  = smooth(eui_map, [w_eui, 1], /edge_mirror)
  
  aia_map = aia_map - aia_bg
  eui_map = eui_map - eui_bg


  ; Time axis on bottom row of Event 1 (slit 4, j=3)
  if j eq 3 then xtitle_str = 'Time [min]' else xtitle_str = ''
  if j eq 3 then xtickformat = "(I0)" else xtickformat="(A1)"

  ; ----- AIA PANEL (Event 1) -----
  aia_lct, wave=171, /load
  plot_image, (aia_map), /nosquare, $
    xtitle=xtitle_str, ytitle='Distance [Mm]', $
    title='AIA, Slit '+strtrim(j+1,2), $
    scale=[aia_cad/60., aia_pixel_size], $
    position=[sgap, y0, sgap+w, y1], $
    xtickformat=xtickformat
  if j eq 0 then let = 'a'
  if j eq 1 then let = 'b'
  if j eq 2 then let = 'c'
  if j eq 3 then let = 'd'
  xyouts, sgap+0.007, y1-0.02, '('+let+',i)', color=255, /normal, charsize=0.9, charthick=3

  loadct, 39, /silent
  oplot, [0, n_elements(aia_map[*,0])], [event1_aia_dls[j], event1_aia_dls[j]], $
         color=255, linestyle=2, thick=4
  x1=0 & x2=n_elements(aia_map[*,0])
  yb=event1_aia_dls[j]-event1_aia_dl_errs[j]
  yt=event1_aia_dls[j]+event1_aia_dl_errs[j]
  oplot, [x1, x2, x2, x1, x1], [yb, yb, yt, yt, yb], color=255, thick=1.5
  

  ; ----- EUI PANEL (Event 1) -----
  aia_lct, wave=171, /load
  plot_image, bright_spot(eui_map, 1600), /nosquare, $
    xtitle=xtitle_str, $
    title='EUI, Slit '+strtrim(j+1,2), $
    scale=[eui_cad/60., eui_pixel_size], $
    position=[sgap+w+cgap, y0, sgap+2*w+cgap, y1], $
    ytickformat="(A1)", xtickformat=xtickformat
  xyouts, sgap+w+cgap+0.007, y1-0.02, '('+let+',ii)', color=255, /normal, charsize=0.9, charthick=3

  loadct, 39, /silent
  oplot, [0, n_elements(eui_map[*,0])], [event1_eui_dls[j], event1_eui_dls[j]], $
         color=255, linestyle=2, thick=4
  x1=0 & x2=n_elements(eui_map[*,0])
  yb=event1_eui_dls[j]-event1_eui_dl_errs[j]
  yt=event1_eui_dls[j]+event1_eui_dl_errs[j]
  oplot, [x1, x2, x2, x1, x1], [yb, yb, yt, yt, yb], color=255, thick=1.5

endfor


; ============================================================
; ====================== EVENT 2 (BOTTOM) ====================
; rows 2..0 (top to bottom), slits 1..3
; ============================================================
eui_cad = 5.
aia_pixel_size=0.431
eui_pixel_size=0.304
for j = 0, 2 do begin

  ; Row index for Event 2 group:
  ; j=0 (slit1) -> row=2 (top of bottom group)
  ; j=2 (slit3) -> row=0 (bottom of figure)
  row = (nrows_event2-1) - j

  ; Event 2 rows do NOT get the group_gap offset
  y0 = bgap + row*(h + gap)
  y1 = y0 + h

  ; Fetch TD maps by name (Event 2) - NO "_crop"
  aia_name = 'event2_aia_td_' + strtrim(j+1,2)
  eui_name = 'event2_eui_td_' + strtrim(j+1,2)
  aia_map = scope_varfetch(aia_name)
  eui_map = scope_varfetch(eui_name)

  ; Crop TD maps to 8 Mm vertically
  aia_crops = [1.9, 1.0, 0] + [0,0.92,0]
  eui_crops = [0, 0, 0.9] + [0,1.52,0]
  aia_map = aia_map[*, aia_crops[j]:aia_crops[j]+8./aia_pixel_size]
  eui_map = eui_map[*, eui_crops[j]:eui_crops[j]+8./eui_pixel_size]
  
  ; ------------------------------------------------------------
  ; DETREND using SMOOTH over 250 s with mirrored edges
  ; (smooth along TIME axis, assumed to be the FIRST dimension)
  ; ------------------------------------------------------------
  win_s = 250.0
  
  ; AIA: window in samples
  w_aia = long(round(win_s / aia_cad))
  if w_aia lt 3 then w_aia = 3
  if (w_aia mod 2) eq 0 then w_aia = w_aia + 1   ; force odd
  
  ; EUI: window in samples
  w_eui = long(round(win_s / eui_cad))
  if w_eui lt 3 then w_eui = 3
  if (w_eui mod 2) eq 0 then w_eui = w_eui + 1   ; force odd
  
  ; Smooth along TIME (dimension 1 in SMOOTH's box size)
  aia_bg  = smooth(aia_map, [w_aia, 1], /edge_mirror)
  eui_bg  = smooth(eui_map, [w_eui, 1], /edge_mirror)
  
  aia_map = aia_map - aia_bg
  eui_map = eui_map - eui_bg


  ; Time axis on bottom row of Event 2 (slit 3, j=2)
  if j eq 2 then xtitle_str = 'Time [min]' else xtitle_str = ''
  if j eq 2 then xtickformat = "(I0)" else xtickformat="(A1)"

  ; ----- AIA PANEL (Event 2) -----
  aia_lct, wave=171, /load
  plot_image, (aia_map), /nosquare, $
    xtitle=xtitle_str, ytitle='Distance [Mm]', $
    title='AIA, Slit '+strtrim(j+1,2), $
    scale=[aia_cad/60., aia_pixel_size], $
    position=[sgap, y0, sgap+w, y1], $
    xtickformat=xtickformat
    
  if j eq 0 then let = 'e'
  if j eq 1 then let = 'f'
  if j eq 2 then let = 'g'
  xyouts, sgap+0.007, y1-0.02, '('+let+',i)', color=255, /normal, charsize=0.9, charthick=3

  loadct, 39, /silent
  oplot, [0, n_elements(aia_map[*,0])], [event2_aia_dls[j], event2_aia_dls[j]], $
         color=255, linestyle=2, thick=4
  x1=0 & x2=n_elements(aia_map[*,0])
  yb=event2_aia_dls[j]-event2_aia_dl_errs[j]
  yt=event2_aia_dls[j]+event2_aia_dl_errs[j]
  oplot, [x1, x2, x2, x1, x1], [yb, yb, yt, yt, yb], color=255, thick=1.5


  ; ----- EUI PANEL (Event 2) -----
  aia_lct, wave=171, /load
  plot_image, (eui_map), /nosquare, $
    xtitle=xtitle_str, $
    title='EUI, Slit '+strtrim(j+1,2), $
    scale=[eui_cad/60., eui_pixel_size], $
    position=[sgap+w+cgap, y0, sgap+2*w+cgap, y1], $
    ytickformat="(A1)", xtickformat=xtickformat
  xyouts, sgap+w+cgap+0.007, y1-0.02, '('+let+',ii)', color=255, /normal, charsize=0.9, charthick=3

  loadct, 39, /silent
  oplot, [0, n_elements(eui_map[*,0])], [event2_eui_dls[j], event2_eui_dls[j]], $
         color=255, linestyle=2, thick=4
  x1=0 & x2=n_elements(eui_map[*,0])
  yb=event2_eui_dls[j]-event2_eui_dl_errs[j]
  yt=event2_eui_dls[j]+event2_eui_dl_errs[j]
  oplot, [x1, x2, x2, x1, x1], [yb, yb, yt, yt, yb], color=255, thick=1.5

endfor


; ----- EVENT HEADERS (centered over both columns) -----
xmid = sgap + (2*w + cgap)/2.0

; Event 1 header: above top row (row=6) which has y1 = bgap + 6*(h+gap) + group_gap + h
y_event1_top = bgap + (nrows_total-1)*(h+gap) + group_gap + h
xyouts, xmid, y_event1_top + 0.025, 'Event 1 - 2024-10-24', /normal, align=0.5, $
       charsize=1.1, charthick=3, color=0

; Event 2 header: above row=2 -> y1 = bgap + 2*(h+gap) + h
y_event2_top = bgap + (nrows_event2-1)*(h+gap) + h
xyouts, xmid, y_event2_top + 0.025, 'Event 2 - 2021-11-05', /normal, align=0.5, $
       charsize=1.1, charthick=3, color=0



device, /close
set_plot, 'X'

END



;
;
;END