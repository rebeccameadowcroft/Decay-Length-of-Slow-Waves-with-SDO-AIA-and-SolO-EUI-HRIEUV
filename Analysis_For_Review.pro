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

function interpol_td_map, aia_map_in, aia_cad, aia_pixel_size, eui_map_in, eui_cad, eui_pixel_size
  ;interpolate for time
  nT_a = n_elements(aia_map_in[*,0])
  nT_e = n_elements(eui_map_in[*,0])
  tA = findgen(nT_a) * aia_cad
  tE = findgen(nT_e) * eui_cad
  ; Use overlapping time range only
  tmax_common = min([max(tA), max(tE)])

  ; Use the finer cadence as the common cadence
  dt_new = min([aia_cad, eui_cad])
  nTnew = floor(tmax_common / dt_new) + 1
  tNew = findgen(nTnew) * dt_new

  aia_inter = fltarr(nTnew, n_elements(aia_map_in[0,*]))
  eui_inter = fltarr(nTnew, n_elements(eui_map_in[0,*]))

  for h = 0, n_elements(aia_map_in[0,*]) - 1 do begin
    aia_inter[*,h] = interpol(aia_map_in[*,h], tA, tNew)
  endfor
  for h = 0, n_elements(eui_map_in[0,*]) - 1 do begin
    eui_inter[*,h] = interpol(eui_map_in[*,h], tE, tNew)
  endfor

  ;interpolate for distance
  nHa = n_elements(aia_map_in[0,*])
  nHe = n_elements(eui_map_in[0,*])
  hA = findgen(nHa) * aia_pixel_size
  hE = findgen(nHe) * eui_pixel_size

  ; Use overlapping height range only
  Hmax_common = min([max(hA), max(hE)])

  ; Use the finer spatial sampling as the common height spacing
  dh = min([aia_pixel_size, eui_pixel_size])
  nHnew = floor(Hmax_common / dh) + 1
  hNew = findgen(nHnew) * dh

  aia_hires = fltarr(nTnew, nHnew)
  eui_hires = fltarr(nTnew, nHnew)

  for t = 0, nTnew - 1 do begin
    aia_hires[t,*] = interpol(aia_inter[t,*], hA, hNew)
    eui_hires[t,*] = interpol(eui_inter[t,*], hE, hNew)
  endfor
  return, {aia:aia_hires, eui:eui_hires, cad:dt_new, pix:dh}
end

;cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response Parallel'
;
;;Redoing TD maps with event 1 aia 7 pixels wide and interpolating to HRIEUV resolution-----------------------------------------
;;1) Load datas and indexs
;path = '/Users/rebecca/Documents/PhD Year 3/Project 5 - Quasi-Fan'
;aia_data=aia_data_cropped
;aia_index=aia_index
;eui_data=eui_data_cropped
;eui_index=eui_index_spacealigned
;
;aia_cad=12.
;eui_cad=5.
;
;; Compute AIA and EUI pixel sizes (in radians per pixel)
;aia_pixel_size = mean(aia_index.dsun_obs)/1e6 * mean(aia_index.cdelt1) * !dpi/180./3600.
;eui_pixel_size = mean(eui_index.dsun_obs)/1e6 * mean(eui_index.cdelt1) * !dpi/180./3600.
;;;EVENT 1
;;aia_xs=[113.113, 116.875, 119.508, 121.013, 120.637,$
;;  115.74627, 121.01315, 126.65624, 131.17071, 137.94241,$
;;  125.15141, 130.79450, 136.06138, 141.32826, 145.84273,$
;;  127.40865, 134.55656, 140.57585, 148.85238, 156.75270]
;;
;;aia_ys=[113.489, 116.875, 122.142, 130.042, 138.695, $
;;  113.48904, 117.62730, 124.39901, 130.41830, 139.82344,$
;;  92.797725, 80.759141, 69.849175, 58.186797, 43.138568,$
;;  96.559782, 89.411873, 83.768787, 73.987438, 64.206089]
;;
;;eui_xs=[165.829, 179.368, 190.327, 200.642, 207.734,$
;;  169.05240, 178.72282, 189.03793, 198.70835, 209.02346,$
;;  170.98649, 183.23568, 198.06366, 212.89163, 226.43021,$
;;  183.88038, 199.35305, 214.18102, 230.94308, 248.34983]
;;
;;eui_ys=[194.840, 201.932, 207.734, 218.049, 232.232,$
;;  194.19550, 194.84019, 198.06366, 203.22122, 212.24694,$
;;  160.67139, 145.19872, 125.85789, 106.51705, 85.242138,$
;;  156.80322, 144.55403, 130.37075, 112.31930, 96.201943]
;;
;;;crop_aia_xs=fixed_aia_xs-aia_crop_x1
;;;crop_aia_ys=fixed_aia_ys-aia_crop_y1
;;;crop_eui_xs=fixed_eui_xs-eui_crop_x1
;;;crop_eui_ys=fixed_eui_ys-eui_crop_y1
;
;;EVENT 2
;aia_xs=[116.41183, 125.28643, 133.38932, 143.03562, 153.45363, $
;  116.41183, 122.19961, 127.21569, 133.38932, 139.94881, $
;  120.10957, 123.22211, 127.50185, 131.00346, 135.28320]
;
;aia_ys=[122.28644, 122.28644, 122.28644, 122.67230, 126.53082, $
;  124.21571, 126.14497, 130.00349, 134.24786, 139.26394, $
;  127.50346, 131.78321, 136.45202, 141.12083, 146.56777]
;
;eui_xs=[168.79228, 183.12990, 198.57041, 211.25369, 220.62829, $
;  165.48360, 174.85820, 182.57845, 191.40160, 200.77620, $
;  167.49162, 173.61059, 180.28583, 186.40480, 195.30512]
;
;eui_ys=[173.21834, 172.66689, 172.66689, 175.97557, 180.93859, $
;  172.66689, 177.62991, 183.69583, 190.31319, 198.58489, $
;  176.51091, 183.18615, 189.30512, 196.53663, 210.44338]
;
;!p.background=255
;!p.color=0
;!p.charthick=1
;aia_lct, wave=171, /load
;for i = 0 , 2 do begin
;  ;create td map
;  xpoints=aia_xs[i*5: i*5+4]
;  ypoints=aia_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, aia_data, int_slice, sm_scl=250./aia_cad, ncut=1, width=3, winnum=0, points=points
;  aia_td=transpose(int_slice)
;  xpoints=eui_xs[i*5: i*5+4]
;  ypoints=eui_ys[i*5: i*5+4]
;  points=[xpoints, ypoints]
;  stplot_curve, eui_data, int_slice, sm_scl=250./eui_cad, ncut=1, width=7, winnum=0, points=points
;  eui_td=transpose(int_slice)
;  save, filename='aia_eui_td_maps_'+strtrim(i+1,2)+'.sav', aia_td, eui_td
;endfor
;
;;restore and save all td map info
;restore, 'aia_eui_td_maps_1.sav'
;aia_td_1=[aia_td] & eui_td_1=[eui_td]
;restore, 'aia_eui_td_maps_2.sav'
;aia_td_2=[aia_td] & eui_td_2=[eui_td]
;restore, 'aia_eui_td_maps_3.sav'
;aia_td_3=[aia_td] & eui_td_3=[eui_td]
;restore, 'aia_eui_td_maps_4.sav'
;aia_td_4=[aia_td] & eui_td_4=[eui_td]
;
;wdef, 10, 800, 800
;!p.multi=[0, 2, 4]
;!p.background = 255
;!p.color=0
;for j = 1, 3 do begin
;  aia_name='aia_td_' + STRTRIM(j, 2) & eui_name='eui_td_'+strtrim(j,2)
;  aia_map=scope_varfetch(aia_name) & eui_map=scope_varfetch(eui_name)
;  plot_image, aia_map, title='AIA_Map_'+strtrim(j,2), xtitle='Time [frames]', ytitle='Distance [Mm]', /nosquare
;  plot_image, eui_map, title='EUI_Map_'+strtrim(j,2), xtitle='Time [frames]', ytitle='Distance [Mm]', /nosquare
;endfor
;write_jpeg, 'All TD Maps 7 EUI.jpg',TVRD(/TRUE), /TRUE, quality=100
;
;save, filename='all_tdog_maps 7 EUI.sav', aia_td_1, aia_td_2, aia_td_3, aia_td_4, $
;      eui_td_1, eui_td_2, eui_td_3, eui_td_4
;
;;restore, 'all_td_maps.sav'
;crop_options=[1,1,2]
;dl_regions=[5,25, 20,40, 5,25, 10,30 ]
;
;for j = 1, 3 do begin
;  aia_name='aia_td_' + STRTRIM(j, 2) & eui_name='eui_td_'+strtrim(j,2)
;  aia_map=scope_varfetch(aia_name) & eui_map=scope_varfetch(eui_name)
;  per_spe_dl, aia_map, aia_index, aia_pixel_size, aia_cad, $
;              eui_map, eui_index, eui_pixel_size, eui_cad, 4, $
;              SAVENAME='slit_' + strtrim(j, 2)+'_7_EUI', $
;              crop_option=crop_options[j-1], dl_region = [0,10]
;;  save, aia_map, eui_map, filename='Slit_'+strtrim(j,2)+'_maps.sav'
;endfor
;
;print, 'ANALYSIS DONE!'

;Interpolating all TD plots to HRIEUV and saving
cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Review Response'


aia_cad=12.
aia_pixel_size=0.432

eui_cad1=3.
eui_cad2=5.
eui_pixel_size1=0.194
eui_pixel_size2=0.304


;tmp=interpol_td_map(event1_aia_td_1, aia_cad, aia_pixel_size, event1_eui_td_1, eui_cad1, eui_pixel_size1)
;event1_aia_hires_1 = tmp.aia & event1_eui_hires_1 = tmp.eui
;tmp=interpol_td_map(event1_aia_td_2, aia_cad, aia_pixel_size, event1_eui_td_2, eui_cad1, eui_pixel_size1)
;event1_aia_hires_2 = tmp.aia & event1_eui_hires_2 = tmp.eui
;tmp=interpol_td_map(event1_aia_td_3, aia_cad, aia_pixel_size, event1_eui_td_3, eui_cad1, eui_pixel_size1)
;event1_aia_hires_3 = tmp.aia & event1_eui_hires_3 = tmp.eui
;tmp=interpol_td_map(event1_aia_td_4, aia_cad, aia_pixel_size, event1_eui_td_4, eui_cad1, eui_pixel_size1)
;event1_aia_hires_4 = tmp.aia & event1_eui_hires_4 = tmp.eui
;
;tmp=interpol_td_map(event2_aia_td_1, aia_cad, aia_pixel_size, event2_eui_td_1, eui_cad2, eui_pixel_size2)
;event2_aia_hires_1 = tmp.aia & event2_eui_hires_1 = tmp.eui
;tmp=interpol_td_map(event2_aia_td_2, aia_cad, aia_pixel_size, event2_eui_td_2, eui_cad2, eui_pixel_size2)
;event2_aia_hires_2 = tmp.aia & event2_eui_hires_2 = tmp.eui
;tmp=interpol_td_map(event2_aia_td_3, aia_cad, aia_pixel_size, event2_eui_td_3, eui_cad2, eui_pixel_size2)
;event2_aia_hires_3 = tmp.aia & event2_eui_hires_3 = tmp.eui
;
;save, filename='interpolated_aia_td_maps.sav', event1_aia_hires_1, event1_aia_hires_2, event1_aia_hires_3, event1_aia_hires_4,$
;  event2_aia_hires_1, event2_aia_hires_2, event2_aia_hires_3

; ----- GENERAL PLOTTING SETUP -----
set_plot,'ps'
!p.multi=[0, 2, 7]      ; 2 columns, 7 rows
!p.charsize=1.4
!p.color=0
!p.background=255
!p.charthick=2.5

width_eps=19
height_eps=19.0        ; ~2x original height (was 9.5)

device,filename='DecayLengthParallelNEW.eps',$
  xsize=width_eps,ysize=height_eps,$
  xoff=(21.0-width_eps)/2.0,yoff=(29.7-height_eps)/2.0,$
  /color,bits=24,/Helvetica, encapsulated=1; ,isolatin1=1

;wdef, 0, 600, 700


; ----- PANEL LAYOUT CONTROL -----
gap=0.025           ; regular spacing between adjacent rows
group_gap=0.07     ; extra gap between Event 1 (top 4) and Event 2 (bottom 3)

cgap=0.02          ; horizontal gap between left and right panels
bgap=0.08          ; bottom margin
tgap=0.06          ; top margin
sgap=0.05          ; left margin
rgap=0.02          ; right margin

nrows_total  = 7
nrows_event2 = 3
nrows_event1 = 4

; panel height from available vertical space
avail = 1.0 - bgap - tgap - (nrows_total-1)*gap - group_gap
h = avail / float(nrows_total)

; panel width for 2 columns
w = (1.0 - sgap - rgap - cgap) / 2.0


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
aia_cad=eui_cad
eui_pixel_size=0.194
aia_pixel_size=eui_pixel_size
for j = 0, 3 do begin
;  print, 'j = '+strtrim(j)
  ; Row index in the 7-row stack:
  ; j=0 (slit1) -> row=6 (top)
  ; j=3 (slit4) -> row=3
  row = (nrows_total-1) - j

  ; y0 calculation with group_gap applied for rows >= nrows_event2
  ; (Event 1 rows are 3..6, so they all get the group_gap offset)
  y0 = bgap + row*(h + gap) + group_gap
  y1 = y0 + h

  ; Fetch TD maps by name (Event 1) - NO "_crop"
  aia_name = 'event1_aia_hires_' + strtrim(j+1,2)
  print, aia_name
  eui_name = 'event1_eui_hires_' + strtrim(j+1,2)
  aia_map = scope_varfetch(aia_name)
  eui_map = scope_varfetch(eui_name)

  ; Crop TD maps to 8 Mm vertically
  aia_crops = [1.5, 1.0, 0, 0]
  eui_crops = [0, 0, 3.2, 0]
  print, size(aia_map)
  print, aia_crops[j]/aia_pixel_size+8./aia_pixel_size
  aia_map = aia_map[*, aia_crops[j]/aia_pixel_size:aia_crops[j]/aia_pixel_size+7./aia_pixel_size]
  eui_map = eui_map[*, eui_crops[j]/eui_pixel_size:eui_crops[j]/eui_pixel_size+7./eui_pixel_size]

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
  plot_image, (aia_map), $
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
  load, 0, /silent
  oplot, [5, 25, 25, 5, 5], [0, 0, 7, 7, 0], color=200, thick=4


  ; ----- EUI PANEL (Event 1) -----
  aia_lct, wave=171, /load
  plot_image, bright_spot(eui_map, 1600), $
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
  load, 0, /silent
  oplot, [5, 25, 25, 5, 5], [0, 0, 7, 7, 0], color=200, thick=4

endfor


; ============================================================
; ====================== EVENT 2 (BOTTOM) ====================
; rows 2..0 (top to bottom), slits 1..3
; ============================================================
eui_cad = 5.
aia_cad=eui_cad
eui_pixel_size=0.304
aia_pixel_size=eui_pixel_size
for j = 0, 2 do begin

  ; Row index for Event 2 group:
  ; j=0 (slit1) -> row=2 (top of bottom group)
  ; j=2 (slit3) -> row=0 (bottom of figure)
  row = (nrows_event2-1) - j

  ; Event 2 rows do NOT get the group_gap offset
  y0 = bgap + row*(h + gap)
  y1 = y0 + h

  ; Fetch TD maps by name (Event 2) - NO "_crop"
  aia_name = 'event2_aia_hires_' + strtrim(j+1,2)
  eui_name = 'event2_eui_hires_' + strtrim(j+1,2)
  aia_map = scope_varfetch(aia_name)
  eui_map = scope_varfetch(eui_name)

  ; Crop TD maps to 8 Mm vertically
  aia_crops = [1.9, 1.0, 0]
  eui_crops = [0, 0, 0.9]
  aia_map = aia_map[*, aia_crops[j]/aia_pixel_size:aia_crops[j]/aia_pixel_size+7./aia_pixel_size]
  eui_map = eui_map[*, eui_crops[j]/eui_pixel_size:eui_crops[j]/eui_pixel_size+7./eui_pixel_size]

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
  plot_image, (aia_map), $
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
  load, 0, /silent
  oplot, [0, 10, 10, 0, 0], [0, 0, 7, 7, 0], color=200, thick=4


  ; ----- EUI PANEL (Event 2) -----
  aia_lct, wave=171, /load
  plot_image, (eui_map), $
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
  load, 0, /silent
  oplot, [0, 10, 10, 0, 0], [0, 0, 7, 7, 0], color=200, thick=4
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
;per_spe_dl, event1_aia_hires_1, aia_index, eui_pixel_size1, 3., $
;  event1_eui_hires_1, eui_index, eui_pixel_size1, 3., 4, $
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


END