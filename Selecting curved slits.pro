;Slow wave properties, parallel line of sight
;
cd, '/Users/rebecca/Documents/PhD Year 3/Project 6 - Decay Length /Paralle LOS'
;restore, 
aia_cad=12.
eui_cad=5.
aia_pixel_size = mean(aia_index.dsun_obs)/1.e6 * mean(aia_index.cdelt1) * !dpi/180./3600.
eui_pixel_size = mean(eui_index.dsun_obs)/1.e6 * mean(eui_index.cdelt1) * !dpi/180./3600.

; Slit selection - curved slits
;-------------------------------------
; Prepare storage
fixed_aia_xs = [] & fixed_aia_ys = []
fixed_eui_xs = [] & fixed_eui_ys = []
aia_td_1=[] & eui_td_1=[] & aia_xs_1=[] & aia_ys_1=[] & eui_xs_1=[] & eui_ys_1=[]
aia_td_2=[] & eui_td_2=[] & aia_xs_2=[] & aia_ys_2=[] & eui_xs_2=[] & eui_ys_2=[]
aia_td_3=[] & eui_td_3=[] & aia_xs_3=[] & aia_ys_3=[] & eui_xs_3=[] & eui_ys_3=[]
aia_td_4=[] & eui_td_4=[] & aia_xs_4=[] & aia_ys_4=[] & eui_xs_4=[] & eui_ys_4=[]
aia_td_5=[] & eui_td_5=[] & aia_xs_5=[] & aia_ys_5=[] & eui_xs_5=[] & eui_ys_5=[]


for i = 0, 2 do begin
  done = 0
  while done eq 0 do begin
    ;-----------------------------
    ; AIA Slit Selection
    ;-----------------------------
    stplot_curve, aia_data_final, sm_scl=250./aia_cad, ncut=1, width=3, winnum=0
    restore, filename='stplot_cut=0.sav'
    restore, filename='stplot_cutslocation.sav'
    aia_xs=xi
    aia_ys=yi
    aia_td=transpose(int_slice)

    ;-----------------------------
    ; EUI Slit Selection
    ;-----------------------------
    stplot_curve, eui_data_final, sm_scl=250./eui_cad, ncut=1, width=5, winnum=2
    restore, filename='stplot_cut=0.sav'
    restore, filename='stplot_cutslocation.sav'
    eui_xs=xi
    eui_ys=yi
    eui_td=transpose(int_slice)

    ;-----------------------------
    ; Perform Analysis
    ;-----------------------------
    per_spe_dl, aia_td, aia_index, aia_pixel_size, aia_cad, $
      eui_td, eui_index_a, eui_pixel_size, eui_cad, 4, $
      SAVENAME='slit_' + strtrim(i+1, 2)

    ;-----------------------------
    ; Check if user is happy
    ;-----------------------------
    print, 'Right-click to confirm this slit. Left-click to reselect.'
    cursor, junkx, junky, /data
    if !mouse.button eq 4 then begin
      ; Right-click = accept
      done = 1
      fixed_aia_xs=[fixed_aia_xs, aia_xs] & fixed_aia_ys=[fixed_aia_ys, aia_ys]
      fixed_eui_xs=[fixed_eui_xs, eui_xs] & fixed_eui_ys=[fixed_eui_ys, eui_ys]
    endif else begin
      print, 'Reselecting slit ', i+1
    endelse
  endwhile
endfor


END