; NAME:
;   Figure_2.pro
;
; PURPOSE:
;   Code used to produce Figure 2 from the Brewin et al. (In review) paper submitted to JGR-Oceans in 2022
;
; CATEGORY:
;   Data plotting
;
; CALLING SEQUENCE:
;   Figure_2
;
; INPUTS:
;   BGC-Argo float data file 6901573_Sprof.nc
;   Model parameters from fits conducted in Python (see main Python code for the paper)
;
; OUTPUTS:
;   Figure 2 to the paper
;
; MODIFICATION HISTORY:
; Bob Brewin 17/02/2022

pro Figure_2

  ;;;Float number
  FLOAT_NUMBER = '6901573'
  ;;;Define path to float data (NB: Needs to be edited depending on where data is)
  PATH_TO_FLOAT = '/Users/bobbrewin/Jupyter/Ideas/Two_comp_vertical_Chl/'
  ;;;Read the BGC-Argo netCDF file into IDL
  FILE_SPE = PATH_TO_FLOAT+FLOAT_NUMBER+'_Sprof.nc'
  
  ;;;SET PLOT UP
  SET_PLOT, 'PS'
  DEVICE, /color, /ENCAPSULATED, FILENAME = PATH_TO_FLOAT+'Figure_1.eps', XSIZE = 40, YSIZE = 50
  !p.font=0
  !X.MARGIN = [10, 3]
  !Y.MARGIN = [5, 2]
  device, /helvetica ;
  loadct, 39
  PSIZE=1.0
  PLOTSYM, 0, PSIZE, /FILL, THICK=1
  tau     = GREEK('tau')
  multiplot,[0,3,4,0,1],  ygap = 0.03, xgap = 0.035, /doyaxis, /doxaxis

  ;;;READ IN DATA FROM FLOAT
  filename_ig  = ncdf_open(FILE_SPE)
  ;;CORE ARGO
  Rrsid        = ncdf_varid(filename_ig,'PRES')
  ncdf_varget,filename_ig,Rrsid,PRES
  ;;;BASIC CONVERSION OF PRESSURE TO DEPTH
  DEPTH_ADJUSTED = ((((-1.82e-15  * PRES + 2.279e-10 ) * PRES - 2.2512e-5 ) * PRES + 9.72659) * PRES) / 9.81
  ;;BGC ARGO
  Rrsid        = ncdf_varid(filename_ig,'CHLA_ADJUSTED')
  ncdf_varget,filename_ig,Rrsid,CHLA_ADJUSTED
  Rrsid        = ncdf_varid(filename_ig,'BBP700_ADJUSTED')
  ncdf_varget,filename_ig,Rrsid,BBP700_ADJUSTED
  
  ;;;EXTRACT DATA FOR PROFILE 167 (NOTE IDL INDEX STARTS ON 0, SO INDEX 166)
  CHL = CHLA_ADJUSTED[*,166]
  AS = where(CHL gt 99990)
  CHL(AS) = !VALUES.F_NAN
  BBP = BBP700_ADJUSTED[*,166]
  AS = where(CHL gt 99990)
  BBP(AS) = !VALUES.F_NAN
  DEPTH  = DEPTH_ADJUSTED[*,166]
  AS = where(DEPTH gt 0 and DEPTH le 300)
  CHL   = CHL(AS)
  BBP   = BBP(AS)
  DEPTH = DEPTH(AS)
  
  ;;;MODEL PARAMATERS FROM FIT TO PROFILE 167 (Extracted from main Python code from the paper)
  Kd_PAR          = 0.048281593118789795
  Zp              = 4.6/Kd_PAR
  MLD             = 139.72
  CHL_SAT         = 0.17640002268665242
  BBP_SAT         = 0.000514139246661216
  P1              = 54.303707753866966
  T1              = 6.691624793773988
  Bm2             = 0
  T2              = 0
  Sig             = 0
  bbpk            = 0.5290324727671847
  bbpS1           = 1 - bbpk
  bbpS2           = 0
  
  ;;;NORMALISE DATA
  CHL_NOM         = CHL / CHL_SAT
  BBP_NOM         = BBP / BBP_SAT
  DEPTH_NOM       = DEPTH * Kd_PAR
     
  ;;;SHADING FOR PLOT
  CHL_TEST       = findgen(2000)*0.00025
  Table_data_SAT = fltarr([n_elements(CHL_TEST),500])
  Table_data_BEL = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
    minS   = 0.0
    maxS   = 1./Kd_PAR
    minB   = 1./Kd_PAR
    maxB   = 200
    Table_data_SAT(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
    Table_data_BEL(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor
  CHL_TEST       = findgen(2000)/1000
  Table_data_SATN = fltarr([n_elements(CHL_TEST),500])
  Table_data_BELN = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
    minS   = 0.0
    maxS   = 1.
    minB   = 1.
    maxB   = 12
    Table_data_SATN(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
    Table_data_BELN(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor      
  
  ;;;PLOTS PROFILE 167
  plot,  CHL_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,2.0], yrange=[10,0], symsize = 1.5, psym=5, $
       xminor = 4, /xstyle, /ystyle , /nodata, title = '2017-01-04', xtitle = '!8B!3!N / !8B!Ds!3!N', ytitle = ''+tau+'!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,2.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,2.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  DCM_pop = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
  POP_ADD = MLD_pop+DCM_pop
  loadct, 39
  oplot, DCM_pop, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, POP_ADD, DEPTH_NOM, psym = 0, color = 255, thick = 7
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.8,7.0]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.8,7.65]
  LEGEND, ['1+2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[255], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[7], position=[0.8,8.3]
  TVLCT, [[0],[100],[0]], 150
  LEGEND, ['!8Z!Dm!3!N','!8Z!Dp!3!N'], /fill,textcolors = 0, psym=[0,0], linestyle = [1,5], $
    usersym=usersym,colors=[150,200], BOX=0, spacing =1.5, charsize =1.5,CHARTHICK=2, $
    thick =[8,8], position=[0.0,2.0]
  legend, ['(a)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[1.1,0]
  AXIS, 2.0,0.0,XAXIS = 1, XRANGE=[0.0,2.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,2.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 2.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doyaxis, /doxaxis
  plot,  BBP_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,2.0], yrange=[10,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N / !8b!Dbp,s!3!N', ytitle = ''+tau+'!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,2.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,2.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop_bbp  =  MLD_pop*bbpS1
  DCM_pop_bbp  =  DCM_pop*bbpS2
  bbpk_pop_bbp = fltarr(n_elements(BBP_NOM))+bbpk
  POP_ADD_bbp = MLD_pop_bbp+DCM_pop_bbp+bbpk_pop_bbp
  loadct, 39
  oplot, DCM_pop_bbp, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp, DEPTH_NOM, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp, DEPTH_NOM, psym = 0, color = 255, thick = 7
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.8,7.0]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.8,7.65]
  LEGEND, ['!N!8b!S!E*!R!Ibp,k!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[20], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.8,8.3]
  LEGEND, ['1+2+!N!8b!S!E*!R!Ibp,k!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[255], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[7], position=[0.8,8.95]
  legend, ['(b)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[1.1,0]
  AXIS, 2.0,0.0,XAXIS = 1, XRANGE=[0.0,2.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,2.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 2.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doyaxis, /doxaxis
  plot, CHL,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.4], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8B!3!N !3(mg m!E-3!N!3)', ytitle = 'Depth (m)';, $
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3  
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.4],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.4],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop*CHL_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop*CHL_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, POP_ADD*CHL_SAT, DEPTH, psym = 0, color = 255, thick = 7
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.168,7.0/Kd_PAR]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.168,7.65/Kd_PAR]
  LEGEND, ['1+2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[255], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[7], position=[0.168,8.3/Kd_PAR]
  legend, ['(c)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.22,0]
  TVLCT, [[0],[100],[0]], 150
  AXIS, 0.4,0.0,XAXIS = 1, XRANGE=[0.0,0.4], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.4], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.4,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doyaxis, /doxaxis
  plot, BBP,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.001], yrange=[200,0], symsize = 1.5, psym=5, $
     xticks = 2, xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)', ytitle = 'Depth (m)';;;;, YTickFormat='exponent';, $
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.4],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.4],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp*BBP_SAT, DEPTH, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp*BBP_SAT, DEPTH, psym = 0, color = 255, thick = 7 
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.00042,7.0/Kd_PAR-5]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.00042,7.65/Kd_PAR-5]
  LEGEND, ['!N!8b!S!E k!R!Ibp!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[20], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.00042,8.3/Kd_PAR-5]
  LEGEND, ['1+2+!N!8b!S!E k!R!Ibp!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[255], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[7], position=[0.00042,8.95/Kd_PAR-5]
  legend, ['(d)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.00055,0]
  AXIS, 0.001,0.0,XAXIS = 1, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.001,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
   
  ;;;EXTRACT DATA FOR PROFILE 59 (NOTE IDL INDEX STARTS ON 0, SO INDEX 58)
  CHL = CHLA_ADJUSTED[*,58]
  AS = where(CHL gt 99990)
  CHL(AS) = !VALUES.F_NAN
  BBP = BBP700_ADJUSTED[*,58]
  AS = where(CHL gt 99990)
  BBP(AS) = !VALUES.F_NAN
  DEPTH  = DEPTH_ADJUSTED[*,58]
  AS = where(DEPTH gt 0 and DEPTH le 300)
  CHL   = CHL(AS)
  BBP   = BBP(AS)
  DEPTH = DEPTH(AS)
   
  ;;;MODEL PARAMATERS FROM FIT TO PROFILE 59 (Extracted from main Python code from the paper)
  Kd_PAR          = 0.060450627910307086
  Zp              = 4.6/Kd_PAR
  MLD             = 12.7
  CHL_SAT         = 0.24839996036730314
  BBP_SAT         = 0.0008170124492608011
  P1              = 7.434495388967195
  T1              = 2.7606623046621968
  Bm2             = 1.6479212523630853
  T2              = 5.259882060893519
  Sig             = 1.5258290179184437
  bbpk            = 0.39143617476027615
  bbpS1           = 1 - bbpk
  bbpS2           = 0.46502291661728246
  CHL_NOM         = CHL / CHL_SAT
  BBP_NOM         = BBP / BBP_SAT
  DEPTH_NOM       = DEPTH * Kd_PAR
   
  ;;;SHADING FOR PLOT
  CHL_TEST       = findgen(2000)*0.00025
  Table_data_SAT = fltarr([n_elements(CHL_TEST),500])
  Table_data_BEL = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
     minS   = 0.0
     maxS   = 1./Kd_PAR
     minB   = 1./Kd_PAR
     maxB   = 200
     Table_data_SAT(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
     Table_data_BEL(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor
  CHL_TEST       = findgen(2000)/500
  Table_data_SATN = fltarr([n_elements(CHL_TEST),500])
  Table_data_BELN = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
     minS   = 0.0
     maxS   = 1.
     minB   = 1.
     maxB   = 12
     Table_data_SATN(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
     Table_data_BELN(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor     
  
  ;;;PLOTS PROFILE 59
  multiplot,  /doxaxis, /doyaxis
  plot,  CHL_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,3.0], yrange=[10,0], symsize = 1.5, psym=5, $
     xminor = 4, /xstyle, /ystyle , /nodata, title = '2016-01-14', xtitle = '!8B!3!N / !8B!Ds!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,3.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,3.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  DCM_pop = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
  POP_ADD = MLD_pop+DCM_pop
  loadct, 39
  oplot, DCM_pop, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, POP_ADD, DEPTH_NOM, psym = 0, color = 255, thick = 7
  legend, ['(e)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
     charsize =4.0,CHARTHICK=4.0, thick = [4], position=[1.65,0]
  AXIS, 3.0,0.0,XAXIS = 1, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
     xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
     xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 3.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
     ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
     ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot,  BBP_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,3.0], yrange=[10,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N / !8b!Dbp,s!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,3.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,3.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop_bbp  =  MLD_pop*bbpS1
  DCM_pop_bbp  =  DCM_pop*bbpS2
  bbpk_pop_bbp = fltarr(n_elements(BBP_NOM))+bbpk
  POP_ADD_bbp = MLD_pop_bbp+DCM_pop_bbp+bbpk_pop_bbp
  loadct, 39
  oplot, DCM_pop_bbp, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp, DEPTH_NOM, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp, DEPTH_NOM, psym = 0, color = 255, thick = 7
  legend, ['(f)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[1.72,0]
  AXIS, 3.0,0.0,XAXIS = 1, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 3.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot, CHL,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.6], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8B!3!N !3(mg m!E-3!N!3)'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.6],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.6],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop*CHL_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop*CHL_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, POP_ADD*CHL_SAT, DEPTH, psym = 0, color = 255, thick = 7
  legend, ['(g)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.33,0]
  AXIS, 0.6,0.0,XAXIS = 1, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.6,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot, BBP,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.0015], yrange=[200,0], symsize = 1.5, psym=5, $
    xticks = 2, xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.0015],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.0015],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp*BBP_SAT, DEPTH, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp*BBP_SAT, DEPTH, psym = 0, color = 255, thick = 7
  legend, ['(h)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
     charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.00082,0]
  AXIS, 0.0015,0.0,XAXIS = 1, XRANGE=[0.0,0.0015], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.0015], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0015,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
     
  ;;;EXTRACT DATA FOR PROFILE 139 (NOTE IDL INDEX STARTS ON 0, SO INDEX 138)
  CHL = CHLA_ADJUSTED[*,138]
  AS = where(CHL gt 99990)
  CHL(AS) = !VALUES.F_NAN
  BBP = BBP700_ADJUSTED[*,138]
  AS = where(CHL gt 99990)
  BBP(AS) = !VALUES.F_NAN
  DEPTH  = DEPTH_ADJUSTED[*,138]
  AS = where(DEPTH gt 0 and DEPTH le 300)
  CHL   = CHL(AS)
  BBP   = BBP(AS)
  DEPTH = DEPTH(AS)

  ;;;MODEL PARAMATERS FROM FIT TO PROFILE 139 (Extracted from main Python code from the paper)
  Kd_PAR          = 0.051146009059312485
  Zp              = 4.6/Kd_PAR
  MLD             = 32.4
  CHL_SAT         = 0.01080000028014183
  BBP_SAT         = 0.0005762166692875326
  P1              = 8.181781454429675
  T1              = 3.3114394893159638
  Bm2             = 27.399422889286427
  T2              = 4.731820450679684
  Sig             = 1.5772477021328102
  bbpk            = 0.5070754821855764
  bbpS1           = 1 - bbpk
  bbpS2           = 0.018887802925762642
  CHL_NOM         = CHL / CHL_SAT
  BBP_NOM         = BBP / BBP_SAT
  DEPTH_NOM       = DEPTH * Kd_PAR
        
  ;;;SHADING FOR PLOT
  CHL_TEST       = findgen(2000)*0.00025
  Table_data_SAT = fltarr([n_elements(CHL_TEST),500])
  Table_data_BEL = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
     minS   = 0.0
     maxS   = 1./Kd_PAR
     minB   = 1./Kd_PAR
     maxB   = 200
     Table_data_SAT(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
     Table_data_BEL(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor
  CHL_TEST       = findgen(2000)/40
  Table_data_SATN = fltarr([n_elements(CHL_TEST),500])
  Table_data_BELN = fltarr([n_elements(CHL_TEST),500])
  for id = 0, n_elements(CHL_TEST)-1 do begin
     minS   = 0.0
     maxS   = 1.
     minB   = 1.
     maxB   = 12
     Table_data_SATN(id,*) = (findgen(500) * ((maxS - minS)/499.))+minS
     Table_data_BELN(id,*) = (findgen(500) * ((maxB - minB)/499.))+minB
  endfor

  ;;;PLOTS PROFILE 139
  multiplot,  /doxaxis, /doyaxis
  plot,  CHL_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,50.0], yrange=[10,0], symsize = 1.5, psym=5, $
     xminor = 4, /xstyle, /ystyle , /nodata, title = '2016-07-13', xtitle = '!8B!3!N / !8B!Ds!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,50.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,50.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  DCM_pop = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
  POP_ADD = MLD_pop+DCM_pop
  loadct, 39
  oplot, DCM_pop, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, POP_ADD, DEPTH_NOM, psym = 0, color = 255, thick = 7
  legend, ['(i)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
     charsize =4.0,CHARTHICK=4.0, thick = [4], position=[30,0]
  AXIS, 50.0,0.0,XAXIS = 1, XRANGE=[0.0,50.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,50.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 50.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot,  BBP_NOM, DEPTH_NOM, background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,3.0], yrange=[10,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N / !8b!Dbp,s!3!N'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SATN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SATN(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BELN(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BELN(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP_NOM, DEPTH_NOM, psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,3.0],[4.6,4.6], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,3.0],[MLD*KD_PAR,MLD*KD_PAR], psym = 0, thick = 8, linestyle = 1, color = 150
  MLD_pop_bbp  =  MLD_pop*bbpS1
  DCM_pop_bbp  =  DCM_pop*bbpS2
  bbpk_pop_bbp = fltarr(n_elements(BBP_NOM))+bbpk
  POP_ADD_bbp = MLD_pop_bbp+DCM_pop_bbp+bbpk_pop_bbp
  loadct, 39
  oplot, DCM_pop_bbp, DEPTH_NOM, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp, DEPTH_NOM, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp, DEPTH_NOM, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp, DEPTH_NOM, psym = 0, color = 255, thick = 7
  legend, ['(j)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[1.75,0]
  AXIS, 3.0,0.0,XAXIS = 1, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 10,XAXIS = 0, XRANGE=[0.0,3.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 3.0,0.0,YAXIS = 1, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 10,YAXIS = 0, YRANGE=[10,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot, CHL,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.4], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8B!3!N !3(mg m!E-3!N!3)'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, CHL,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.4],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.4],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop*CHL_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop*CHL_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, POP_ADD*CHL_SAT, DEPTH, psym = 0, color = 255, thick = 7
  legend, ['(k)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.22,0]
  TVLCT, [[0],[100],[0]], 150
  AXIS, 0.4,0.0,XAXIS = 1, XRANGE=[0.0,0.46], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.4], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.4,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doxaxis, /doyaxis
  plot, BBP,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.001], yrange=[200,0], symsize = 1.5, psym=5, $
    xticks = 2, xminor = 4, /xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)'
  TVLCT, [[176],[226],[255]], 9
  for id = 0l, n_elements(Table_data_SAT(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_SAT(*,id), psym = 0, color = 9, thick = 2
  endfor
  TVLCT, [[0],[178],[238]], 10
  for id = 0l, n_elements(Table_data_BEL(0,*))-1 do begin
    oplot, CHL_TEST, Table_data_BEL(*,id), psym = 0, color = 10, thick = 2
  endfor
  loadct, 0
  oplot, BBP,DEPTH,  psym = 1, color = 70, symsize = 2.0, thick = 3
  loadct, 39
  TVLCT, [[0],[100],[0]], 150
  oplot, [0.0,0.001],[4.6/Kd_PAR,4.6/Kd_PAR], psym = 0, thick = 8, linestyle = 5, color = 200
  oplot, [0.0,0.001],[MLD,MLD], psym = 0, thick = 8, linestyle = 1, color = 150
  oplot, DCM_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15,color = 50
  oplot, MLD_pop_bbp*BBP_SAT, DEPTH, psym = 0, thick = 15, color = 250
  oplot, bbpk_pop_bbp*BBP_SAT, DEPTH, psym = 0, color = 20, thick = 15
  oplot, POP_ADD_bbp*BBP_SAT, DEPTH, psym = 0, color = 255, thick = 7
  legend, ['(l)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.0006,0]
  AXIS, 0.001,0.0,XAXIS = 1, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.001,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,   Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save

  DEVICE, /CLOSE
  multiplot,/reset
  !X.MARGIN = 0
  !Y.MARGIN = 0
  !P.multi = 0
  Cleanplot

end
