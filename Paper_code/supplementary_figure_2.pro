; NAME:
;   Supplementary_Figure_2.pro
;
; PURPOSE:
;   Code used to produce Supplementary_Figure_2 from the Brewin et al. (In review) paper submitted to JGR-Oceans in 2022
;
; CATEGORY:
;   Data plotting
;
; CALLING SEQUENCE:
;   Supplementary_Figure_2
;
; INPUTS:
;   Model parameters and uncertainties from fits conducted in Python (see main Python code for the paper)
;
; OUTPUTS:
;   Supplementary Figure 2 to the paper
;
; MODIFICATION HISTORY:
; Bob Brewin 17/02/2022

pro Supplementary_Figure_2

  ;;;Float number
  FLOAT_NUMBER = '6901573'
  ;;;Define path to float data (NB: Needs to be edited depending on where data is)
  PATH_TO_FLOAT = '/Users/bobbrewin/Jupyter/Ideas/Two_comp_vertical_Chl/'
  ;;;Read the BGC-Argo netCDF file into IDL
  FILE_SPE = PATH_TO_FLOAT+FLOAT_NUMBER+'_Sprof.nc'
  
  ;;;SET PLOT UP
  SET_PLOT, 'PS'
  DEVICE, /color, /ENCAPSULATED, FILENAME = PATH_TO_FLOAT+'Figure_SENS_1.eps', XSIZE = 40, YSIZE = 25
  !p.font=0
  !X.MARGIN = [10, 3]
  !Y.MARGIN = [5, 2]
  device, /helvetica ;
  loadct, 39
  PSIZE=1.0
  PLOTSYM, 0, PSIZE, /FILL, THICK=1
  tau     = GREEK('tau')
  multiplot,[0,3,2,0,1],  ygap = 0.06, xgap = 0.035, /doyaxis, /doxaxis
  
  ;;;READ IN DATA FROM FLOAT
  filename_ig  = ncdf_open(FILE_SPE)
  ;;CORE ARGO
  Rrsid        = ncdf_varid(filename_ig,'PRES')
  ncdf_varget,filename_ig,Rrsid,PRES
  ;;;BASIC CONVERSION OF PRESSURE TO DEPTH
  DEPTH_ADJUSTED = ((((-1.82e-15  * PRES + 2.279e-10 ) * PRES - 2.2512e-5 ) * PRES + 9.72659) * PRES) / 9.81
  
  ;;;EXTRACT DEPTH FOR PROFILE 167 (NOTE IDL INDEX STARTS ON 0, SO INDEX 166)
  DEPTH  = DEPTH_ADJUSTED[*,166]
  AS = where(DEPTH gt 0 and DEPTH le 300)
  DEPTH = DEPTH(AS)

  ;;;MODEL PARAMATERS FROM FIT TO PROFILE 167 (Extracted from main Python code from the paper)
  Kd_PAR          = 0.048281593118789795
  Zp              = 4.6/Kd_PAR
  MLD             = 139.72
  CHL_SAT         = 0.17640002268665242
  BBP_SAT         = 0.000514139246661216
  P1              = 54.303707753866966
  T1              = 6.691624793773988
  Bm2             = 0.
  T2              = 0.
  Sig             = 0.
  bbpk            = 0.5290324727671847
  bbpS1           = 1. - bbpk
  bbpS2           = 0.

  ;;;MODEL RESULTS (BULK LINES IN PLOT) FOR PROFILE 167
  DEPTH_NOM          = DEPTH * Kd_PAR
  MLD_pop_OR         = 1. - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  MLD_pop_B_OR       = MLD_pop_OR*CHL_SAT
  MLD_pop_bbp_BBP_OR = MLD_pop_OR*bbpS1*BBP_SAT
  MLD_pop_bbp_BBK_OR = fltarr(n_elements(DEPTH))+bbpk*BBP_SAT

  ;;;UNCERTAINTIY IN MODEL PARAMATERS FOR PROFILE 167 (Some xtracted from main Python code from the paper)
  Kd_PAR_SD       = 0.048281593118789795 * 0.1 ;10% error on Kd
  MLD_SD          = 139.72 * 0.1               ;10% error on MLD
  CHL_SAT_SD      = CHL_SAT * 0.1              ;10% error on Surf CHl 
  BBP_SAT_SD      = 0.0001657995605532726      ; SD Surf bbp
 
  ;;;UPPER AND LOWER BOUNDS TO UNCERTAINTIY IN MODEL PARAMATERS FOR PROFILE 167
  Kd_PAR_ARR  = [Kd_PAR-Kd_PAR_SD,Kd_PAR+Kd_PAR_SD]
  MLD_ARR     = [MLD-MLD_SD,MLD+MLD_SD]
  CHL_SAT_ARR = [CHL_SAT-CHL_SAT_SD,CHL_SAT+CHL_SAT_SD]
  BBP_SAT_ARR = [BBP_SAT-BBP_SAT_SD,BBP_SAT+BBP_SAT_SD]
  P1_ARR      = [46.781311610770466,65.82597068307717]   ;LB and UB Bootstap
  T1_ARR      = [6.667055718233579,6.725841894532686]    ;LB and UB Bootstap
  bbpk_ARR    = [0.5043385312255136,0.5609642779672706]  ;LB and UB Bootstap

  ;;ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 167
  MLD_pop_B     = fltarr([128,n_elements(DEPTH)])
  MLD_pop_bbp   = fltarr([128,n_elements(DEPTH)])
  bbpk_pop_bbp  = fltarr([128,n_elements(DEPTH)])
  cnt = 0
  for id = 0, n_elements(Kd_PAR_ARR)-1 do begin
    Kd_PAR = Kd_PAR_ARR[id]
    for ie = 0, n_elements(MLD_ARR)-1 do begin
      MLD = MLD_ARR[ie]
      for ig = 0, n_elements(CHL_SAT_ARR)-1 do begin
        CHL_SAT = CHL_SAT_ARR[ig]
        for ih = 0, n_elements(BBP_SAT_ARR)-1 do begin
          BBP_SAT = BBP_SAT_ARR[ih]
           for ii = 0, n_elements(P1_ARR)-1 do begin
            P1 = P1_ARR[ii]
            for ij = 0, n_elements(T1_ARR)-1 do begin
              T1 = T1_ARR[ij]
              for ik = 0, n_elements(bbpk_ARR)-1 do begin
                bbpk = bbpk_ARR[ik] 
                ;Start the processing
                bbpS1                = 1 - bbpk
                DEPTH_NOM            = DEPTH * Kd_PAR
                MLD_pop              =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
                MLD_pop_B[cnt,*]     =  MLD_pop*CHL_SAT
                MLD_pop_bbp[cnt,*]   =  MLD_pop*bbpS1 * BBP_SAT
                bbpk_pop_bbp[cnt,*]  = fltarr(n_elements(DEPTH_NOM))+bbpk*BBP_SAT
                cnt = cnt + 1
              endfor
            endfor
           endfor
        endfor
      endfor
    endfor
  endfor
  
  ;;;SHADING FOR ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 167 (MIN AND MAX OF ENSEMBLE)
  Table_data_MLD_B    = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbp  = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbpk = fltarr([100,n_elements(DEPTH)])
  for id = 0, n_elements(DEPTH)-1 do begin
    minS   = min(MLD_pop_B[*,id])
    maxS   = max(MLD_pop_B[*,id])
    Table_data_MLD_B(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(MLD_pop_bbp[*,id])
    maxS   = max(MLD_pop_bbp[*,id])
    Table_data_MLD_bbp(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(bbpk_pop_bbp[*,id])
    maxS   = max(bbpk_pop_bbp[*,id])
    Table_data_MLD_bbpk(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
  endfor
  
  ;;;PLOTS PROFILE 167
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.4], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4,yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata, title = '2017-01-04', xtitle = '!8B!3!N !3(mg m!E-3!N!3)', ytitle = 'Depth (m)';, $
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_B(*,0))-1 do begin
    oplot, Table_data_MLD_B(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  loadct, 39
  oplot, MLD_pop_B_OR, DEPTH, psym = 0, thick = 15, color = 250
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.2,160]
  legend, ['(a)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.22,0]
  TVLCT, [[0],[100],[0]], 150
  AXIS, 0.4,0.0,XAXIS = 1, XRANGE=[0.0,0.4], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.4], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.4,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  yminor = 4,yticks = 4, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save 
  multiplot, /doyaxis, /doxaxis
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.001], yrange=[200,0], symsize = 1.5, psym=5, $
    xticks = 2, xminor = 4, yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata,  xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)', ytitle = 'Depth (m)';;;;, YTickFormat='exponent';, $
  TVLCT, [[147],[112],[219]], 150
  for id = 0l, n_elements(Table_data_MLD_bbpk(*,0))-1 do begin
    oplot, Table_data_MLD_bbpk(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBK_OR, DEPTH, psym = 0, color = 20, thick = 15
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_bbp(*,0))-1 do begin
    oplot, Table_data_MLD_bbp(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBP_OR, DEPTH, psym = 0, thick = 15, color = 250
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,145]
  LEGEND, ['!N!8b!S!E k!R!Ibp!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[20], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,160]
  legend, ['(b)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.00055,0]
  AXIS, 0.001,0.0,XAXIS = 1, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.001,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,  Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save

  ;;;EXTRACT DEPTH FOR PROFILE 59 (NOTE IDL INDEX STARTS ON 0, SO INDEX 58)
  DEPTH  = DEPTH_ADJUSTED[*,58]
  AS = where(DEPTH gt 0 and DEPTH le 300)
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

  ;;;MODEL RESULTS (BULK LINES IN PLOT) FOR PROFILE 59
  DEPTH_NOM          = DEPTH * Kd_PAR
  MLD_pop_OR         =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  MLD_pop_B_OR       =  MLD_pop_OR*CHL_SAT
  DCM_pop_OR         = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
  DCM_pop_B_OR       = DCM_pop_OR*CHL_SAT
  MLD_pop_bbp_BBP_OR = MLD_pop_OR*bbpS1*BBP_SAT
  DCM_pop_bbp_BBP_OR = DCM_pop_OR*bbpS2*BBP_SAT
  MLD_pop_bbp_BBK_OR =  fltarr(n_elements(DEPTH))+bbpk*BBP_SAT
  
  ;;;UNCERTAINTIY IN MODEL PARAMATERS FOR PROFILE 59 (Some xtracted from main Python code from the paper)
  Kd_PAR_SD       = Kd_PAR * 0.1          ;10% error on Kd
  MLD_SD          = MLD * 0.1             ;10% error on MLD
  CHL_SAT_SD      = CHL_SAT * 0.1         ;10% error on Surf Chl
  BBP_SAT_SD      = 0.00011046154325102882; SD Surf bbp
  
  ;;;UPPER AND LOWER BOUNDS TO UNCERTAINTIY IN MODEL PARAMATERS FOR PROFILE 59
  Kd_PAR_ARR      = [Kd_PAR-Kd_PAR_SD,Kd_PAR+Kd_PAR_SD]
  MLD_ARR         = [MLD-MLD_SD,MLD+MLD_SD]
  CHL_SAT_ARR     = [CHL_SAT-CHL_SAT_SD,CHL_SAT+CHL_SAT_SD]
  BBP_SAT_ARR     = [BBP_SAT-BBP_SAT_SD,BBP_SAT+BBP_SAT_SD]
  Bm2_ARR         = [1.513945586079651,1.7592450891025422]   ;LB and UB Bootstap
  T2_ARR          = [5.187569587873273,5.326843331650464]    ;LB and UB Bootstap
  Sig_ARR         = [1.4051091567498055,1.741277834561453]   ;LB and UB Bootstap
  T1_COFF1_ARR    = [0.53,0.71]                              ;confidence intervals of regression
  T1_COFF2_ARR    = [1.52000,3.06000]                        ;confidence intervals of regression
  P1_COFF1_ARR    = [0.07,0.09]                              ;confidence intervals of regression
  P1_COFF2_ARR    = [0.65,0.66]                              ;confidence intervals of regression
  bbpk_ARR        = [0.3786210724209813,0.40551147105812624] ;LB and UB Bootstap
  bbpS2_ARR       = [0.4134414327530018,0.5303675665572615]  ;LB and UB Bootstap
  
  ;;;ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 59
  MLD_pop_B     = fltarr([8192,n_elements(DEPTH)])
  DCM_pop_B     = fltarr([8192,n_elements(DEPTH)])
  MLD_pop_bbp   = fltarr([8192,n_elements(DEPTH)])
  DCM_pop_bbp   = fltarr([8192,n_elements(DEPTH)])
  bbpk_pop_bbp  = fltarr([8192,n_elements(DEPTH)])
  cnt = 0
  for id = 0, n_elements(Kd_PAR_ARR)-1 do begin
    Kd_PAR = Kd_PAR_ARR[id]
    for ie = 0, n_elements(MLD_ARR)-1 do begin
      MLD = MLD_ARR[ie]
      for ig = 0, n_elements(CHL_SAT_ARR)-1 do begin
        CHL_SAT = CHL_SAT_ARR[ig]
        for ih = 0, n_elements(BBP_SAT_ARR)-1 do begin
          BBP_SAT = BBP_SAT_ARR[ih]
          for ii = 0, n_elements(P1_COFF1_ARR)-1 do begin
            P1_COFF1 = P1_COFF1_ARR[ii]
            for ij = 0, n_elements(P1_COFF2_ARR)-1 do begin
              P1_COFF2 = P1_COFF2_ARR[ij]              
              for ik = 0, n_elements(T1_COFF1_ARR)-1 do begin
                T1_COFF1 = T1_COFF1_ARR[ik]
                for il = 0, n_elements(T1_COFF2_ARR)-1 do begin
                  T1_COFF2 = T1_COFF2_ARR[il]
                  for im = 0, n_elements(Bm2_ARR)-1 do begin
                    Bm2 = Bm2_ARR[im]
                    for io = 0, n_elements(T2_ARR)-1 do begin
                      T2 = T2_ARR[io]
                      for ip = 0, n_elements(Sig_ARR)-1 do begin
                        Sig = Sig_ARR[ip]
                        for iq = 0, n_elements(bbpk_ARR)-1 do begin
                          bbpk = bbpk_ARR[iq]
                          for ir = 0, n_elements(bbpS2_ARR)-1 do begin
                            bbpS2 = bbpS2_ARR[ir]
                            ;Start the processing
                            bbpS1                = 1 - bbpk
                            DEPTH_NOM            = DEPTH * Kd_PAR
                            T1 = T1_COFF1*MLD*Kd_PAR+T1_COFF2
                            P1 = 10.^(P1_COFF1*T1+P1_COFF2)
                            MLD_pop              =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
                            MLD_pop_B[cnt,*]     =  MLD_pop*CHL_SAT
                            MLD_pop_bbp[cnt,*]   =  MLD_pop*bbpS1 * BBP_SAT
                            bbpk_pop_bbp[cnt,*]  = fltarr(n_elements(DEPTH_NOM))+bbpk*BBP_SAT
                            DCM_pop              = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
                            DCM_pop_B[cnt,*]     = DCM_pop*CHL_SAT
                            DCM_pop_bbp[cnt,*]   = DCM_pop*bbpS2 * BBP_SAT
                            cnt = cnt + 1
                          endfor
                        endfor
                      endfor
                    endfor
                  endfor
                endfor
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor

  ;;;SHADING FOR ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 59 (MIN AND MAX OF ENSEMBLE)
  Table_data_MLD_B    = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbp  = fltarr([100,n_elements(DEPTH)])
  Table_data_DCM_B    = fltarr([100,n_elements(DEPTH)])
  Table_data_DCM_bbp  = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbpk = fltarr([100,n_elements(DEPTH)])
  for id = 0, n_elements(DEPTH)-1 do begin
    minS   = min(MLD_pop_B[*,id])
    maxS   = max(MLD_pop_B[*,id])
    Table_data_MLD_B(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(MLD_pop_bbp[*,id])
    maxS   = max(MLD_pop_bbp[*,id])
    Table_data_MLD_bbp(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(DCM_pop_B[*,id])
    maxS   = max(DCM_pop_B[*,id])
    Table_data_DCM_B(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(DCM_pop_bbp[*,id])
    maxS   = max(DCM_pop_bbp[*,id])
    Table_data_DCM_bbp(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS    
    minS   = min(bbpk_pop_bbp[*,id])
    maxS   = max(bbpk_pop_bbp[*,id])
    Table_data_MLD_bbpk(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
  endfor

  ;;;PLOTS PROFILE 59
  multiplot, /doyaxis, /doxaxis
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.6], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4,yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata, title = '2016-01-14',  xtitle = '!8B!3!N !3(mg m!E-3!N!3)'
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_B(*,0))-1 do begin
    oplot, Table_data_MLD_B(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  loadct, 39
  oplot, MLD_pop_B_OR, DEPTH, psym = 0, thick = 15, color = 250
  TVLCT, [[135],[206],[255]], 150
  for id = 0l, n_elements(Table_data_DCM_B(*,0))-1 do begin
    oplot, Table_data_DCM_B(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  loadct, 39
  oplot, DCM_pop_B_OR, DEPTH, psym = 0, thick = 15, color = 50
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.3,160]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
      usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
      thick =[15], position=[0.3,175]
  legend, ['(c)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.33,0]
  TVLCT, [[0],[100],[0]], 150
  AXIS, 0.6,0.0,XAXIS = 1, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.6,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  yminor = 4,yticks = 4, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doyaxis, /doxaxis
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.001], yrange=[200,0], symsize = 1.5, psym=5, $
    xticks = 2, xminor = 4, yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)'
  TVLCT, [[147],[112],[219]], 150
  for id = 0l, n_elements(Table_data_MLD_bbpk(*,0))-1 do begin
    oplot, Table_data_MLD_bbpk(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBK_OR, DEPTH, psym = 0, color = 20, thick = 15
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_bbp(*,0))-1 do begin
    oplot, Table_data_MLD_bbp(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBP_OR, DEPTH, psym = 0, thick = 15, color = 250
  TVLCT, [[135],[206],[255]], 150
  for id = 0l, n_elements(Table_data_DCM_bbp(*,0))-1 do begin
    oplot, Table_data_DCM_bbp(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, DCM_pop_bbp_BBP_OR, DEPTH, psym = 0, color = 50, thick = 15
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,145]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,160]
  LEGEND, ['!N!8b!S!E k!R!Ibp!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[20], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,175]
  legend, ['(d)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.00055,0]
  AXIS, 0.001,0.0,XAXIS = 1, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.001,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,  Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
    
  ;;;EXTRACT DEPTH FOR PROFILE 139 (NOTE IDL INDEX STARTS ON 0, SO INDEX 138)
  DEPTH  = DEPTH_ADJUSTED[*,138]
  AS = where(DEPTH gt 0 and DEPTH le 300)
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

  ;;;MODEL RESULTS (BULK LINES IN PLOT) FOR PROFILE 139
  DEPTH_NOM       = DEPTH * Kd_PAR
  MLD_pop_OR      =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
  MLD_pop_B_OR    =  MLD_pop_OR*CHL_SAT
  DCM_pop_OR      = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
  DCM_pop_B_OR    = DCM_pop_OR*CHL_SAT
  MLD_pop_bbp_BBP_OR = MLD_pop_OR*bbpS1*BBP_SAT
  DCM_pop_bbp_BBP_OR = DCM_pop_OR*bbpS2*BBP_SAT
  MLD_pop_bbp_BBK_OR =  fltarr(n_elements(DEPTH))+bbpk*BBP_SAT
    
  ;;;UNCERTAINTIY IN MODEL PARAMATERS FOR PROFILE 139 (Some xtracted from main Python code from the paper)
  Kd_PAR_SD       = Kd_PAR * 0.1    ;10% error on Kd
  MLD_SD          = MLD * 0.1       ;10% error on MLD
  CHL_SAT_SD      = CHL_SAT * 0.1   ;10% error on Surf Chl-a
  BBP_SAT_SD      = 3.6293092546053326e-05 ; SD Surf bbp
  
  ;;;UPPER AND LOWER BOUNDS TO UNCERTAINTIY IN MODELPARAMATERS FOR PROFILE 59
  Kd_PAR_ARR      = [Kd_PAR-Kd_PAR_SD,Kd_PAR+Kd_PAR_SD]
  MLD_ARR         = [MLD-MLD_SD,MLD+MLD_SD]
  CHL_SAT_ARR     = [CHL_SAT-CHL_SAT_SD,CHL_SAT+CHL_SAT_SD]
  BBP_SAT_ARR     = [BBP_SAT-BBP_SAT_SD,BBP_SAT+BBP_SAT_SD]
  Bm2_ARR         = [26.906756120268298,27.926025006580804] ;LB and UB Bootstap
  T2_ARR          = [4.697027941750102,4.774006693248174]   ;LB and UB Bootstap
  Sig_ARR         = [1.5656651053111381,1.591257476486047]  ;LB and UB Bootstap
  T1_COFF1_ARR    = [0.53,0.71]                             ;confidence intervals of regression
  T1_COFF2_ARR    = [1.52000,3.06000]                       ;confidence intervals of regression
  P1_COFF1_ARR    = [0.07,0.09]                             ;confidence intervals of regression
  P1_COFF2_ARR    = [0.65,0.66]                             ;confidence intervals of regression
  bbpk_ARR        = [0.4909940920094005,0.527436495923125]  ;LB and UB Bootstap
  bbpS2_ARR       = [0.017283298106479106,0.020667469567314904];LB and UB Bootstap

  ;;;ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 139
  MLD_pop_B     = fltarr([8192,n_elements(DEPTH)])
  DCM_pop_B     = fltarr([8192,n_elements(DEPTH)])
  MLD_pop_bbp   = fltarr([8192,n_elements(DEPTH)])
  DCM_pop_bbp   = fltarr([8192,n_elements(DEPTH)])
  bbpk_pop_bbp  = fltarr([8192,n_elements(DEPTH)])
  cnt = 0
  for id = 0, n_elements(Kd_PAR_ARR)-1 do begin
    Kd_PAR = Kd_PAR_ARR[id]
    for ie = 0, n_elements(MLD_ARR)-1 do begin
      MLD = MLD_ARR[ie]
      for ig = 0, n_elements(CHL_SAT_ARR)-1 do begin
        CHL_SAT = CHL_SAT_ARR[ig]
        for ih = 0, n_elements(BBP_SAT_ARR)-1 do begin
          BBP_SAT = BBP_SAT_ARR[ih]
          for ii = 0, n_elements(P1_COFF1_ARR)-1 do begin
            P1_COFF1 = P1_COFF1_ARR[ii]
            for ij = 0, n_elements(P1_COFF2_ARR)-1 do begin
              P1_COFF2 = P1_COFF2_ARR[ij]
              for ik = 0, n_elements(T1_COFF1_ARR)-1 do begin
                T1_COFF1 = T1_COFF1_ARR[ik]
                for il = 0, n_elements(T1_COFF2_ARR)-1 do begin
                  T1_COFF2 = T1_COFF2_ARR[il]
                  for im = 0, n_elements(Bm2_ARR)-1 do begin
                    Bm2 = Bm2_ARR[im]
                    for io = 0, n_elements(T2_ARR)-1 do begin
                      T2 = T2_ARR[io]
                      for ip = 0, n_elements(Sig_ARR)-1 do begin
                        Sig = Sig_ARR[ip]
                        for iq = 0, n_elements(bbpk_ARR)-1 do begin
                          bbpk = bbpk_ARR[iq]
                          for ir = 0, n_elements(bbpS2_ARR)-1 do begin
                            bbpS2 = bbpS2_ARR[ir]
                            ;Start the processing
                            bbpS1                = 1 - bbpk
                            DEPTH_NOM            = DEPTH * Kd_PAR
                            T1 = T1_COFF1*MLD*Kd_PAR+T1_COFF2
                            P1 = 10.^(P1_COFF1*T1+P1_COFF2)
                            MLD_pop              =  1 - (1./(1+exp(-(P1/T1)*(DEPTH_NOM-T1))))
                            MLD_pop_B[cnt,*]     =  MLD_pop*CHL_SAT
                            MLD_pop_bbp[cnt,*]   =  MLD_pop*bbpS1 * BBP_SAT
                            bbpk_pop_bbp[cnt,*]  = fltarr(n_elements(DEPTH_NOM))+bbpk*BBP_SAT
                            DCM_pop              = (Bm2*exp(-((DEPTH_NOM - T2)/(Sig))^2.))
                            DCM_pop_B[cnt,*]     = DCM_pop*CHL_SAT
                            DCM_pop_bbp[cnt,*]   = DCM_pop*bbpS2 * BBP_SAT
                            cnt = cnt + 1
                          endfor
                        endfor
                      endfor
                    endfor
                  endfor
                endfor
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor

  ;;;SHADING FOR ENSEMBLE SIMULATION OF UNCERTAINTIY IN MODEL FOR PROFILE 139 (MIN AND MAX OF ENSEMBLE)
  Table_data_MLD_B    = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbp  = fltarr([100,n_elements(DEPTH)])
  Table_data_DCM_B    = fltarr([100,n_elements(DEPTH)])
  Table_data_DCM_bbp  = fltarr([100,n_elements(DEPTH)])
  Table_data_MLD_bbpk = fltarr([100,n_elements(DEPTH)])
  for id = 0, n_elements(DEPTH)-1 do begin
    minS   = min(MLD_pop_B[*,id])
    maxS   = max(MLD_pop_B[*,id])
    Table_data_MLD_B(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(MLD_pop_bbp[*,id])
    maxS   = max(MLD_pop_bbp[*,id])
    Table_data_MLD_bbp(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(DCM_pop_B[*,id])
    maxS   = max(DCM_pop_B[*,id])
    Table_data_DCM_B(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(DCM_pop_bbp[*,id])
    maxS   = max(DCM_pop_bbp[*,id])
    Table_data_DCM_bbp(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
    minS   = min(bbpk_pop_bbp[*,id])
    maxS   = max(bbpk_pop_bbp[*,id])
    Table_data_MLD_bbpk(*,id) = (findgen(100) * ((maxS - minS)/99.))+minS
  endfor  
  
  ;;;PLOTS PROFILE 139
  multiplot, /doyaxis, /doxaxis
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.6], yrange=[200,0], symsize = 1.5, psym=5, $
    xminor = 4,yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata, title = '2016-07-13', xtitle = '!8B!3!N !3(mg m!E-3!N!3)'
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_B(*,0))-1 do begin
    oplot, Table_data_MLD_B(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  loadct, 39
  oplot, MLD_pop_B_OR, DEPTH, psym = 0, thick = 15, color = 250
  TVLCT, [[135],[206],[255]], 150
  for id = 0l, n_elements(Table_data_DCM_B(*,0))-1 do begin
    oplot, Table_data_DCM_B(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  loadct, 39
  oplot, DCM_pop_B_OR, DEPTH, psym = 0, thick = 15, color = 50
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.3,160]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.3,175]
  legend, ['(e)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.33,0]
  TVLCT, [[0],[100],[0]], 150
  AXIS, 0.6,0.0,XAXIS = 1, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.6], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.6,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3,  yminor = 4,yticks = 4, Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  multiplot, /doyaxis, /doxaxis
  plot, MLD_pop_B,DEPTH,  background = 255, color = 0, charsize =2.,charthick = 2.5,  XTHICK =4, YTHICK =4, xrange = [0.0,0.001], yrange=[200,0], symsize = 1.5, psym=5, $
    xticks = 2, xminor = 4, yminor = 4,yticks = 4,/xstyle, /ystyle , /nodata, xtitle = '!8b!Dbp!3!N !3(m!E-1!N!3)'
  TVLCT, [[147],[112],[219]], 150
  for id = 0l, n_elements(Table_data_MLD_bbpk(*,0))-1 do begin
    oplot, Table_data_MLD_bbpk(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBK_OR, DEPTH, psym = 0, color = 20, thick = 15
  TVLCT, [[255],[181],[197]], 150
  for id = 0l, n_elements(Table_data_MLD_bbp(*,0))-1 do begin
    oplot, Table_data_MLD_bbp(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, MLD_pop_bbp_BBP_OR, DEPTH, psym = 0, thick = 15, color = 250
  TVLCT, [[135],[206],[255]], 150
  for id = 0l, n_elements(Table_data_DCM_bbp(*,0))-1 do begin
    oplot, Table_data_DCM_bbp(id,*), DEPTH, psym = 0, color = 150, thick = 5
  endfor
  oplot, DCM_pop_bbp_BBP_OR, DEPTH, psym = 0, color = 50, thick = 15
  LEGEND, ['1'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[250], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,145]
  LEGEND, ['2'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[50], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,160]
  LEGEND, ['!N!8b!S!E k!R!Ibp!N!3'], /fill,textcolors = 0, psym=[0], linestyle = [0], $
    usersym=usersym,colors=[20], BOX=0, spacing =1.0, charsize =1.5,CHARTHICK=2, $
    thick =[15], position=[0.0005,175]
  legend, ['(f)'],/fill,psym=[1],linestyle = [-1], box = 0, textcolor = 0,$
    charsize =4.0,CHARTHICK=4.0, thick = [4], position=[0.00055,0]
  AXIS, 0.001,0.0,XAXIS = 1, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.0, 200,XAXIS = 0, XRANGE=[0.0,0.001], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, xticks = 2, xminor = 4, xthick = 4,  $
    xtickname = replicate(' ',8), /XSTYLE , /save
  AXIS, 0.001,0.0,YAXIS = 1, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,  Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save
  AXIS, 0.0, 200,YAXIS = 0, YRANGE=[150,0.0], color = 0, FONT=-1, CHARSIZE=3.,CHARTHICK=3, yminor = 4,yticks = 4,Ythick = 4,  $
    ytickname = replicate(' ',9), /YSTYLE , /save   
    
  DEVICE, /CLOSE
  multiplot,/reset
  !X.MARGIN = 0
  !Y.MARGIN = 0
  !P.multi = 0
  Cleanplot

end
  

