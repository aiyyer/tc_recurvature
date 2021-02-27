;***************************************************************************
;**********                                                    *************
;**********  A "quick and dirty" program to determine the      *************
;**********  recurvature location of tropical cyclones         *************
;**********  1) reads ibtracs data for user specified dates    *************
;**********  2) output netcdf file with recurvature info       *************
;**********  3) output optional images of each tc track        *************
;**********                                                    *************
;***************************************************************************
;***************************************************************************
; A. Aiyyer 
;
; Citations for this code:
;
; Aiyyer(2015): Recurving western North Pacific tropical cyclones
;               and midlatitude predictability.  GRL, AGU.
;               https://doi.org/10.1002/2015GL065082
;
; Aiyyer and Wade (2021): Acceleration of Tropical Cyclones As a Proxy For
;               Extratropical Interactions: Synoptic-Scale Patterns
;               and Long-Term Trends, Weather Clim. Dynam. Discuss.
;               https://doi.org/10.5194/wcd-2021-4,
;
;----------------------------------------------------------------------------
;
;
; Notes:
; 1. All recurvatures are defined to be constrained by westward movement
;    prior to, and eastward movement after
;
;
; 2. Some storms have multiple recurvature points owing to kinks in the track
;    or "reasonable" recurvatures. I recommend visual inspection and weeding out
;    spurious cases
;;***************************************************************************
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


UNDER DEVELOPMENT....

Need to tweak the recurvature detection to ensure minor kinks are disregarded



;


load "rad_curv.ncl"

begin
;
;
;===========================================================================  
; User inputs
             
; define the ibtracs file and its location

  ibtracsFile = "IBTrACS.since1980.v04r00.nc"
  dir         = "./"
  startYear   = 2019
  endYear     = 2019

; Want to create an image with all recurving tracks (each storm in separate frame)  
  iPlot = 1 ; 1 = yes or 0 = no
  
  basinToConsider = "WP"

  LonSign = 0 ; set to 1 if you need to convert Negative Longitudes (i.e., -180 to 0 ) to
              ; positive (i.e. 180 - 360)
;  
; basin choices: EP=East_Pacific NA=North_Atlantic NI=North_Indian
;                SA=South_Atlantic SI=South_Indian SP=South_Pacific
;                 WP=Western_Pacific
;===========================================================================
; END USER INPUTS  
;===========================================================================
;
;
  
  fileName = dir + ibtracsFile
  fileIn = addfile(fileName, "r" )


; read time stamp for the first record for each TC from ibtracs
  time        = fileIn->time(:,0)
  basin       = chartostring(fileIn->basin)

  dateBeg = cd_inv_calendar (startYear,01,01,00,00,00,time@units,0)
  dateEnd = cd_inv_calendar (endYear,12,31,23,00,00,time@units,0)
  
  indices   :=  ind(basin(:,0) .eq. basinToConsider .and. time .ge. dateBeg .and. time .le. dateEnd) 
  delete(time)
  
;-----------------------------------------------------------------------  
  if ( iPlot .eq. 1 ) then  
    wks  = gsn_open_wks("ncgm","recurv")      ; open as postscript file
    res = True
    res@mpOutlineOn         = True
    res@mpMaxLatF           = 50. ;60.0
    res@mpMinLatF           = 0. ;0.0
    res@mpMinLonF           = 110. ;100. ;200  ;100.
    res@mpMaxLonF           = 180. ;220. ;340. ;220
    res@mpCenterLonF = .5*( res@mpMinLonF +  res@mpMaxLonF)
    res@gsnMajorLatSpacing  = 10
    res@gsnMajorLonSpacing  = 20
    res@mpGeophysicalLineThicknessF = 1.   
    gsres               = True
    gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
    gsres@gsMarkerSizeF     = 20.          ; polymarker size
    res@mpFillOn            = False
    recur = True                          ; start of traj resources
    res@gsnDraw          = False                   
    res@gsnFrame         = False
    dumG  = new(20,"graphic")
    colors  = (/"orange","navy blue","red","dark green"/)
  end if
  
;-----------------------------------------------------------------------    
;
; Read all names,dates and basin type in the file

  name        = charactertostring (fileIn->name)
  Filedates   = fileIn->time
  Lat         = flt2dble(fileIn->lat)
  Lon         = flt2dble(fileIn->lon)
  numObs      = tointeger(short2flt (fileIn->numobs))
  Sid         = chartostring(fileIn->sid) ; stormid

  dimw  = dimsizes(indices)
  print ("No of storms in database for this basin within specified dates = " + dimw(0) )
  NSTORMS  = dimw(0) ; NSTORMS is the total number of storms in this basin that match the dates


  maxRecur = 6

; Prefix Rec for arrays of data for recurvature points  
  RecLat    =  new ( (/NSTORMS,maxRecur/), typeof(Lat))
  RecLon    =  new ( (/NSTORMS,maxRecur/), typeof(Lon))
  RecDate   =  new ( (/NSTORMS,maxRecur/), typeof(Filedates) )

; radius of recuvatures  
  RecRad    =  new ( (/NSTORMS,maxRecur/), typeof(Lat))


  
  RecNam    =  new ( (/NSTORMS/), typeof(name) )
  RecNum    =  new ( (/NSTORMS/), "integer" )
  RecSid    =  new ( (/NSTORMS/), typeof(Sid) )
  RecLat@_FillValue = Lat@_FillValue
  RecLon@_FillValue = Lon@_FillValue
  RecDate@_FillValue = Filedates@_FillValue
  RecLat  = RecLat@_FillValue
  RecLon  = RecLon@_FillValue
  RecDate = RecDate@_FillValue
  
  copy_VarAtts(Lat,RecLat)
  copy_VarAtts(Lon,RecLon)
  copy_VarAtts(Filedates,RecDate)
  copy_VarAtts(Sid,RecSid)
  copy_VarAtts(name,RecNam)
  copy_VarAtts(RecSid,Sid)

  
  NumStormRecurving = 0 ; no of recurving storms found by the code


;  
; Now loop over each storm and check if it recurves

  do i = 0,NSTORMS-1
    
    istorm     = indices(i)
    ntime     = numObs(istorm)    
    sName     = name(istorm)
    RecNam(i) = sName
    RecSid(i) = Sid(istorm)
    
    print ( "Now processing Storm = " + sName + " " + cd_calendar(Filedates(istorm,0),-3)/1000000 + " StormId = " + RecSid(i) )
   
    
;......................................................................    
; occassionaly track data has additional reports. For consistency
; first extract only 3 hourly reports at synoptic hours 0,6,12,18
    
    
    print ( "ntimes: all track data = " + ntime)
    
; create new temporary arrays for dates,lat,lon by copying old arrays
    LonSP  := new( ntime, typeof(Lon) )
    LatSP  := new( ntime, typeof(Lat) )
    TimeSP := new( ntime, typeof(Filedates) )
    TimeSP@units = Filedates@units
    LonSP = LonSP@_FillValue
    LatSP = LatSP@_FillValue
    TimeSP = TimeSP@_FillValue
    is = 0   
    do j = 0,ntime-1      
      curTime = cd_calendar(Filedates(istorm,j),-5)
      hour    = curTime(0,3)
; only consider the three hourly reports
      if ( mod(hour,3) .eq. 0 ) then       
        TimeSP(is)     = (/Filedates(istorm,j)/)
         LatSP(is)      = (/Lat(istorm,j)/)
         LonSP(is)      = (/Lon(istorm,j)/)
        is = is + 1
      end if      
    end do

; reset the ntime to reflect only the 3 hourly reports
; shrink the arrays (happens only if ntime changes)    
    ntime   = is
    LonSP  := LonSP(0:ntime-1)
    LatSP  := LatSP(0:ntime-1)
    TimeSP := TimeSP(0:ntime-1)
;......................................................................
 
    print ( "ntimes 3h track data = " + ntime)
    
; if needed convert negative lons to positive
    if (LonSign .eq. 1 ) then
      LonSP = where(LonSP.lt.0, LonSP+360, LonSP) 
    end if
    
; lets smooth the lat lon a little     
    LL :=  LonSP
    LA :=  LatSP
    
; simple 1-2-1 filter    
;    do ii =  2,ntime-2
;      LL(ii) = .25*(LonSP(ii-1) + 2.*LonSP(ii) + LonSP(ii+1))
;      LA(ii) = .25*(LatSP(ii-1) + 2.*LatSP(ii) + LatSP(ii+1))
;    end do

;Savitzky–Golay filter 1
    do ii =  3,ntime-3
      LL(ii) = (-3.*LonSP(ii-2) + 12.*LonSP(ii-1) + 17.*LonSP(ii) +  12.*LonSP(ii+1) - 3.*LonSP(ii+2))/35.0
      LA(ii) = (-3.*LatSP(ii-2) + 12.*LatSP(ii-1) + 17.*LatSP(ii) +  12.*LatSP(ii+1) - 3.*LatSP(ii+2))/35.0
    end do

;Savitzky–Golay filter 2
;    do ii =  4,ntime-4
;      LL(ii) = (5.*LonSP(ii-3) + -30.*LonSP(ii-2) + 75.*LonSP(ii-1) + 131.*LonSP(ii) +  75.*LonSP(ii+1) - 30.*LonSP(ii+2) + 5.*LonSP(ii+3)  )/231.
;      LA(ii) = (5.*LatSP(ii-3) + -30.*LatSP(ii-2) + 75.*LatSP(ii-1) + 131.*LatSP(ii) +  75.*LatSP(ii+1) - 30.*LatSP(ii+2) + 5.*LatSP(ii+3)  )/231.
;    end do


;

; calculate the radius of curvature of all points along the track
    
    RadCurv := new(ntime, typeof(LA) )
    RadCurv = rad_curv(LA,LL,ntime)
    ; convert to km
    RadCurv = RadCurv*1.e-3 
    RadCurv@units = "km"
    RadCurv@name =  "Radius of curvature"

    jFound  = 0
    irecurv = 0
    
    do it = 4,ntime-4 ; start a little further into the track  
      val1 =  LL(it+1)-LL(it)
      val2 =  LL(it)-LL(it-1)
      val3 =  LL(it+2) - 2.*LL(it) + LL(it-2)
      val4 =  LA(it+3)  - LA(it)
      val5 =  LonSP(it+3) - LonSP(it)

      if ( (val1 .ge. 0. ) .and. ( val2 .le. 0. ) .and. (val3 .ge. 0. ) ) then              
        if (  (val5   .ge. 0.) .and.(val4   .ge. 0.) ) then                                      
          jFound = 1                   
          RecLat (i,irecurv)   =   LatSP(it)
          RecLon (i,irecurv)   =   LonSP(it)           
          RecDate (i,irecurv)  =  TimeSP(it)
          RecRad(i,irecurv)    = RadCurv(it)
          LonR  = LonSP(it)
          LatR =  LatSP(it)
          irecurv     = irecurv + 1
        end if
      end if 
    end do ;it
  
    RecNum(i) = irecurv


    if (   jFound .eq. 1 ) then
      NumStormRecurving = NumStormRecurving + 1
    end if
    
    if ( iPlot .eq. 1) then
      dateTit = "No recurvature found"
      if ( irecurv .ge.  1 ) then
        dateTit = " Recurves: " + cd_string( RecDate(i,0),"%H Z %d %c")
        do ir = 1,irecurv-1
          dateTit = dateTit + ", " + cd_string( RecDate(i,ir),"%H Z %d %c")
        end do
      end if
      
      res@tiMainString         =  sName + " " + cd_calendar(Filedates(istorm,0),-3)/1000000
      res@gsnCenterString      =   dateTit  
      plot = gsn_csm_map(wks,res)        

      do ir = 0,irecurv-1
        gsres@gsMarkerColor  =  colors(ir)
        dumG(ir) = gsn_add_polymarker(wks,plot,RecLon(i,ir),RecLat(i,ir),gsres)
        print ("recured at/on " + RecLon(i,ir) + " " + RecLat(i,ir) + " Date =  " + cd_string (RecDate (i,ir),"%H Z %d %c") + " " +  RecRad(i,ir))
      end do

      gsres@gsLineDashPattern = 0
      gsres@gsLineThicknessF = 2.0
      ; plot the unsmoothed track
      dum1 := gsn_add_polyline(wks,plot,LonSP,LatSP,gsres)
      gsres@gsLineDashPattern = 10
      ; plot the smoothed track
      ;dum2 := gsn_add_polyline(wks,plot,LL,LA,gsres)
      draw(plot)
      frame(wks)
      delete(plot)
 
      print ("Storm found to recurve " + RecNum(i) + " Times " )
    else
      print ( "No recuvature found  for this storm")
    end if
    print ( " " )
  end do


; 
;
; Write out to file

  fileName = "tcRecurves_" + basinToConsider +".nc"

  if (fileexists(fileName) ) then
    system( "rm " + fileName )
  end if
  setfileoption("nc", "format", "NETCDF4") 
  outFile = addfile( fileName, "c" )
  outFile->RecLat  = RecLat
  outFile->RecLon  = RecLon
  outFile->RecDate = RecDate
  outFile->RecNum  = RecNum
  outFile->RecNam  = RecNam
  outFile->RecSid  = RecSid




  print ( "--------------------------------------------------------")
  print ( "Total Number of Storms Processed = " + NSTORMS)
  print ( "Number of Storms that recurved   = " + NumStormRecurving)
  print ( "--------------------------------------------------------")
  print ( "--------------------------------------------------------")
  print ( "--------Normal Code Execution and End-------------------")
  print ( "--------Check output files------------------------------")
  print ( "--------------------------------------------------------")

 
end