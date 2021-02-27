; This program reads the recurvature data produced by: recurve_ibtracs.ncl
;
;
;
;---------------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


begin

  fileName = "tcRecurves_WP.nc"

  inFile = addfile( fileName, "r" )

  ; names of each storm 
  ; not all of the storms recurv
  RecNam  = inFile->RecNam
  ; ibtracs storm id
  RecSid= inFile->RecSid
  
  ; Number of times each storm recurved 0 to 4
  RecNum  = inFile->RecNum

  ; Lat/Lon of recurvature (default = missing)
  RecLat  = inFile->RecLat
  RecLon  = inFile->RecLon

  ; date of each recurvature if it happened
  RecDate = inFile->RecDate



  NSTORMS = dimsizes(RecNam)

  do i = 0, NSTORMS-1

    if ( RecNum(i) .eq. 0 ) then
      print ( "Did not recurve: Storm ID, Name = " + RecSid(i) + " " + RecNam(i) )
    else
      print ( "Recurved " + RecNum(i) + " times: Storm ID, Name = " + RecSid(i) + " " + RecNam(i) )
      do ir = 0,RecNum(i)-1
        print ("recured at/on " + RecLon(i,ir) + " " + RecLat(i,ir) + " Date =  " + cd_string (RecDate (i,ir),"%H Z %d %c") )
      end do

      
    end if

    print ( " " ) 
  end do



end







