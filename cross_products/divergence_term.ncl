;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate vorticity flux divergence
; Sandro Lubis (sandro.lubis@pnnl.gov)
; February 2023
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function divergence_term(i_u,i_v,i_x)
begin

 ; Constants
 rad = get_d2r("float")
 re  = 6.371e+06
 g   = 9.80665 
 pi  = atan(1.0)*4.

 lon = dble2flt(i_u&lon)
 lat = dble2flt(i_u&lat)

 print("Compute horizontal gradient of wind " )

 dUdLon = center_finite_diff_n(i_u, lon * rad, True, 0, 2) / conform(i_u, re * cos(lat * rad), 1)
 dVdLat = center_finite_diff_n(i_v, lat * rad, False, 0, 1) / re

 copy_VarCoords(i_u,dUdLon)
 copy_VarCoords(i_u,dVdLat)

  
 print("Compute q divergence ")

 qDivx  = i_u
 qDivy  = i_u
 qDiv   = i_u
 qDivx  = (/-i_x*dUdLon/)
 qDivy  = (/-i_x*dVdLat/)
 qDiv   = qDivx +qDivy
 
 ;printMinMax(qDiv, True) 

 ;print("Compute vorticity divergence Directly. ")

 ;div       = uv2dv_cfd(i_u,i_v,lat,lon,3)
 ;copy_VarCoords(i_u,div)

 ;vortDivergence  = i_x
 ;vortDivergence  = (/-i_x*div/)

 ;printMinMax(vortDivergence, True)

 return((/qDivx,qDivy/))

end 

