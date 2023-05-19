;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate vorticity flux divergence
; Sandro Lubis (sandro.lubis@pnnl.gov)
; February 2023
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function  vort_flux_divergence(i_u,i_v,i_x)
begin

 ; Constants
 rad = get_d2r("float")
 re  = 6.371e+06
 g   = 9.80665 
 pi  = atan(1.0)*4.

 lon = dble2flt(i_u&lon)
 lat = dble2flt(i_u&lat)

 print("Compute horizontal gradient of vorticity. " )

 dVortdLon = center_finite_diff_n(i_x, lon * rad, True, 0, 2) / conform(i_u, re * cos(lat * rad), 1)
 dVortdLat = center_finite_diff_n(i_x, lat * rad, False, 0, 1) / re

 copy_VarCoords(i_u,dVortdLat)
 copy_VarCoords(i_u,dVortdLon)
  
 print("Compute divergence of wind. ")

 div       = uv2dv_cfd(i_u,i_v,lat,lon,3)
 copy_VarCoords(i_u,div)
  
 print("Compute vorticity advection. ")

 vortAdvectionx  = i_u
 vortAdvectiony  = i_u
 vortAdvection   = i_u
 vortAdvectionx  = (/-i_u*dVortdLon/)
 vortAdvectiony  = (/-i_v*dVortdLat/)
 vortAdvection   = vortAdvectionx + vortAdvectiony
  
 print("Compute vorticity divergence. ")
  
 vortDivergence  = i_x
 vortDivergence  = (/-i_x*div/)

 tot=i_u
 tot=vortAdvection+vortDivergence

 ;printMinMax(vortAdvection, True)
 ;printMinMax(vortDivergence, True)
 ;printMinMax(tot, True)

 return(tot)

end 

