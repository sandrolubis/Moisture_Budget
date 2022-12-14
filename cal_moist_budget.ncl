; Author: Dr. Sandro Lubis (Jan 2020)
;         Moisture Budget as in Lubis and Respati (2020)
;         Lubis, SW, Respati, MR. Impacts of convectively coupled 
;         equatorial waves on rainfall extremes in Java, Indonesia. 
;         Int J Climatol. 2021; 41: 2418– 2440
; Contact: slubis.geomar@gmail.com
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; input [time][lev][lat][lon]
; time    - "seconds "
; p       - Pressure [Pa]
; u,v     - zonal, meridional wind components[m/s]
; q       - specific humidity [kg/kg]
; T       - temperature [K]  or [C]
; omega   - vertical velocity [Pa/s]
; latitude is in ascending order
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Constants
 rad = get_d2r("float")
 re = 6.371e+06
 psfc = 1000.
 ptop = 300.
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

 do year=2001,2001
     print(year)


print("Reading the input file...")

f = addfile("input/input_"+year+".nc", "r")

time 		= f->time
lev 		= int2flt(f->level({psfc:ptop})) * 100. ; hPa to Pa
lat 		= f->latitude(::-1)
lon	 	= f->longitude

q = short2flt(f->q(:,{psfc:ptop},::-1,:))
u = short2flt(f->u(:,{psfc:ptop},::-1,:))
v = short2flt(f->v(:,{psfc:ptop},::-1,:))
w = short2flt(f->w(:,{psfc:ptop},::-1,:))

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating temporal tendency of moisture...")

dq_dt = center_finite_diff_n(q, time, False, 0, 0) / 3600.0
copy_VarCoords(q, dq_dt)

dq_dt@long_name = "Temporal tendency of moisture"
dq_dt@units = "kg/kg s**-1"

printVarSummary(dq_dt)
printMinMax(dq_dt, False)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating advection of moisture...")


dq_dx = center_finite_diff_n(q, lon * rad, True, 0, 3) / conform(q, re * cos(lat * rad), 2)
dq_dy = center_finite_diff_n(q, lat * rad, False, 0, 2) / re
dq_dp = center_finite_diff_n(q, lev, False, 0, 1) 

adv_q_zonal 		= -u * dq_dx
adv_q_meridional 	= -v * dq_dy
adv_q_vertical 		= -w * dq_dp

adv_q = adv_q_zonal + adv_q_meridional
copy_VarCoords(q, adv_q)

adv_q@long_name = "Advection of moisture"
adv_q@units = "kg/kg s**-1"

printVarSummary(adv_q)
printMinMax(adv_q, False)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating moisture loss (gain) due to the condensational heating")
q2 = dq_dt + (u * dq_dx) + (v * dq_dy) + (w * dq_dp)

copy_VarCoords(q, q2)
printVarSummary(q2)
printMinMax(q2, False)
q2@long_name = "moisture sink"
q2@units = "kg/kg s**-1"


delete(dq_dx)
delete(dq_dy)
delete(dq_dp)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating convergence of moisture...")

du_dx = center_finite_diff_n(u, lon * rad, True, 0, 3) / conform(q, re * cos(lat * rad), 2)
dv_dy = center_finite_diff_n(v, lat * rad, False, 0, 2) / re
dw_dp = center_finite_diff_n(w, lev, False, 0, 1) 

conv_q_zonal = -q * du_dx
conv_q_meridional = -q * dv_dy
conv_q_vertical = -q * dw_dp

conv_q = conv_q_zonal + conv_q_meridional
copy_VarCoords(q, conv_q)

conv_q@long_name = "Convergence of moisture"
conv_q@units = "kg/kg s**-1"

printVarSummary(conv_q)
printMinMax(conv_q, False)

delete(du_dx)
delete(dv_dy)
delete(dw_dp)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating horizontal convergence of moisture flux...")

mfc = adv_q + conv_q
copy_VarCoords(q, mfc)

mfc@long_name = "Convergence of moisture flux"
mfc@units = "kg/kg s**-1"

printVarSummary(mfc)
printMinMax(mfc, False)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating vertical convergence of moisture flux...")

mfc_vertical = adv_q_vertical + conv_q_vertical
copy_VarCoords(q, mfc_vertical)

mfc_vertical@long_name = "Vertical convergence of moisture flux"
mfc_vertical@units = "kg/kg s**-1"

printVarSummary(mfc_vertical)
;printMinMax(adv_q_vertical, False)
;printMinMax(conv_q_vertical, False)


;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Save output

setfileoption("nc","Format","LargeFile")
system( "rm " + "moisture_budget."+year+".nc")
output = addfile("moisture_budget."+year+".nc","c")
output->dq_dt  = dq_dt
output->adv_q  = adv_q
output->conv_q = conv_q
output->mfc    = mfc
output->mfc_ver    = mfc_vertical
output->q2   = q2

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; delete temporary files

	delete(f)
	delete(q)
	delete(u)
	delete(v)
	delete(w)

	delete(dq_dt)
	delete(adv_q_zonal)
	delete(adv_q_meridional)
	delete(adv_q_vertical)
	delete(adv_q)

	delete(conv_q_zonal)
	delete(conv_q_meridional)
	delete(conv_q_vertical)
	delete(conv_q)

	delete(mfc)
	delete(mfc_vertical)
	delete(q2)

end do



end

