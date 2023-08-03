; Calculate mass-weighted column integrated moisture budget
; Author: Dr. Sandro Lubis (Jan 2020)
;         Moisture Budget as in Lubis et al., (2023) and Lubis and Respati (2020)
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
 psfc = 100000. ; in Pa
 ptop = 10000.
 g=9.81
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

 do year=1979,2020
     print(year)


print("Reading the input file...")

  f1      = addfile("../../U/daily/u."+year+".nc","r")
  f2      = addfile("../../V/daily/v."+year+".nc","r")
  f3      = addfile("../../W/daily/w."+year+".nc","r")
  f4      = addfile("../../SH/daily/q."+year+".nc","r")

time 		= f1->time
ntime           = dimsizes(time)
lev 		= dble2flt(f1->plev({psfc:ptop}))* 1. ; hPa to Pa
lat 		= dble2flt(f1->lat(::-1))
lon	 	= dble2flt(f1->lon)

u = f1->var131(:,{psfc:ptop},::-1,:)
v = f2->var132(:,{psfc:ptop},::-1,:)
w = f3->var135(:,{psfc:ptop},::-1,:)
q = f4->var133(:,{psfc:ptop},::-1,:)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating temporal tendency of moisture...")

  q1 = q
  q1(1:ntime-1,:,:,:) = (/q(0:ntime-2,:,:,:)/)

  dq_dt = (/(q-q1)/86400.0/)
  copy_VarCoords(q, dq_dt)

dq_dt@long_name = "Temporal tendency of moisture"
dq_dt@units = "kg/kg s**-1"

printVarSummary(dq_dt)
printMinMax(dq_dt, False)

delete(q1)

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
print("Calculating vertical integrals...")

psfc = max(lev)
ptop = min(lev)
dp = dpres_plevel(lev, psfc, ptop, 0)
printVarSummary(dp)
print(dp)


vi_q = (1.0 / g) * wgt_vertical_n(q, conform(q, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_q)
vi_q@long_name = "Vertical integral of moisture (a.k.a column water vapour)"
vi_q@units = "kg m**-2 (mm)"
printVarSummary(vi_q)
printMinMax(vi_q, False)


vi_dq_dt = (1.0 / g) * wgt_vertical_n(dq_dt, conform(dq_dt, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_dq_dt)
vi_dq_dt@long_name = "Vertical integral of temporal tendency of moisture"
vi_dq_dt@units = "kg m**-2 s**-1"
printVarSummary(vi_dq_dt)
printMinMax(vi_dq_dt, False)

vi_adv_q = (1.0 / g) * wgt_vertical_n(adv_q, conform(adv_q, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_adv_q)
vi_adv_q@long_name = "Vertical integral of advection of moisture"
vi_adv_q@units = "kg m**-2 s**-1"
printVarSummary(vi_adv_q)
printMinMax(vi_adv_q, False)

vi_conv_q = (1.0 / g) * wgt_vertical_n(conv_q, conform(conv_q, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_conv_q)
vi_conv_q@long_name = "Vertical integral of convergence of moisture"
vi_conv_q@units = "kg m**-2 s**-1"
printVarSummary(vi_conv_q)
printMinMax(vi_conv_q, False)

vi_conv_q_zonal = (1.0 / g) * wgt_vertical_n(conv_q_zonal, conform(conv_q_zonal, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_conv_q_zonal)
vi_conv_q_zonal@long_name = "Vertical integral of zonal convergence of moisture"
vi_conv_q_zonal@units = "kg m**-2 s**-1"
printVarSummary(vi_conv_q_zonal)
printMinMax(vi_conv_q_zonal, False)

vi_conv_q_meridional = (1.0 / g) * wgt_vertical_n(conv_q_meridional, conform(conv_q_meridional, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_conv_q_meridional)
vi_conv_q_meridional@long_name = "Vertical integral of meridional convergence of moisture"
vi_conv_q_meridional@units = "kg m**-2 s**-1"
printMinMax(vi_conv_q_meridional, False)

vi_mfc = (1.0 / g) * wgt_vertical_n(mfc, conform(mfc, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_mfc)
vi_mfc@long_name = "Vertical integral of moisture flux convergence"
vi_mfc@units = "kg m**-2 s**-1"
printVarSummary(vi_mfc)
printMinMax(vi_mfc, False)

vi_mfc_vertical = (1.0 / g) * wgt_vertical_n(mfc_vertical, conform(mfc_vertical, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_mfc_vertical)
vi_mfc_vertical@long_name = "Vertical integral of vertical convergence of moisture flux"
vi_mfc_vertical@units = "kg m**-2 s**-1"
printVarSummary(vi_mfc_vertical)
printMinMax(vi_mfc_vertical, False)

vi_adv_q_vertical= (1.0 / g) * wgt_vertical_n(adv_q_vertical, conform(adv_q_vertical, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_adv_q_vertical)
vi_adv_q_vertical@long_name = "Vertical integral of adv_q_vertical"
vi_adv_q_vertical@units = "kg m**-2 s**-1"
printVarSummary(vi_adv_q_vertical)
printMinMax(vi_adv_q_vertical, False) 

vi_conv_q_vertical= (1.0 / g) * wgt_vertical_n(conv_q_vertical, conform(conv_q_vertical, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_conv_q_vertical)
vi_conv_q_vertical@long_name = "Vertical integral of conv_q_vertical"
vi_conv_q_vertical@units = "kg m**-2 s**-1"
printVarSummary(vi_conv_q_vertical)
printMinMax(vi_conv_q_vertical, False) 

vi_q2 = (1.0 / g) * wgt_vertical_n(q2, conform(q2, dp, 1), 1, 1)
copy_VarCoords(q(:, 0, :, :), vi_q2)
vi_q2@long_name = "Vertical integral of moisture sink/source"
vi_q2@units = "kg m**-2 s**-1"
printVarSummary(vi_q2)
printMinMax(vi_q2, False)


;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Save output

setfileoption("nc","Format","LargeFile")
system( "rm " + "moisture_budget.vint."+year+".nc")
output = addfile("moisture_budget.vint."+year+".nc","c")
output->dq_dt  = vi_dq_dt
output->adv_q  = vi_adv_q
output->conv_q = vi_conv_q
output->conv_q_zonal = vi_conv_q_zonal
output->conv_q_meridional = vi_conv_q_meridional
output->mfc        = vi_mfc
output->mfc_ver    = vi_mfc_vertical
output->adv_q_vertical   = vi_adv_q_vertical
output->conv_q_vertical   = vi_conv_q_vertical
output->q2         = vi_q2
output->q           = vi_q

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; delete temporary files

	delete(time)
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

	delete(vi_dq_dt)
	delete(vi_adv_q)
	delete(vi_conv_q)
	delete(vi_conv_q_zonal)
	delete(vi_conv_q_meridional)
	delete(vi_mfc)
	delete(vi_mfc_vertical)
	delete(vi_adv_q_vertical)
	delete(vi_conv_q_vertical)
	delete(vi_q2)
	delete(vi_q)
end do



end

