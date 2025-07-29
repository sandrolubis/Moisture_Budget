begin

;=======================
; Constants
;=======================
rad = get_d2r("float")
re  = 6.371e6  ; Earth radius in meters

;=======================
; Read files
;=======================
f1 = addfile("u850.nc", "r")
f2 = addfile("v850.nc", "r")
f3 = addfile("q850.nc", "r")

time = f1->time
lat  = f1->lat
lon  = f1->lon

; Read single pressure level: index 0 (850 hPa)
u = f1->var131(:,0,:,:)  ; [time, lat, lon]
v = f2->var132(:,0,:,:)  ; [time, lat, lon]
q = f3->var133(:,0,:,:)  ; [time, lat, lon]

print("Variables loaded:")
printVarSummary(u)
printVarSummary(v)
printVarSummary(q)

;=======================
; Compute horizontal moisture advection
; adv_q = - (u ∂q/∂x + v ∂q/∂y)
;=======================
dq_dx = center_finite_diff_n(q, lon * rad, True, 0, 2) / conform(q, re * cos(lat * rad), 1)
dq_dy = center_finite_diff_n(q, lat * rad, False, 0, 1) / re

adv_q_zonal      = -u * dq_dx
adv_q_meridional = -v * dq_dy
adv_q = adv_q_zonal + adv_q_meridional
copy_VarCoords(q, adv_q)

adv_q@long_name = "Horizontal moisture advection"
adv_q@units     = "kg/kg/s"

;=======================
; Compute convergence of moisture from wind divergence
; conv_q = - (q ∂u/∂x + q ∂v/∂y)
;=======================
du_dx = center_finite_diff_n(u, lon * rad, True, 0, 2) / conform(u, re * cos(lat * rad), 1)
dv_dy = center_finite_diff_n(v, lat * rad, False, 0, 1) / re

conv_q_zonal      = -q * du_dx
conv_q_meridional = -q * dv_dy
conv_q = conv_q_zonal + conv_q_meridional
copy_VarCoords(q, conv_q)

conv_q@long_name = "Convergence of moisture"
conv_q@units     = "kg/kg/s"

;=======================
; Compute total moisture flux convergence
; MFC = adv_q + conv_q
;=======================
mfc = adv_q + conv_q
copy_VarCoords(q, mfc)
mfc@long_name = "Moisture Flux Convergence"
mfc@units     = "kg/kg/s"

;=======================
; Save to NetCDF
;=======================
setfileoption("nc", "Format", "LargeFile")
system("rm -f mfc_850.nc")
fout = addfile("mfc_850.nc", "c")

fout->mfc               = mfc
;fout->adv_q             = adv_q
;fout->conv_q            = conv_q
;fout->adv_q_zonal       = adv_q_zonal
;fout->adv_q_meridional  = adv_q_meridional
;fout->conv_q_zonal      = conv_q_zonal
;fout->conv_q_meridional = conv_q_meridional

print(" Moisture flux convergence saved to mfc_850.nc")

end
