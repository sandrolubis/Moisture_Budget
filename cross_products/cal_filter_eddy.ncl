;begin
;=====================================================================
f1 = addfile ("u700.anom.nc","r")    ;u
f2 = addfile ("v700.anom.nc","r")    ;v
f3 = addfile ("q700.anom.nc","r")    ;SH

;=====================================================================

ga=9.80665
pi = atan(1.0)*4.

u   = f1->var131  (:,0,:,:)
v   = f2->var132  (:,0,:,:)
q   = f3->var133  (:,0,:,:)


   printVarSummary(u)
   printVarSummary(v)
   printVarSummary(q)

; Bandpass filter:
;===========================================================================

    nwt = 141
    fca = 1./90.
    fcb = 1./20.
    ihp = 2           ;ihp = 0; high-pass ihp = 1; band-pass ihp = 2 
    nsigma = 1.
    wgt = filwgts_lanczos (nwt, ihp, fca, fcb, nsigma)

u_bf   = wgt_runave_n_Wrap ( u, wgt, 0, 0 )
v_bf   = wgt_runave_n_Wrap ( v, wgt, 0, 0 )
q_bf   = wgt_runave_n_Wrap ( q, wgt, 0, 0 )


printVarSummary (u_bf)
printVarSummary (v_bf)
printVarSummary (q_bf)

; Highpass filter:
;===========================================================================

    nwt2 = 141
    fca2 = 1./20.
    fcb2 = -999.
    ihp2 = 1           ;ihp = 0; high-pass ihp = 1; band-pass ihp = 2 
    nsigma2 = 1.
    wgt2 = filwgts_lanczos (nwt2, ihp2, fca2, fcb2, nsigma2)

u_hf   = wgt_runave_n_Wrap ( u, wgt2, 0, 0 )
v_hf   = wgt_runave_n_Wrap ( v, wgt2, 0, 0 )
q_hf   = wgt_runave_n_Wrap ( q, wgt2, 0, 0 )


printVarSummary (u_hf)
printVarSummary (v_hf)
printVarSummary (q_hf)


; Lowpass filter:
;===========================================================================

u_lf   = u - (u_hf+u_bf)
v_lf   = v - (v_hf+v_bf)
q_lf   = q - (q_hf+q_bf)

copy_VarCoords(u, u_lf)
copy_VarCoords(v, v_lf)
copy_VarCoords(q, q_lf)


;Save
;*************************************
  setfileoption( "nc", "Format", "LargeFile" )
  system   ("/bin/rm ./eddy/eddy_bf.nc")        ; remove any pre-existing file
  fout     = addfile("./eddy/eddy_bf.nc","c")      ; open output netCDF file
  fout->u_bf= u_bf
  fout->v_bf= v_bf
  fout->q_bf= q_bf

  system   ("/bin/rm ./eddy/eddy_hf.nc")        ; remove any pre-existing file
  fout2     = addfile("./eddy/eddy_hf.nc","c")      ; open output netCDF file
  fout2->u_hf= u_hf
  fout2->v_hf= v_hf
  fout2->q_hf= q_hf

  system   ("/bin/rm ./eddy/eddy_lf.nc")        ; remove any pre-existing file
  fout3     = addfile("./eddy/eddy_lf.nc","c")      ; open output netCDF file
  fout3->u_lf= u_lf
  fout3->v_lf= v_lf
  fout3->q_lf= q_lf



;end
