load "divergence_term.ncl"

begin

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Constants
 rad = get_d2r("float")
 re  = 6.371e+06
 g   = 9.80665 
 pi  = atan(1.0)*4.
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

 do year=1979,2021
     print(year)


print("Reading the input file...")

  f1      = addfile("./climo/climo.nc","r")
  f2      = addfile("./eddy/eddy_bf."+year+".nc","r")
  f3      = addfile("./eddy/eddy_hf."+year+".nc","r")

  uc 	  = lonFlip(f1->var131(:,0,::-1,:))
  vc 	  = lonFlip(f1->var132(:,0,::-1,:))
  xc 	  = lonFlip(f1->var133(:,0,::-1,:))

  printVarSummary(uc)
  printVarSummary(vc)
  printVarSummary(xc)


  ub 	 = lonFlip(f2->u_bf(:,::-1,:))
  uh 	 = lonFlip(f3->u_hf(:,::-1,:))

  vb 	 = lonFlip(f2->v_bf(:,::-1,:))
  vh 	 = lonFlip(f3->v_hf(:,::-1,:))

  xb 	 = lonFlip(f2->q_bf(:,::-1,:))
  xh 	 = lonFlip(f3->q_hf(:,::-1,:))

  printVarSummary(ub)
  printVarSummary(uh)


 time 		 	= f2->time
 ntime			= dimsizes(time)
 lat 		 	= dble2flt(ub&lat)
 lon	 	   	= dble2flt(ub&lon)

 ; Coriolis parameter
  f    =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat)
  ftmp = conform_dims(dimsizes(ub),f,1)


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating terms...")


;; interaction v_clim and u_bf


ub_xc = divergence_term (ub, vb, xc)
ub_xh = divergence_term (ub, vb, xh)
ub_xb = divergence_term (ub, vb, xb)

uc_xb = divergence_term (uc, vc, xb)
uc_xh = divergence_term (uc, vc, xh)

uh_xh = divergence_term (uh, vh, xh)
uh_xb = divergence_term (uh, vh, xb)
uh_xc = divergence_term (uh, vh, xc)

copy_VarCoords(ub,ub_xc(0,:,:,:))
copy_VarCoords(ub,ub_xb(0,:,:,:))
copy_VarCoords(ub,ub_xh(0,:,:,:))

copy_VarCoords(ub,uc_xb(0,:,:,:))
copy_VarCoords(ub,uc_xh(0,:,:,:))

copy_VarCoords(ub,uh_xh(0,:,:,:))
copy_VarCoords(ub,uh_xb(0,:,:,:))
copy_VarCoords(ub,uh_xc(0,:,:,:))



printVarSummary(ub_xc)


;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Save output

setfileoption("nc","Format","LargeFile")
system( "rm " + "./output/xqbudget."+year+".nc")
output = addfile("./output/xqbudget."+year+".nc","c")
output->conv_qcub  = ub_xc(0,:,:,:)
output->conv_qbub  = ub_xb(0,:,:,:)
output->conv_qhub  = ub_xh(0,:,:,:)
output->conv_qbuc  = uc_xb(0,:,:,:)
output->conv_qhuc  = uc_xh(0,:,:,:)
output->conv_qcuh  = uh_xc(0,:,:,:)
output->conv_qbuh  = uh_xb(0,:,:,:)
output->conv_qhuh  = uh_xh(0,:,:,:)

system( "rm " + "./output/yqbudget."+year+".nc")
output2 = addfile("./output/yqbudget."+year+".nc","c")
output2->conv_qcub  = ub_xc(1,:,:,:)
output2->conv_qbub  = ub_xb(1,:,:,:)
output2->conv_qhub  = ub_xh(1,:,:,:)
output2->conv_qbuc  = uc_xb(1,:,:,:)
output2->conv_qhuc  = uc_xh(1,:,:,:)
output2->conv_qcuh  = uh_xc(1,:,:,:)
output2->conv_qbuh  = uh_xb(1,:,:,:)
output2->conv_qhuh  = uh_xh(1,:,:,:)


;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; delete temporary files

	delete(time)
	delete(ntime)
	delete(ub)
	delete(uh)
	delete(vb)
	delete(vh)
	delete(xb)
	delete(xh)
	delete(ftmp)
	delete(f)

	delete(ub_xc)
	delete(ub_xb)
	delete(ub_xh)

	delete(uc_xb)
	delete(uc_xh)

	delete(uh_xc)
	delete(uh_xb)
	delete(uh_xh)



	


end do



end

