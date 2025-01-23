./sirecon /Users/hesper/local_data/20100521_Cho_lifeact_PPD/20100521_164304_NLSIM_9dirs_7phases_20sat-region3-1.mrc ./debug3.proc /Users/hesper/local_data/20100518_TIRF_Obj_Cell_coverslips/20100519_145033_PSF-9-256-20.otf -nphases 7 -ndirs 9 -ls 0.172 -na 1.2 -nimm 1.514 -angle0 -1.918 -wiener .01 -zoomfact 2 -nordersout 2 -driftfix -noff 1 -combineoff -background 385 -fix2Ddrift -nosuppress -fitonephase -forcemodamp 0.2


#-fixphasestep 

#./sirecon /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20100310_DronpaLifeAct/20100310_172313_NLSIM_7phases-9dirs-region5-2.mrc ./debug.proc /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20091027_203218_2D-green-TIRF_PSF.otf -nphases 7 -ndirs 9 -ls 0.166 -na 1.3 -nimm 1.514 -angle0 -1.918 -wiener .01 -zoomfact 2 -driftfix -noff 1 -background 280 -fitonephase -fix2Ddrift -fixphasestep

#./sirecon /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20100421_Dr40beads_PPD_GODCAT/20100421_133935_218uW-NLSIM_9angles-7phases-sat15-region4-1.mrc /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20100421_Dr40beads_PPD_GODCAT/20100421_133935_218uW-NLSIM_9angles-7phases-sat15-region4-1-nordersout3.proc /Users/hesper/local_data/20100406_1DSIM/20100402_163223_PSF-1-256.otf -nphases 9 -ndirs 7 -ls 0.168 -na 1.3 -nimm 1.514 -angle0 -1.921 -wiener .005 -zoomfact 3 -driftfix -noff 1 -background 470 -fitonephase -nordersout 3 -forcemodamp 0.5 0.2


# -fixphasestep  -fix2Ddrift  -noequalize


#./sirecon /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20100311_DronpaLifeAct/20100311_151625_NLSIM-7ph-9dirs-region3-1.mrc ./debug /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20091027_203218_2D-green-TIRF_PSF.otf -nphases 7 -ndirs 9 -ls 0.166 -na 1.3 -nimm 1.514 -angle0 -1.918 -wiener .01 -zoomfact 2 -driftfix -noff 1 -background 280 -fitonephase -fixphasestep  -fix2Ddrift -saveprefiltered /Volumes/gustafssonlab/MG2/hesper/OM_NLSI/20100311_DronpaLifeAct/20100311_151625_NLSIM-7ph-9dirs-region3-1-phicorr.sep
