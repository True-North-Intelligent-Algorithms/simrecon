./sirecon /Users/shaol/SIMrecon/testdata/20090724_YG100_2DSIM_200ms_exp_3phases_3angles-2-512.mrc /Users/shaol/SIMrecon/testdata/2D.proc /Users/shaol/SIMrecon/testdata/20090724_YG100_2DPSF_Zeisswater_500ms_Conv-512.otf -nphases 3 -ndirs 3 -ls 0.224 -na 1.2 -nimm 1.333 -angle0 -1.271 -wiener .005 -zoomfact 2
