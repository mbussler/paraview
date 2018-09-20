mencoder mf://*.tga -o %1.avi -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2

rem mencoder mf://*.tga -o %1.avi -mf fps=15:type=tga -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell

mkdir %1
move *.tga %1\