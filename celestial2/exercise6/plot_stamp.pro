pro plot_stamp,filename

hostname=hostname()
whoami=mywhoami()
date=systime()
cd,curr=cwd
space='  '
header1=cwd+'/'+filename
;header2=whoami+'@'+hostname+space+date
header2=whoami+space+date



xyouts,0.01,0.005,header1,/nor,chars=.7
xyouts,.99,0.005,header2,/nor,chars=.7,ali=1
end
