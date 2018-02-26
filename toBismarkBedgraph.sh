CG='.CG.bedGraph'
gzcat $1 | awk '$4 == "CG" {print $1,'\t',$3,'\t',$3,'\t',$6}' | awk '$3+=1' > $1$CG

CHG='.CHG.bedGraph'
gzcat $1 | awk '$4 == "CHG" {print $1,'\t',$3,'\t',$3,'\t',$6}' | awk '$3+=1' > $1$CHG

CHH='.CHH.bedGraph'
gzcat $1 | awk '$4 == "CHH" {print $1,'\t',$3,'\t',$3,'\t',$6}' | awk '$3+=1' > $1$CHH

ALL='.bedGraph'
gzcat $1 | awk '{print $1,'\t',$3,'\t',$3,'\t',$6}' | awk '$3+=1' > $1$ALL
