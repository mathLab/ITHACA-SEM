

#!/bin/sh
for i in {61771..62770}
do
  FieldConvert cavity_poi_130_stage1.xml cavity_poi_130_stage1_$i.chk  mov_$i.vtu
done



