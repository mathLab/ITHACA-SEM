

#!/bin/sh
for i in {57414..58413}
do
  FieldConvert cavity_poi_28_stage1.xml cavity_poi_28_stage1_$i.chk  mov_$i.vtu
done



