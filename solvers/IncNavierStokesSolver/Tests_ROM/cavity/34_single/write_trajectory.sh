

#!/bin/sh
for i in {57767..58766}
do
  FieldConvert cavity_poi_34_stage1.xml cavity_poi_34_stage1_$i.chk  mov_$i.vtu
done



