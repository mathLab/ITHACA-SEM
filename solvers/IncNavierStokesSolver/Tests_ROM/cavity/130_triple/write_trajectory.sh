

#!/bin/sh
for i in {61771..62770}
do
  FieldConvert cavity_poi_130_stage4.xml cavity_poi_130_stage4_$i.chk  mov4_$i.vtu
done



