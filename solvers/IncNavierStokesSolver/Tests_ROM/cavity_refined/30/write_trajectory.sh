

#!/bin/sh
for i in {1001..2001}
do
  FieldConvert cavity_poi_30.xml cavity_poi_30_$i.chk  mov_$i.vtu
done



