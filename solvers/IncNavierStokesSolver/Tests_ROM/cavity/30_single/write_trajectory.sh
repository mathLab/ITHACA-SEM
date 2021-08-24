

#!/bin/sh
for i in {49744..52744}
do
  FieldConvert cavity_poi.xml cavity_poi_$i.chk  mov_$i.vtu
done



