

#!/bin/sh
for i in {56765..57764}
do
  FieldConvert cavity_poi.xml cavity_poi_$i.chk  mov_$i.vtu
done



