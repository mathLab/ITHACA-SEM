<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate a .csv onto a mesh </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e  --nparts 4 -m interppointdatatofld:frompts=chan_interpfield_nparts.csv chan_interpfield_nparts_xml:xml chan_interpfield_nparts.fld</parameters>
    <files>
        <file description="Session File">chan_interpfield_nparts_xml</file>
        <file description="Points File">chan_interpfield_nparts.csv</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">1.20866</value>
            <value variable="y" tolerance="1e-4">0.641719</value>
            <value variable="u" tolerance="1e-3">1.73853</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">-0.3</value>
            <value variable="y" tolerance="1e-4">0.309017</value>
            <value variable="u" tolerance="1e-3">3.23248</value>
        </metric>
    </metrics>
</test>