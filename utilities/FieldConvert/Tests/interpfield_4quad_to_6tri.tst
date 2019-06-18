<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate field from 4 quads to 6 triangles </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=interpfield_4quad.xml:fromfld=interpfield_4quad.fld interpfield_6tri.xml interpfield_6tri.fld
    </parameters>
    <files>
        <file description="Mesh File 4 ele">interpfield_4quad.xml</file>
        <file description="Mesh File 1 ele">interpfield_6tri.xml</file>
        <file description="Field File 4 ele">interpfield_4quad.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.000683002</value>
            <value variable="y" tolerance="1e-6"> 0.00222435</value>
            <value variable="rho" tolerance="1e-6">0.0112035</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">0.0790222</value>
            <value variable="y" tolerance="1e-6">-0.233356</value>
            <value variable="rho" tolerance="1e-6">1.20061</value>
        </metric>
    </metrics>
</test>