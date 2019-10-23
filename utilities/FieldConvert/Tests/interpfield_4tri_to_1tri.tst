<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>  Interpolate field from 4 triangles to 1 triangle </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=interpfield_4tri.xml:fromfld=interpfield_4tri.fld interpfield_1tri.xml interpfield_1tri.fld
    </parameters>
    <files>
        <file description="Mesh File 4 ele">interpfield_4tri.xml</file>
        <file description="Mesh File 1 ele">interpfield_1tri.xml</file>
        <file description="Field File 4 ele">interpfield_4tri.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.21973</value>
            <value variable="y" tolerance="1e-6"> 0.0326321</value>
            <value variable="u" tolerance="1e-6">0.28525</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">4.48679</value>
            <value variable="y" tolerance="1e-6">0.199915</value>
            <value variable="u" tolerance="1e-6">1.01355</value>
        </metric>
    </metrics>
</test>