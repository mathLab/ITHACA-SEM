<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=interpfield_4el.xml:fromfld=interpfield_4el.fld interpfield_1el.xml interpfield_1el.fld
    </parameters>
    <files>
        <file description="Mesh File 4 ele">interpfield_4el.xml</file>
        <file description="Mesh File 1 ele">interpfield_1el.xml</file>
        <file description="Field File 4 ele">interpfield_4el.fld</file>
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