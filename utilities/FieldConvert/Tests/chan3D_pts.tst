<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable>FieldConvert</executable>
    <parameters>-e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_pts.csv out.csv</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_pts.csv</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.98603e-15</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">8</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-6">1.77636e-15</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">4</value>
        </metric>
    </metrics>
</test>

