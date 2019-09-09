<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp 2D field to line </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interppoints:fromxml=interpfieldline.xml:fromfld=interpfieldline.fld:line=10,0.024,0.0,0.16,0.0 interpfieldline.dat</parameters>
    <files>
        <file description="Session File">interpfieldline.xml</file>
        <file description="Session File">interpfieldline.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.101724</value>
            <value variable="y" tolerance="1e-4">0.0</value>
            <value variable="Shear_mag" tolerance="1e-4">4.68391</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">0.16</value>
            <value variable="y" tolerance="1e-4">0.0</value>
            <value variable="Shear_mag" tolerance="1e-4">7.74933</value>
        </metric>
    </metrics>
</test>