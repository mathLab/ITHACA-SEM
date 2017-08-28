<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .dat </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld:topts=chan3D_probe.pts chan3D_probe.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_probe.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.5</value>
            <value variable="y" tolerance="1e-6">0.5</value>
            <value variable="z" tolerance="1e-6">0.5</value>
            <value variable="u" tolerance="1e-6">0.75</value>
            <value variable="v" tolerance="1e-5">0</value>
            <value variable="w" tolerance="1e-5">0</value>
            <value variable="p" tolerance="1e-5">2.23607</value>
        </metric>
    </metrics>
</test>

