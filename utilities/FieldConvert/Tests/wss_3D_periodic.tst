<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process wss with periodic boundary </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m wss:bnd=1 wss_3D_periodic.xml wss_3D_periodic.fld wss_3D_periodic-wss.fld </parameters>
    <files>
        <file description="Session File">wss_3D_periodic.xml</file>
        <file description="Field File">wss_3D_periodic.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x"       tolerance="1e-6">0.23206</value>
            <value variable="y"       tolerance="1e-6">0</value>
            <value variable="z"       tolerance="1e-7">0.023206</value>
            <value variable="Shear_x" tolerance="1e-9">0.00380052</value>
            <value variable="Shear_y" tolerance="1e-9">0.00152021</value>
        </metric>
    </metrics>
</test>