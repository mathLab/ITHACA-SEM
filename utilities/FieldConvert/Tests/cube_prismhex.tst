<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Surface distance calculation on 3D hex/prism mesh</description>
    <executable>FieldConvert</executable>
    <parameters>-e -m surfdistance:bnd=0 cube_prismhex.xml out.fld</parameters>
    <files>
        <file description="Session File">cube_prismhex.xml</file>
    </files>
    <metrics>
        <metric type="l2" id="1">
            <value variable="dist" tolerance="1e-12">0.5</value>
        </metric>
    </metrics>
</test>
