<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> take a list of points at quadrature points and produce a fld file </description>
    <executable>FieldConvert</executable>
    <parameters>-m pointdatatofld -n 5  -e ceiling_velocity.pts ceiling.xml out.fld </parameters>
    <files>
        <file description="Session File">ceiling_velocity.pts</file>
        <file description="Session File">ceiling.xml</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-4">3.09428</value>
            <value variable="v" tolerance="1e-4">0.073833</value>
            <value variable="w" tolerance="1e-4">0.11834</value>
        </metric>
    </metrics>
</test>

