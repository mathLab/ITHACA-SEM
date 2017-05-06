<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> take a list of points at quadrature points and produce a fld file </description>
    <executable>FieldConvert</executable>
    <parameters> -f -m pointdatatofld:frompts=ceiling_velocity.pts -n 5 -e ceiling.xml out.fld </parameters>
    <files>
        <file description="Session File">ceiling_velocity.pts</file>
        <file description="Session File">ceiling.xml</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">22.4127</value>
            <value variable="y" tolerance="1e-6">0.807287</value>
            <value variable="z" tolerance="1e-6">8.27412</value>
            <value variable="u" tolerance="1e-4">3.09428</value>
            <value variable="v" tolerance="1e-4">0.073833</value>
            <value variable="w" tolerance="1e-4">0.11834</value>
        </metric>
    </metrics>
</test>

