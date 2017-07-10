<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Project point data to a fld file </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e --noequispaced -m pointdatatofld outflow.pts outflow.xml outflow.fld</parameters>
    <files>
        <file description="Session File">outflow.pts</file>
        <file description="Session File">outflow.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-6"> 0.0959102</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="p" tolerance="1e-6"> -0.0178023</value>
        </metric>
    </metrics>
</test>

