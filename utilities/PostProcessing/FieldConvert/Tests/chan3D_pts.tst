<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable>FieldConvert</executable>
    <parameters>-e -m interppointdatatofld chan3D.xml chan3D_pts.pts chan3D_pts.fld</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D_pts.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="10">385025</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="10">175606</value>
        </metric>
    </metrics>
</test>

