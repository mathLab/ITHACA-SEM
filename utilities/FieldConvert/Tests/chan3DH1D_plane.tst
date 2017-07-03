<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D plane output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -m homplane:planeid=4 -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_plane.fld</parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
	<file description="Session File">chan3DH1D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.57735</value>
            <value variable="y" tolerance="1e-6">0.57735</value>
            <value variable="u" tolerance="1e-6">0.182574</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.1547</value>
        </metric>
    </metrics>
</test>

