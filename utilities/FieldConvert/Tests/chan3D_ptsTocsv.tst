<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert pts to csv </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e chan3D_pts.pts chan3D_pts.csv</parameters>
    <files>
        <file description="Session File">chan3D_pts.pts</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">1</value>
            <value variable="p" tolerance="1e-6">136931</value>
        </metric>
    </metrics>
</test>

