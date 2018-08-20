<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interppointdatatofld:frompts=chan3D_pts1D.pts:interpcoord=1 chan3D.xml chan3D_pts.fld</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D_pts1D.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.63299</value>
            <value variable="y" tolerance="1e-6">1.63299</value>
            <value variable="z" tolerance="1e-6">1.63299</value>
            <value variable="p" tolerance="20">  384817</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">1</value>
            <value variable="p" tolerance="10">  170000</value>
        </metric>
    </metrics>
</test>
