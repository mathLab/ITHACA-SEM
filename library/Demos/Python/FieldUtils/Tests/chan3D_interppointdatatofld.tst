<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .pts file to .fld </description>
    <executable python="true">chan3D_interppointdatatofld.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D_pts.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.63299</value>
            <value variable="y" tolerance="1e-6">1.63299</value>
            <value variable="z" tolerance="1e-6">1.63299</value>
            <value variable="p" tolerance="10">385062</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">1</value>
            <value variable="p" tolerance="20">170000</value>
        </metric>
    </metrics>
</test>

