<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a .fld file to .pts </description>
    <executable python="true">chan3D_pts.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
        <file description="Session File">chan3D_pts.pts</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">1.0</value>
            <value variable="y" tolerance="1e-4">1.0</value>
            <value variable="z" tolerance="1e-4">1.0</value>
            <value variable="u" tolerance="1e-6">1.98603e-15</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">2.82843</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="x" tolerance="1e-4">1.0</value>
            <value variable="y" tolerance="1e-4">1.0</value>
            <value variable="z" tolerance="1e-4">1.0</value>
            <value variable="u" tolerance="1e-6">1.77636e-15</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">4</value>
        </metric>
    </metrics>
</test>
