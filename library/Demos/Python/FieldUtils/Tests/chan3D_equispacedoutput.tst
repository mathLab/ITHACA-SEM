<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert a field into a equi-spaced tecplot file</description>
    <executable python="true">chan3D_equispacedoutput.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.707107</value>
            <value variable="y" tolerance="1e-4">0.707107</value>
            <value variable="z" tolerance="1e-4">0.707107</value>
            <value variable="u" tolerance="1e-4">0.65192</value>
            <value variable="v" tolerance="1e-4">0</value>
            <value variable="w" tolerance="1e-4">0</value>
            <value variable="p" tolerance="1e-4">2.44949</value>
        </metric>
    </metrics>
</test>

