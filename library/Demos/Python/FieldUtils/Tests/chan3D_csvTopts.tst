<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Convert csv to pts </description>
    <executable python="true">chan3D_csvTopts.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D_pts.csv</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">1</value>
        </metric>
    </metrics>
</test>

