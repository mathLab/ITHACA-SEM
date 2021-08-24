<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output with n=10 physical points</description>
    <executable python="true">chan3D_tec_n10.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.63299</value>
            <value variable="y" tolerance="1e-6">1.63299</value>
            <value variable="z" tolerance="1e-6">1.63299</value>
            <value variable="u" tolerance="1e-6">2.06559</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">6.53197</value>
        </metric>
    </metrics>
</test>

