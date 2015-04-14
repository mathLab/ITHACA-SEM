<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m3_dealiasing.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m3_dealiasing.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-12">3.57809e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.36905e-16</value>
            <value variable="v" tolerance="1e-12">1.00726e-16</value>
            <value variable="w" tolerance="1e-12">6.93889e-16</value>
	    <value variable="p" tolerance="1e-12">1.77636e-14</value>
        </metric>
    </metrics>
</test>

