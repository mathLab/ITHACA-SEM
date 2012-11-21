<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m3.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-12">1.17064e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.04578e-16</value>
            <value variable="v" tolerance="1e-12">2.27081e-16</value>
            <value variable="w" tolerance="1e-12">1.58207e-15</value>
	    <value variable="p" tolerance="1e-12">3.9746e-14</value>
        </metric>
    </metrics>
</test>

