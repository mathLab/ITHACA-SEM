<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_robin_m3.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_robin_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">1.55986e-13</value>
	    <value variable="p" tolerance="1e-12">1.86517e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">.15815e-14</value>
	    <value variable="p" tolerance="1e-12">4.97947e-15</value>
        </metric>
    </metrics>
</test>

