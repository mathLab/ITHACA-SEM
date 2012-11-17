<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8.xml</parameters>
    <files>
        <file description="Session File">Hex_channel_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">8.2406e-16</value>
            <value variable="v" tolerance="1e-12">8.90602e-16</value>
            <value variable="w" tolerance="1e-12">3.35358e-15</value>
	    <value variable="p" tolerance="1e-12">1.08744e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.02458e-15</value>
            <value variable="v" tolerance="1e-12">4.58854e-15</value>
            <value variable="w" tolerance="1e-12">1.77636e-13</value>
	    <value variable="p" tolerance="1e-12">5.15477e-13</value>
        </metric>
    </metrics>
</test>

