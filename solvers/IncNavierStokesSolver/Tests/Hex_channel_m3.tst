<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m3.xml</parameters>
    <files>
        <file description="Session File">Hex_channel_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">1.90094e-16</value>
            <value variable="w" tolerance="1e-12">7.91195e-16</value>
	    <value variable="p" tolerance="1e-12">2.01754e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.7391e-16</value>
            <value variable="v" tolerance="1e-12">6.11777e-16</value>
            <value variable="w" tolerance="1e-12">9.9365e-15</value>
	    <value variable="p" tolerance="1e-12">1.04361e-13</value>
        </metric>
    </metrics>
</test>

