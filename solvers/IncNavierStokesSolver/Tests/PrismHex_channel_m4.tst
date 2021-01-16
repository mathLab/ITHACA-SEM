<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral and Prismatic elements, P=4</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PrismHex_channel_m4.xml</parameters>
    <files>
        <file description="Session File">PrismHex_channel_m4.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.79242e-15</value>
            <value variable="v" tolerance="1e-12">1.76764e-15</value>
            <value variable="w" tolerance="1e-12">7.97721e-15</value>
	    <value variable="p" tolerance="1e-12">1.11992e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">9.44016e-15</value>
            <value variable="v" tolerance="1e-12">9.24349e-15</value>
            <value variable="w" tolerance="1e-12">4.39926e-14</value>
	    <value variable="p" tolerance="1e-12">6.00853e-13</value>
        </metric>
    </metrics>
</test>

