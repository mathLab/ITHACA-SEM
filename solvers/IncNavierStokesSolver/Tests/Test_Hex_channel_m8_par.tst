<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hex elements, par(2), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_Hex_channel_m8_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Test_Hex_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">8.0442e-16</value>
            <value variable="v" tolerance="1e-12">6.74106e-16</value>
            <value variable="w" tolerance="1e-12">2.86244e-15</value>
	    <value variable="p" tolerance="1e-12">6.48523e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">5.23153e-15</value>
            <value variable="v" tolerance="1e-12">3.86119e-15</value>
            <value variable="w" tolerance="1e-12">4.83391e-13</value>
	    <value variable="p" tolerance="1e-12">6.32827e-13</value>
        </metric>
    </metrics>
</test>
