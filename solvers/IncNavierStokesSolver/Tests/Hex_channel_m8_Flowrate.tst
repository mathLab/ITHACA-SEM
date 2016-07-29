<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hexahedral elements, P=8, flowrate</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_channel_m8_Flowrate.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Hex_channel_m8_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.49933e-12</value>
            <value variable="v" tolerance="1e-12">1.53987e-12</value>
            <value variable="w" tolerance="1e-12">1.09127e-11</value>
            <value variable="p" tolerance="1e-8">1.2548e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.38902e-12</value>
            <value variable="v" tolerance="1e-12">7.15068e-12</value>
            <value variable="w" tolerance="1e-12">8.04076e-11</value>
            <value variable="p" tolerance="1e-8">2.87263e-09</value>
        </metric>
    </metrics>
</test>
