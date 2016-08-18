<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hexahedral elements, P=3, flowrate</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=XxtStaticCond Hex_channel_m3_Flowrate.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Hex_channel_m3_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">1.13277e-15</value>
            <value variable="p" tolerance="1e-8">4.56093e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.48905e-16</value>
            <value variable="v" tolerance="1e-12">3.12187e-16</value>
            <value variable="w" tolerance="1e-12">1.90958e-14</value>
            <value variable="p" tolerance="1e-8">5.11546e-14</value>
        </metric>
    </metrics>
</test>
