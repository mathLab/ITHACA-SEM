<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hex elements, par(2), varP</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=XxtMultiLevelStaticCond Hex_channel_varP.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Hex_channel_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-11">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-11">0</value>
        </metric>
    </metrics>
</test>
