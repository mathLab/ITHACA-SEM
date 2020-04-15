<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Pyramidic elements, variable P</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Pyr_channel_varP.xml</parameters>
    <files>
        <file description="Session File">Pyr_channel_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
        <value variable="p" tolerance="7e-12">0</value>
        </metric>
    </metrics>
</test>
