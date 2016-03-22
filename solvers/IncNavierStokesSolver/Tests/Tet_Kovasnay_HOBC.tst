<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using High Order Outflow BCsd</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_HOBC.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_HOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00116189</value>
            <value variable="v" tolerance="1e-6">0.000235289</value>
            <value variable="w" tolerance="1e-6">0.000191207</value>
            <value variable="p" tolerance="1e-6">0.00112908</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00929358</value>
            <value variable="v" tolerance="1e-6">0.0014178</value>
            <value variable="w" tolerance="1e-6">0.000961437</value>
            <value variable="p" tolerance="1e-6">0.00692496</value>
        </metric>
    </metrics>
</test>
