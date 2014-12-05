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
            <value variable="u" tolerance="1e-6">0.00064288</value>
            <value variable="v" tolerance="1e-6">0.000157378</value>
            <value variable="w" tolerance="1e-6">0.000166338</value>
            <value variable="p" tolerance="1e-6">0.000563287</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0022434</value>
            <value variable="v" tolerance="1e-6">0.000779688</value>
            <value variable="w" tolerance="1e-6">0.000658035</value>
            <value variable="p" tolerance="1e-6">0.00317511</value>
        </metric>
    </metrics>
</test>
