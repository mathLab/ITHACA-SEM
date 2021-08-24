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
            <value variable="u" tolerance="1e-6">0.00116277</value>
            <value variable="v" tolerance="1e-6">0.000311973</value>
            <value variable="w" tolerance="1e-6">0.000231987</value>
            <value variable="p" tolerance="1e-6">0.00119195</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00900795</value>
            <value variable="v" tolerance="1e-6">0.0014178</value>
            <value variable="w" tolerance="1e-6">0.00102937</value>
            <value variable="p" tolerance="1e-6">0.00603119</value>
        </metric>
    </metrics>
</test>
