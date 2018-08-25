<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.000791971</value>
            <value variable="v" tolerance="1e-6">0.000253692</value>
            <value variable="w" tolerance="1e-6">0.000182633</value>
            <value variable="p" tolerance="1e-6">0.000829782</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00900795</value>
            <value variable="v" tolerance="1e-6">0.00146155</value>
            <value variable="w" tolerance="1e-6">0.00111463</value>
            <value variable="p" tolerance="2e-6">0.00892642</value>
        </metric>
    </metrics>
</test>
