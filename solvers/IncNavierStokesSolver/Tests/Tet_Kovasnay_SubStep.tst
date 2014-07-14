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
            <value variable="u" tolerance="1e-6">0.00155108</value>
            <value variable="v" tolerance="1e-6">0.00142855</value>
            <value variable="w" tolerance="1e-6">0.000589637</value>
            <value variable="p" tolerance="1e-6">0.00658042</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00721426</value>
            <value variable="v" tolerance="1e-6">0.00609401</value>
            <value variable="w" tolerance="1e-6">0.00684261</value>
            <value variable="p" tolerance="1e-6">0.0456315</value>
        </metric>
    </metrics>
</test>
