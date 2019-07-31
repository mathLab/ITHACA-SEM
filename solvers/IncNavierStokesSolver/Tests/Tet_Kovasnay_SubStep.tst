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
            <value variable="u" tolerance="1e-6">0.00105429</value>
            <value variable="v" tolerance="1e-6">0.000275093</value>
            <value variable="w" tolerance="1e-6">0.000218578</value>
            <value variable="p" tolerance="1e-6">0.001029114</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00900795</value>
            <value variable="v" tolerance="1e-6">0.0014178</value>
            <value variable="w" tolerance="1e-6">0.0010477</value>
            <value variable="p" tolerance="2e-6">0.00540266</value>
        </metric>
    </metrics>
</test>
