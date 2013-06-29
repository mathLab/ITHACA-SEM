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
            <value variable="u" tolerance="1e-6">0.00065488</value>
            <value variable="v" tolerance="1e-6">0.00017177</value>
            <value variable="w" tolerance="1e-6">0.000155218</value>
            <value variable="p" tolerance="1e-6">0.00065224</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0022794</value>
            <value variable="v" tolerance="1e-6">0.000758697</value>
            <value variable="w" tolerance="1e-6">0.000551977</value>
            <value variable="p" tolerance="1e-6">0.00887438</value>
        </metric>
    </metrics>
</test>
