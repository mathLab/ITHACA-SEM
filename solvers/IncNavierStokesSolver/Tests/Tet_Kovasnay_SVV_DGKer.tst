<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using DG SVV Kerneal and dealiasing</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_SVV_DGKer.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_SVV_DGKer.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00476290</value>
            <value variable="v" tolerance="1e-6">0.000917981</value>
            <value variable="w" tolerance="1e-6">0.000637714</value>
            <value variable="p" tolerance="1e-6">0.00409401</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0410753</value>
            <value variable="v" tolerance="1e-6">0.00407471</value>
            <value variable="w" tolerance="1e-6">0.00346771</value>
            <value variable="p" tolerance="1e-6">0.02489</value>
        </metric>
    </metrics>
</test>
