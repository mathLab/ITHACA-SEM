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
            <value variable="u" tolerance="1e-6">0.00476409</value>
            <value variable="v" tolerance="1e-6">0.000908482</value>
            <value variable="w" tolerance="1e-6">0.000616045</value>
            <value variable="p" tolerance="1e-6">0.00422569</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0418501</value>
            <value variable="v" tolerance="1e-6">0.00387299</value>
            <value variable="w" tolerance="1e-6">0.00349122</value>
            <value variable="p" tolerance="1e-6">0.0254931</value>
        </metric>
    </metrics>
</test>
