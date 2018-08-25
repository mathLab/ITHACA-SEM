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
            <value variable="u" tolerance="1e-6">0.00475271</value>
            <value variable="v" tolerance="1e-6">0.000915189</value>
            <value variable="w" tolerance="1e-6">0.000609439</value>
            <value variable="p" tolerance="1e-6">0.00391764</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.041269</value>
            <value variable="v" tolerance="1e-6">0.00385715</value>
            <value variable="w" tolerance="1e-6">0.00355733</value>
            <value variable="p" tolerance="1e-6">0.0243688</value>
        </metric>
    </metrics>
</test>
