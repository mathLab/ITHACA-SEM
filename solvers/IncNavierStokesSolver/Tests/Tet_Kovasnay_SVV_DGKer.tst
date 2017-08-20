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
            <value variable="u" tolerance="1e-6">0.00970459</value>
            <value variable="v" tolerance="1e-6">0.001718</value>
            <value variable="w" tolerance="1e-6">0.00120868</value>
            <value variable="p" tolerance="1e-6">0.00977338</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0770578</value>
            <value variable="v" tolerance="1e-6">0.00803026</value>
            <value variable="w" tolerance="1e-6">0.00678364</value>
            <value variable="p" tolerance="1e-6">0.0679505</value>
        </metric>
    </metrics>
</test>
