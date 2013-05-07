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
            <value variable="u" tolerance="1e-6">0.000798217</value>
            <value variable="v" tolerance="1e-6">0.000342217</value>
            <value variable="w" tolerance="1e-6">0.00018671</value>
            <value variable="p" tolerance="1e-6">0.000872105</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00287857</value>
            <value variable="v" tolerance="1e-6">0.00152883</value>
            <value variable="w" tolerance="1e-6">0.000611055</value>
            <value variable="p" tolerance="1e-6">0.00999607</value>
        </metric>
    </metrics>
</test>
