<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rossby modon, CG, P=9</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_RossbyModon_CG_P9.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_RossbyModon_CG_P9.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h" tolerance="1e-12">1.00464</value>
            <value variable="hu" tolerance="1e-12">0.0216812</value>
            <value variable="hv" tolerance="1e-12">0.00553698</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-12">1.15917</value>
            <value variable="hu" tolerance="1e-12">0.295374</value>
            <value variable="hv" tolerance="1e-12">0.0483798</value>
        </metric>
    </metrics>
</test>


