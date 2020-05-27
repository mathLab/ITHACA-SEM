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
            <value variable="h" tolerance="1e-8">1.00464</value>
            <value variable="hu" tolerance="1e-8">0.0216787</value>
            <value variable="hv" tolerance="1e-8">0.00551759</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-8">1.15925</value>
            <value variable="hu" tolerance="1e-8">0.29531</value>
            <value variable="hv" tolerance="1e-8">0.0483422</value>
        </metric>
    </metrics>
</test>


