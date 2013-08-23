<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rossby modon, DG, P=9</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_RossbyModon_DG_P9.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_RossbyModon_DG_P9.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h" tolerance="1e-12">1.00464</value>
            <value variable="hu" tolerance="1e-12">0.0218085</value>
            <value variable="hv" tolerance="1e-12">0.00563633</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-12">1.15518</value>
            <value variable="hu" tolerance="1e-12">0.322833</value>
            <value variable="hv" tolerance="1e-12">0.0491446</value>
        </metric>
    </metrics>
</test>


