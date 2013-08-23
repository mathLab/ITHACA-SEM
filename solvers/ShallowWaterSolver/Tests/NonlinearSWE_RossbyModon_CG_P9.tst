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
            <value variable="hu" tolerance="1e-12">0.0218119</value>
            <value variable="hv" tolerance="1e-12">0.00563029</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-12">1.15544</value>
            <value variable="hu" tolerance="1e-12">0.282566</value>
            <value variable="hv" tolerance="1e-12">0.0480312</value>
        </metric>
    </metrics>
</test>


