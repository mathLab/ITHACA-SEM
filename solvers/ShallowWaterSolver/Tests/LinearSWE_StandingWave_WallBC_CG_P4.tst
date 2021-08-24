<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Standing Wave, CG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>LinearSWE_StandingWave_WallBC_CG_P4.xml</parameters>
    <files>
        <file description="Session File">LinearSWE_StandingWave_WallBC_CG_P4.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-12">4.34843e-05</value>
            <value variable="u" tolerance="1e-12">2.56707e-05</value>
            <value variable="v" tolerance="1e-12">2.56707e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-12">0.000207527</value>
            <value variable="u" tolerance="1e-12">9.25444e-05</value>
            <value variable="v" tolerance="1e-12">9.25444e-05</value>
        </metric>
    </metrics>
</test>


