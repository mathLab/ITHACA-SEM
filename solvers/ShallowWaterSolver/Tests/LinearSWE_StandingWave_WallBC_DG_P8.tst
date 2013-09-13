<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Standing Wave, DG, P=8</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>LinearSWE_StandingWave_WallBC_DG_P8.xml</parameters>
    <files>
        <file description="Session File">LinearSWE_StandingWave_WallBC_DG_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-12">1.66047e-11</value>
            <value variable="u" tolerance="1e-12">1.59717e-09</value>
            <value variable="v" tolerance="1e-12">1.59717e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-12">5.0495e-11</value>
            <value variable="u" tolerance="1e-12">2.27038e-09</value>
            <value variable="v" tolerance="1e-12">2.27038e-09</value>
        </metric>
    </metrics>
</test>


