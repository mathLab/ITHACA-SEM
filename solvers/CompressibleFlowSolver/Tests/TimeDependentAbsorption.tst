<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Time dependent absorption filter</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>TimeDependentAbsorption.xml</parameters>
    <files>
        <file description="Session File">TimeDependentAbsorption.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00969206</value>
            <value variable="rhou" tolerance="1e-12">0.0315752</value>
            <value variable="rhov" tolerance="1e-12">0.030466</value>
            <value variable="E" tolerance="1e-12">0.105786</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.000974088</value>
            <value variable="rhou" tolerance="1e-12">0.00282339</value>
            <value variable="rhov" tolerance="1e-12">0.00270836</value>
            <value variable="E" tolerance="1e-12">0.0105478</value>
        </metric>
    </metrics>
</test>

