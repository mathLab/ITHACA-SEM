<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Helmholtz 3D with hybrid element, mixed order and low energy precondition with linear space </description>
    <executable>Helmholtz3D</executable>
    <parameters>-v -I GlobalSysSoln=IterativeStaticCond -I Preconditioner=FullLinearSpaceWithLowEnergyBlock CubeAllElements.xml</parameters>
    <files>
        <file description="Session File">CubeAllElements.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-6">0.000833058</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.024387</value>
        </metric>
        <metric type="Precon" id="3">
            <value tolerance="2">46</value>
        </metric>
    </metrics>
</test>
