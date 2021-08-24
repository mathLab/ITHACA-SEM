<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder, Dirichlet bcs, FRHU, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_FRHU_SEM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_FRHU_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">10.7744</value>
            <value variable="rhou" tolerance="1e-12">1.07625</value>
            <value variable="rhov" tolerance="1e-8">0.0228178</value>
            <value variable="E" tolerance="1e-12">2.228e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.22537</value>
            <value variable="rhou" tolerance="1e-12">0.148574</value>
            <value variable="rhov" tolerance="1e-8">0.0874203</value>
            <value variable="E" tolerance="1e-12">253419</value>
        </metric>
    </metrics>
</test>


