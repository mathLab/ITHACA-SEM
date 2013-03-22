<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder Mixed</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonicMix.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonicMix.xml</file>
        <file description="Restart File">CylinderSubsonicMix.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">42.8503</value>
            <value variable="rhou" tolerance="1e-12">4.28809</value>
            <value variable="rhov" tolerance="1e-12">0.0812521</value>
            <value variable="E" tolerance="1e-12">5.24289</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.32074</value>
            <value variable="rhou" tolerance="1e-12">0.200515</value>
            <value variable="rhov" tolerance="1e-12">0.358861</value>
            <value variable="E" tolerance="1e-12">0.159723</value>
        </metric>
    </metrics>
</test>


