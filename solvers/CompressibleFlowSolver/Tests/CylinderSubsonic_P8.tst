<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder P=8</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_P8.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_P8.xml</file>
        <file description="Restart File">CylinderSubsonic_P8.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">42.7337</value>
            <value variable="rhou" tolerance="1e-12">4.27378</value>
            <value variable="rhov" tolerance="1e-12">0.151735</value>
            <value variable="E" tolerance="1e-12">5.22544</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.39917</value>
            <value variable="rhou" tolerance="1e-12">0.198933</value>
            <value variable="rhov" tolerance="1e-12">0.114663</value>
            <value variable="E" tolerance="1e-12">0.175697</value>
        </metric>
    </metrics>
</test>


