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
            <value variable="rho" tolerance="1e-12">42.8519</value>
            <value variable="rhou" tolerance="1e-12">4.28839</value>
            <value variable="rhov" tolerance="1e-12">0.0809356</value>
            <value variable="E" tolerance="1e-12">5.24313</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.3205</value>
            <value variable="rhou" tolerance="1e-12">0.198933</value>
            <value variable="rhov" tolerance="1e-12">0.114663</value>
            <value variable="E" tolerance="1e-12">0.159681</value>
        </metric>
    </metrics>
</test>


