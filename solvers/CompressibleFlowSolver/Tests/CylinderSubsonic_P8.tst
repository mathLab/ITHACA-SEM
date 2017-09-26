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
            <value variable="rho" tolerance="1e-12">42.7136</value>
            <value variable="rhou" tolerance="1e-12">4.2717</value>
            <value variable="rhov" tolerance="1e-12">0.168074</value>
            <value variable="E" tolerance="1e-12">5.22247</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.51366</value>
            <value variable="rhou" tolerance="1e-12">0.205826</value>
            <value variable="rhov" tolerance="1e-12">0.156464</value>
            <value variable="E" tolerance="1e-12">0.192841</value>
        </metric>
    </metrics>
</test>


