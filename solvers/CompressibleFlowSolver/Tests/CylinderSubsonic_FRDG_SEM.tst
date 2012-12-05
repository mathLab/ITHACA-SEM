<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder, Dirichlet bcs, FRDG, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_FRDG_SEM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_FRDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">10.7729</value>
            <value variable="rhou" tolerance="1e-12">1.07719</value>
            <value variable="rhov" tolerance="1e-12">0.00234645</value>
            <value variable="E" tolerance="1e-12">1.23115</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.28604</value>
            <value variable="rhou" tolerance="1e-12">0.124984</value>
            <value variable="rhov" tolerance="1e-12">0.0090395</value>
            <value variable="E" tolerance="1e-12">0.149695</value>
        </metric>
    </metrics>
</test>


