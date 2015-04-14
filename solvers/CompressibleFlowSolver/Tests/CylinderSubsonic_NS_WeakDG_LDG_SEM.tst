<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_SEM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">5.50797</value>
            <value variable="rhou" tolerance="1e-12">2016.24</value>
            <value variable="rhov" tolerance="1e-8">25.872</value>
            <value variable="E" tolerance="1e-12">6.27025e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.53102</value>
            <value variable="rhou" tolerance="1e-12">92.5537</value>
            <value variable="rhov" tolerance="1e-8">58.4422</value>
            <value variable="E" tolerance="1e-12">347067</value>
        </metric>
    </metrics>
</test>


