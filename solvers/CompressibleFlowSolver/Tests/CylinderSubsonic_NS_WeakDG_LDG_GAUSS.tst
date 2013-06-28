<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">5.50797</value>
            <value variable="rhou" tolerance="1e-12">2016.24</value>
            <value variable="rhov" tolerance="1e-8">25.87</value>
            <value variable="E" tolerance="1e-12">6.27025e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.459344</value>
            <value variable="rhou" tolerance="1e-12">89.6381</value>
            <value variable="rhov" tolerance="1e-8">43.8967</value>
            <value variable="E" tolerance="1e-12">324671</value>
        </metric>
    </metrics>
</test>


