<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, variable viscosity</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">5.50797</value>
            <value variable="rhou" tolerance="1e-6">2016.24</value>
            <value variable="rhov" tolerance="1e-4">25.8675</value>
            <value variable="E" tolerance="1e-1">6.27025e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">0.530811</value>
            <value variable="rhou" tolerance="1e-4">91.9676</value>
            <value variable="rhov" tolerance="1e-4">58.1355</value>
            <value variable="E" tolerance="1e-1">346847</value>
        </metric>
    </metrics>
</test>


