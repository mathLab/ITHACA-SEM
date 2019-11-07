<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, variable viscosity, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="2e-6">5.50797</value>
            <value variable="rhou" tolerance="2e-6">2016.24</value>
            <value variable="rhov" tolerance="2e-4">25.8675</value>
            <value variable="E" tolerance="2e-6">6.27025e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="2e-6">0.530811</value>
            <value variable="rhou" tolerance="2e-4">91.9676</value>
            <value variable="rhov" tolerance="2e-4">58.1354</value>
            <value variable="E" tolerance="2e-1">346847</value>
        </metric>
    </metrics>
</test>


