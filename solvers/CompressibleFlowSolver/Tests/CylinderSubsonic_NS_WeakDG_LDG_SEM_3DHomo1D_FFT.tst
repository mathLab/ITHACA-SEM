<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.89459</value>
            <value variable="rhou" tolerance="1e-12">1426.6</value>
            <value variable="rhov" tolerance="1e-12">58.1957</value>
            <value variable="rhow" tolerance="1e-12">17.3092</value>
            <value variable="E" tolerance="1e-12">4.43372e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.343891</value>
            <value variable="rhou" tolerance="1e-12">165.462</value>
            <value variable="rhov" tolerance="1e-12">84.0672</value>
            <value variable="rhow" tolerance="1e-12">1.00005</value>
            <value variable="E" tolerance="1e-12">291493</value>
        </metric>
    </metrics>
</test>


