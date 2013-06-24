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
            <value variable="rho" tolerance="1e-12">3.89561</value>
            <value variable="rhou" tolerance="1e-12">1433.1</value>
            <value variable="rhov" tolerance="1e-12">167.108</value>
            <value variable="rhow" tolerance="1e-12">17.3092</value>
            <value variable="E" tolerance="1e-12">4.43381e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.905233</value>
            <value variable="rhou" tolerance="1e-12">364.429</value>
            <value variable="rhov" tolerance="1e-12">280.234</value>
            <value variable="rhow" tolerance="1e-12">1.00016</value>
            <value variable="E" tolerance="1e-12">467054</value>
        </metric>
    </metrics>
</test>


