<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_adiabatic.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_adiabatic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.0927202</value>
            <value variable="rhou" tolerance="1e-8">64.2723</value>
            <value variable="rhov" tolerance="1e-8">0.215845</value>
            <value variable="E" tolerance="1e-8">4952.65</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.0727326</value>
            <value variable="rhou" tolerance="1e-8">70.1616</value>
            <value variable="rhov" tolerance="2e-8">0.365804</value>
            <value variable="E" tolerance="1e-8">4505.77</value>
        </metric>
    </metrics>
</test>


