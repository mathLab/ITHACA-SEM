<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000397131</value>
            <value variable="rhou" tolerance="1e-12">48.1295</value>
            <value variable="rhov" tolerance="1e-12">0.145835</value>
            <value variable="rhow" tolerance="1e-8">9.33488e-06</value>
            <value variable="E" tolerance="1e-12">17519.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00139966</value>
            <value variable="rhou" tolerance="1e-12">83.3517</value>
            <value variable="rhov" tolerance="1e-12">0.50519</value>
            <value variable="rhow" tolerance="1e-8">3.11565e-05</value>
            <value variable="E" tolerance="1e-12">18953</value>
        </metric>
    </metrics>
</test>


