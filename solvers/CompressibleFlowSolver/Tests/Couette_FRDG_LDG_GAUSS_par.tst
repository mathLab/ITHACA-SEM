<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LDG diffusion, SEM, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LDG_GAUSS_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Couette_FRDG_LDG_GAUSS_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0885112</value>
            <value variable="rhou" tolerance="1e-12">60.3078</value>
            <value variable="rhov" tolerance="1e-8">0.212138</value>
            <value variable="E" tolerance="1e-12">4868.44</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0738533</value>
            <value variable="rhou" tolerance="1e-12">60.9698</value>
            <value variable="rhov" tolerance="2e-6">0.226637</value>
            <value variable="E" tolerance="1e-12">4360.52</value>
        </metric>
    </metrics>
</test>


