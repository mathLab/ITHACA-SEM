<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LDG diffusion, SEM, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-metis Couette_FRDG_LDG_GAUSS_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Couette_FRDG_LDG_GAUSS_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0885176</value>
            <value variable="rhou" tolerance="1e-12">60.3077</value>
            <value variable="rhov" tolerance="1e-8">0.2117</value>
            <value variable="E" tolerance="1e-12">4867.55</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0740721</value>
            <value variable="rhou" tolerance="1e-12">60.9805</value>
            <value variable="rhov" tolerance="2e-6">0.225975</value>
            <value variable="E" tolerance="1e-12">4357.9</value>
        </metric>
    </metrics>
</test>


