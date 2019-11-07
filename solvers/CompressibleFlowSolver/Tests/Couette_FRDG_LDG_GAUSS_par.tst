<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LDG diffusion, SEM, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Couette_FRDG_LDG_GAUSS_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Couette_FRDG_LDG_GAUSS_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-7">0.0885179</value>
            <value variable="rhou" tolerance="1e-4">60.3076</value>
            <value variable="rhov" tolerance="1e-6">0.211622</value>
            <value variable="E" tolerance="1e-2">4866.11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0738853</value>
            <value variable="rhou" tolerance="1e-4">60.9894</value>
            <value variable="rhov" tolerance="2e-6">0.226922</value>
            <value variable="E" tolerance="1e-2">4357.43</value>
        </metric>
    </metrics>
</test>


