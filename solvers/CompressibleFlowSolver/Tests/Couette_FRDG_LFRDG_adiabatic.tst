<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRDG advection and LFRDG diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LFRDG_adiabatic.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LFRDG_adiabatic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.0926743</value>
            <value variable="rhou" tolerance="1e-8">64.2553</value>
            <value variable="rhov" tolerance="1e-8">0.210065</value>
            <value variable="E" tolerance="1e-8">4937.29</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.0753016</value>
            <value variable="rhou" tolerance="1e-8">70.0111</value>
            <value variable="rhov" tolerance="2e-8">0.358313</value>
            <value variable="E" tolerance="1e-8">4483.59</value>
        </metric>
    </metrics>
</test>


