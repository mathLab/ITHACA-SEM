<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FrDG advection and LFRDG diffusion, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LFRDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LFRDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0885287</value>
            <value variable="rhou" tolerance="1e-12">60.3493</value>
            <value variable="rhov" tolerance="1e-8">0.21423</value>
            <value variable="E" tolerance="1e-12">4868.18</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0745423</value>
            <value variable="rhou" tolerance="1e-12">61.761</value>
            <value variable="rhov" tolerance="2e-6">0.23867</value>
            <value variable="E" tolerance="1e-12">4361.86</value>
        </metric>
    </metrics>
</test>


