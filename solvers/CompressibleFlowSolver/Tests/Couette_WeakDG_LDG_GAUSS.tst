<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0878864</value>
            <value variable="rhou" tolerance="1e-12">60.3318</value>
            <value variable="rhov" tolerance="1e-8">0.229384</value>
            <value variable="E" tolerance="1e-12">4925.53</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0743007</value>
            <value variable="rhou" tolerance="1e-12">61.7707</value>
            <value variable="rhov" tolerance="2e-6">0.278381</value>
            <value variable="E" tolerance="1e-12">4406.17</value>
        </metric>
    </metrics>
</test>


