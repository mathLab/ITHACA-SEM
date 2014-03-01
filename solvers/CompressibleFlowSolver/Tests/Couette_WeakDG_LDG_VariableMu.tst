<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, Variable Viscosity</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_VariableMu.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_VariableMu.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0807273</value>
            <value variable="rhou" tolerance="1e-12">51.8835</value>
            <value variable="rhov" tolerance="1e-8">0.223588</value>
            <value variable="E" tolerance="1e-12">4415.29</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0721994</value>
            <value variable="rhou" tolerance="1e-12">55.3784</value>
            <value variable="rhov" tolerance="1e-8">0.427537</value>
            <value variable="E" tolerance="1e-12">4222.3</value>
        </metric>
    </metrics>
</test>


