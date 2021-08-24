<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, Variable Viscosity</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM_VariableMu.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_VariableMu.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0805614</value>
            <value variable="rhou" tolerance="1e-12">51.8867</value>
            <value variable="rhov" tolerance="1e-8">0.222213</value>
            <value variable="E" tolerance="1e-12">4415.64</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0716895</value>
            <value variable="rhou" tolerance="1e-12">55.3017</value>
            <value variable="rhov" tolerance="1e-8">0.424398</value>
            <value variable="E" tolerance="1e-12">4227.23</value>
        </metric>
    </metrics>
</test>


