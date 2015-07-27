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
            <value variable="rho" tolerance="1e-8">0.0927225</value>
            <value variable="rhou" tolerance="1e-8">64.2725</value>
            <value variable="rhov" tolerance="1e-8">0.215208</value>
            <value variable="E" tolerance="1e-8">4953.71</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.0729093</value>
            <value variable="rhou" tolerance="1e-8">70.1654</value>
            <value variable="rhov" tolerance="2e-8">0.366318</value>
            <value variable="E" tolerance="1e-8">4507.27</value>
        </metric>
    </metrics>
</test>


