<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0878644</value>
            <value variable="rhou" tolerance="1e-12">60.3309</value>
            <value variable="rhov" tolerance="1e-8">0.230213</value>
            <value variable="E" tolerance="1e-12">4926.55</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0743376</value>
            <value variable="rhou" tolerance="1e-12">60.8075</value>
            <value variable="rhov" tolerance="1e-8">0.302839</value>
            <value variable="E" tolerance="1e-12">4448.92</value>
        </metric>
    </metrics>
</test>


