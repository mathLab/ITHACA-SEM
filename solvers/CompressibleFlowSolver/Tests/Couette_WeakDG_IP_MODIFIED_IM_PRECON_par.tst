<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">3.40454e-04</value>
            <value variable="rhou" tolerance="1e-6">1.57601e-01</value>
            <value variable="rhov" tolerance="1e-7">3.47048e-02</value>
            <value variable="E" tolerance="1e-4">3.24624e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-9">2.78071e-04</value>
            <value variable="rhou" tolerance="1e-6">1.32527e-01</value>
            <value variable="rhov" tolerance="1e-7">4.35076e-02</value>
            <value variable="E" tolerance="1e-4">3.37727e+01</value>
        </metric>

    </metrics>
</test>
