<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0882059</value>
            <value variable="rhou" tolerance="1e-12">62.096</value>
            <value variable="rhov" tolerance="1e-8">0.193418</value>
            <value variable="E" tolerance="1e-12">4962.78</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0773241</value>
            <value variable="rhou" tolerance="1e-12">56.0591</value>
            <value variable="rhov" tolerance="2e-6">0.313404</value>
            <value variable="E" tolerance="1e-12">4424.91</value>
        </metric>
    </metrics>
</test>


