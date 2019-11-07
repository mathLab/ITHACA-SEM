<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRHU_LDG_SEM.xml</parameters>
    <files>
        <file description="Session File">Couette_FRHU_LDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="2e-7">0.0889264</value>
            <value variable="rhou" tolerance="1e-12">62.0963</value>
            <value variable="rhov" tolerance="1e-8">0.175956</value>
            <value variable="E" tolerance="1e-12">4903.95</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0760421</value>
            <value variable="rhou" tolerance="1e-12">56.0835</value>
            <value variable="rhov" tolerance="2e-6">0.265899</value>
            <value variable="E" tolerance="1e-12">4333.92</value>
        </metric>
    </metrics>
</test>


