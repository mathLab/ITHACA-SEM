<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LFRHU diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRHU_LFRHU_SEM.xml</parameters>
    <files>
        <file description="Session File">Couette_FRHU_LFRHU_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.0889266</value>
            <value variable="rhou" tolerance="1e-12">62.1102</value>
            <value variable="rhov" tolerance="1e-8">0.175939</value>
            <value variable="E" tolerance="1e-12">4904.77</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0760033</value>
            <value variable="rhou" tolerance="1e-12">56.0587</value>
            <value variable="rhov" tolerance="2e-6">0.26569</value>
            <value variable="E" tolerance="1e-12">4376.14</value>
        </metric>
    </metrics>
</test>


