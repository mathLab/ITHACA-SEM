<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LFRHU diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LFRHU_SEM.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LFRHU_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0889271</value>
            <value variable="rhou" tolerance="1e-12">62.1128</value>
            <value variable="rhov" tolerance="1e-8">0.175956</value>
            <value variable="E" tolerance="1e-12">4905.04</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0760154</value>
            <value variable="rhou" tolerance="1e-12">56.0464</value>
            <value variable="rhov" tolerance="2e-6">0.265763</value>
            <value variable="E" tolerance="1e-12">4381.12</value>
        </metric>
    </metrics>
</test>


