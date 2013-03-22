<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=11, periodic bcs, FRDG, GLL_LAGRANGE_SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex_FRDG_SEM.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex_FRDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0447202</value>
            <value variable="rhou" tolerance="1e-12">0.0789786</value>
            <value variable="rhov" tolerance="1e-12">0.0618054</value>
            <value variable="E" tolerance="1e-12">0.175818</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0454094</value>
            <value variable="rhou" tolerance="1e-12">0.0898222</value>
            <value variable="rhov" tolerance="1e-12">0.0728035</value>
            <value variable="E" tolerance="1e-12">0.180614</value>
        </metric>
    </metrics>
</test>


