<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=11, periodic bcs, FRHU, GLL_LAGRANGE_SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex_FRHU_SEM.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex_FRHU_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0547039</value>
            <value variable="rhou" tolerance="1e-12">0.0795769</value>
            <value variable="rhov" tolerance="1e-12">1.09364</value>
            <value variable="E" tolerance="1e-12">0.751095</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.144182</value>
            <value variable="rhou" tolerance="1e-12">0.155464</value>
            <value variable="rhov" tolerance="1e-12">0.689679</value>
            <value variable="E" tolerance="1e-12">0.898678</value>
        </metric>
    </metrics>
</test>


