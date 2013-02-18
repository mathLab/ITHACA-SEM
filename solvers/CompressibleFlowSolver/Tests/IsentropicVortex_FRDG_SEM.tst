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
            <value variable="rho" tolerance="1e-12">0.0525223</value>
            <value variable="rhou" tolerance="1e-12">0.0806749</value>
            <value variable="rhov" tolerance="1e-12">1.09338</value>
            <value variable="E" tolerance="1e-12">0.749846</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0924581</value>
            <value variable="rhou" tolerance="1e-12">0.0917078</value>
            <value variable="rhov" tolerance="1e-12">0.651416</value>
            <value variable="E" tolerance="1e-12">0.881846</value>
        </metric>
    </metrics>
</test>


