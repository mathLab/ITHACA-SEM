<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=11, periodic bcs, FRSD, GLL_LAGRANGE_SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex_FRSD_SEM.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex_FRSD_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0455572</value>
            <value variable="rhou" tolerance="1e-12">0.0790375</value>
            <value variable="rhov" tolerance="1e-12">0.0714988</value>
            <value variable="E" tolerance="1e-12">0.179773</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0494984</value>
            <value variable="rhou" tolerance="1e-12">0.0881977</value>
            <value variable="rhov" tolerance="1e-12">0.129721</value>
            <value variable="E" tolerance="1e-12">0.257346</value>
        </metric>
    </metrics>
</test>


