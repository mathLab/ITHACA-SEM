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
            <value variable="rho" tolerance="1e-12">0.0543871</value>
            <value variable="rhou" tolerance="1e-12">0.0794009</value>
            <value variable="rhov" tolerance="1e-12">1.09357</value>
            <value variable="E" tolerance="1e-12">0.750935</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.137667</value>
            <value variable="rhou" tolerance="1e-12">0.146687</value>
            <value variable="rhov" tolerance="1e-12">0.685895</value>
            <value variable="E" tolerance="1e-12">0.898705</value>
        </metric>
    </metrics>
</test>


