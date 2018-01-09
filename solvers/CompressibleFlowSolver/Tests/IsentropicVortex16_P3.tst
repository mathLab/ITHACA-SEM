<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P3.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00180046</value>
            <value variable="rhou" tolerance="1e-12">0.0032838</value>
            <value variable="rhov" tolerance="1e-12">0.00512593</value>
            <value variable="E" tolerance="1e-12">0.0114334</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00328633</value>
            <value variable="rhou" tolerance="1e-12">0.00639972</value>
            <value variable="rhov" tolerance="1e-12">0.0118186</value>
            <value variable="E" tolerance="1e-12">0.0283689</value>
        </metric>
    </metrics>
</test>


