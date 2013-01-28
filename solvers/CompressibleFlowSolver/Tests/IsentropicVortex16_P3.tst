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
            <value variable="rho" tolerance="1e-12">0.00266419</value>
            <value variable="rhou" tolerance="1e-12">0.00361032</value>
            <value variable="rhov" tolerance="1e-12">0.181602</value>
            <value variable="E" tolerance="1e-12">0.102132</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00555232</value>
            <value variable="rhou" tolerance="1e-12">0.00814604</value>
            <value variable="rhov" tolerance="1e-12">0.113819</value>
            <value variable="E" tolerance="1e-12">0.12761</value>
        </metric>
    </metrics>
</test>


