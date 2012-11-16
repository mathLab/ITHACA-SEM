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
            <value variable="rho" tolerance="1e-12">0.137274</value>
            <value variable="rhou" tolerance="1e-12">0.222511</value>
            <value variable="rhov" tolerance="1e-12">0.303939</value>
            <value variable="E" tolerance="1e-12">0.39357</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.130329</value>
            <value variable="rhou" tolerance="1e-12">0.153607</value>
            <value variable="rhov" tolerance="1e-12">0.229933</value>
            <value variable="E" tolerance="1e-12">0.431074</value>
        </metric>
    </metrics>
</test>


