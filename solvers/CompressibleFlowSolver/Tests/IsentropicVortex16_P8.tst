<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P8.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0351889</value>
            <value variable="rhou" tolerance="1e-12">0.0569376</value>
            <value variable="rhov" tolerance="1e-12">0.0775479</value>
            <value variable="E" tolerance="1e-12">0.100569</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0315969</value>
            <value variable="rhou" tolerance="1e-12">0.0362279</value>
            <value variable="rhov" tolerance="1e-12">0.0547527</value>
            <value variable="E" tolerance="1e-12">0.0964763</value>
        </metric>
    </metrics>
</test>


