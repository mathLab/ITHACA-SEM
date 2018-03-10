<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=1</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P1.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0844975</value>
            <value variable="rhou" tolerance="1e-12">0.168375</value>
            <value variable="rhov" tolerance="1e-12">0.160536</value>
            <value variable="E" tolerance="1e-12">0.4432</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0814545</value>
            <value variable="rhou" tolerance="1e-12">0.105266</value>
            <value variable="rhov" tolerance="1e-12">0.127603</value>
            <value variable="E" tolerance="1e-12">0.379398</value>
        </metric>
    </metrics>
</test>

