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
            <value variable="rho" tolerance="1e-12">0.00271826</value>
            <value variable="rhou" tolerance="1e-12">0.00495649</value>
            <value variable="rhov" tolerance="1e-12">0.00780554</value>
            <value variable="E" tolerance="1e-12">0.0173238</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00554565</value>
            <value variable="rhou" tolerance="1e-12">0.00971562</value>
            <value variable="rhov" tolerance="1e-12">0.0180984</value>
            <value variable="E" tolerance="1e-12">0.0408197</value>
        </metric>
    </metrics>
</test>


