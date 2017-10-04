<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=3 GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P3_GAUSS.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P3_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00395533</value>
            <value variable="rhou" tolerance="1e-12">0.00784074</value>
            <value variable="rhov" tolerance="1e-12">0.00735103</value>
            <value variable="E" tolerance="1e-12">0.0219855</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00451892</value>
            <value variable="rhou" tolerance="1e-12">0.0103558</value>
            <value variable="rhov" tolerance="1e-12">0.0104917</value>
            <value variable="E" tolerance="1e-12">0.033977</value>
        </metric>
    </metrics>
</test>


