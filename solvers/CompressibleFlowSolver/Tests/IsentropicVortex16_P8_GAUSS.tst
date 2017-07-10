<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8 GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P8_GAUSS.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P8_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.34785e-07</value>
            <value variable="rhou" tolerance="1e-12">3.43651e-07</value>
            <value variable="rhov" tolerance="1e-12">3.63171e-07</value>
            <value variable="E" tolerance="1e-12">1.2163e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">3.6015e-07</value>
            <value variable="rhou" tolerance="1e-12">1.29106e-06</value>
            <value variable="rhov" tolerance="1e-12">1.3375e-06</value>
            <value variable="E" tolerance="1e-12">4.99373e-06</value>
        </metric>
    </metrics>
</test>


