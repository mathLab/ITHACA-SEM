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
            <value variable="rho" tolerance="1e-12">1.83107e-08</value>
            <value variable="rhou" tolerance="1e-12">4.88852e-08</value>
            <value variable="rhov" tolerance="1e-12">7.90372e-08</value>
            <value variable="E" tolerance="1e-12">2.72206e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">7.2128e-08</value>
            <value variable="rhou" tolerance="1e-12">1.63335e-07</value>
            <value variable="rhov" tolerance="1e-12">3.42155e-07</value>
            <value variable="E" tolerance="1e-12">8.66884e-07</value>
        </metric>
    </metrics>
</test>


